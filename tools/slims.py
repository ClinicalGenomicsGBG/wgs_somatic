import os
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from slims.slims import Slims
from slims.criteria import is_one_of, equals, conjunction, not_equals

from tools.helpers import read_config
from definitions import WRAPPER_CONFIG_PATH, ROOT_DIR


class slims_credentials:
    def __init__(self, slims_credentials_path):
        config = read_config(slims_credentials_path)
        self.url = config['slims']['url']
        self.user = config['slims']['user']
        self.password = config['slims']['password']


config = read_config(WRAPPER_CONFIG_PATH)
slims_credentials_path = config['slims_credentials_path']
credentials = slims_credentials(slims_credentials_path)

slims_connection = Slims('wgs-somatic_query',
                         credentials.url,
                         credentials.user,
                         credentials.password)


class SlimsSample:
    def __init__(self, sample_name, run_tag=''):
        self.sample_name = sample_name
        self.run_tag = run_tag

        self._dna = None
        self._fastq = None
        self._bioinformatics = None

    @property
    def dna(self):
        if not self._dna:
            records = slims_connection.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 6)))

            if len(records) > 1:
                raise Exception('More than 1 DNA somehow.')

            if records:
                #print(records)
                self._dna = records[0]

        return self._dna

    @property
    def fastq(self):
        if not self.run_tag:
            raise Exception('Can not fetch fastq without a set run tag.')
        if not self._fastq:
            records = slims_connection.fetch('Content', conjunction()
                                  .add(equals('cntn_id', self.sample_name))
                                  .add(equals('cntn_fk_contentType', 22))
                                  .add(equals('cntn_cstm_runTag', self.run_tag)))
            if len(records) > 1:
                raise Exception('More than 1 fastq somehow.')

            if records:
                self._fastq = records[0]

        return self._fastq


def translate_slims_info(record):
    sample_name = record.cntn_id.value

    # 29 = WOPR
    # 54 = wgs_somatic

    pipeline_pks = record.cntn_cstm_secondaryAnalysis.value
    pcr = record.cntn_cstm_pcr.value

    investigator = 'CGG'  # NOTE: Needed?

    department = None
    responder_emails = []
    if record.cntn_cstm_department.value is not None:
        department_record = slims_connection.fetch_by_pk('ReferenceDataRecord', record.cntn_cstm_department.value)
        department = department_record.rdrc_name.value  
        responder_records = [slims_connection.fetch_by_pk('ReferenceDataRecord', pk) for
                            pk in department_record.rdrc_cstm_responder.value]
        responder_emails = [rec.rdrc_cstm_email.value for rec in responder_records]

    is_research = record.cntn_cstm_research.value
    research_project = record.cntn_cstm_researchProject.value

    is_priority = True if record.cntn_cstm_priority.value else False

    gender = record.gender.value

    tumorNormalType = record.cntn_cstm_tumorNormalType.value
    tumorNormalID = record.cntn_cstm_tumorNormalID.value

    tertiary_analysis = record.cntn_cstm_tertiaryAnalysis.value

    master = {
        'content_id': sample_name,
        'investigator': investigator,
        'department': department,
        'responder_mails': responder_emails,
        'is_research': is_research,
        'research_project': research_project,
        'gender': gender,
        'is_priority': is_priority,
        'pcr': pcr,
        'tumorNormalType': tumorNormalType,
        'tumorNormalID': tumorNormalID,
        'secondary_analysis': pipeline_pks,
        'tertiary_analysis': tertiary_analysis
    }
    return master


def get_sample_slims_info(Sctx, run_tag):
    """Query the slims API for relevant metadata given sample_name in samplesheet."""
    SSample = SlimsSample(Sctx.sample_name, run_tag)

    if not SSample.dna:
        Sctx.slims_info = {}
        return
    return translate_slims_info(SSample.dna)


def download_hcp_fq(remote_key, logger, hcp_runtag):
    """Find and download fqs from HCP to fastqdir on seqstore for run"""
    config = read_config(WRAPPER_CONFIG_PATH)

    # Read all download locations from the config
    download_locations = config["hcp"]["download_locations"]

    hcp_download_runpath = f'{config["hcp_download_dir"]}/{hcp_runtag}'  # Directory for downloaded files
    hcp_path = f'{hcp_download_runpath}/{os.path.basename(remote_key)}'  # Full path of the downloaded file

    if os.path.exists(hcp_path):
        logger.info(f'{os.path.basename(remote_key)} already exists in {hcp_download_runpath}')
        return hcp_path

    os.makedirs(hcp_download_runpath, exist_ok=True)

    # Iterate through the locations in the order specified in the config
    for location_name, location_details in download_locations.items():
        try:
            bucket = location_details["bucket"]  # Bucket is specified in the config

            logger.info(f"Trying to download from {location_name})")

            # Construct the download command
            qrsh = [
                "qrsh",
                "-q", config["hcp"]["queue"],
                "-N", f"hcp_download_{os.path.basename(remote_key)}",
                "-pe", "mpi", "1",
                "-now", "no",
                "-cwd", "-V"
            ]

            main_args = [
                "python", os.path.abspath(config["hcp"]["download_script"]),
                "-l", hcp_path,
                "-r", remote_key,
                "-c", location_details["credentials_file"],
                "-b", bucket
            ]

            optional_args = [
                "--connect_timeout", str(config["hcp"]["connect_timeout"]),
                "--read_timeout", str(config["hcp"]["read_timeout"]),
                "--retries", str(config["hcp"]["retries"])
            ]

            cmd = qrsh + main_args + optional_args
            logger.info(f"Running hcp_download.py with args: {cmd}")

            # Run the download command
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            if stdout:
                logger.info(stdout.decode('utf-8'))
            if stderr:
                logger.error(stderr.decode('utf-8'))

            if process.returncode == 0:
                # Verify the file exists after download
                if os.path.exists(hcp_path):
                    logger.info(f"Successfully downloaded {os.path.basename(remote_key)} from {url}")
                    return hcp_path
                else:
                    logger.error(f"Download succeeded but file not found at {hcp_path}")
                    raise RuntimeError(f"File not found after download: {hcp_path}")
            else:
                logger.warning(f"Failed to download from {location_name}")

        except Exception as e:
            logger.warning(f"Error while trying to download from {location_name}: {e}")

    # If none of the locations worked, raise an error
    logger.error(f"Failed to download {remote_key} from all locations.")
    raise RuntimeError(f"Remote key {remote_key} could not be found on any of the specified locations.")


def decompress_downloaded_fastq(complete_file_path, logger):
    config = read_config(WRAPPER_CONFIG_PATH)

    wrapper_log_path = config["wrapper_log_path"]

    filename = os.path.basename(complete_file_path) # This is the filename of the downloaded file
    standardout_decompress = os.path.join(wrapper_log_path, f"decompress_{filename}.stdout")
    standarderr_decompress = os.path.join(wrapper_log_path, f"decompress_{filename}.stderr")

    queue = config["hcp"]["queue"]
    threads = config["hcp"]["threads"]
    compression_type = filename.split('.')[-1] # This is the compression type of the downloaded file, could be either 'spring' or 'fasterq'

    if compression_type == 'spring':
        # Decompress the file using spring
        complete_decompressed_file_path = complete_file_path.replace('.spring', '.fastq.gz')
        if not os.path.exists(complete_decompressed_file_path):
            decompress_script = os.path.join(ROOT_DIR, config["hcp"]["spring_script"])

            qsub_args = ["qsub", "-N", f"decompressing_{filename}", "-q", queue, "-sync", "y", "-pe", "mpi", f"{threads}",
                        "-o", standardout_decompress, "-e", standarderr_decompress, "-v", f"THREADS={threads}",
                        decompress_script, complete_file_path, complete_decompressed_file_path, str(threads)]

            logger.info(f"Decompressing {filename} using spring with args: {qsub_args}")
            subprocess.call(qsub_args)
            logger.info(f"Done decompressing {filename}")
        return complete_decompressed_file_path

    elif compression_type == 'fasterq':
        # Decompress the file using petasuite
        complete_decompressed_file_path = complete_file_path.replace('.fasterq', '.fastq.gz')
        if not os.path.exists(complete_decompressed_file_path):
            decompress_script = os.path.join(ROOT_DIR, config["hcp"]["peta_script"])

            peta_args = ["qsub", "-N", f"decompressing_file_{filename}", "-q", queue, "-sync", "y",
                        "-pe", "mpi", f"{threads}", "-o", standardout_decompress, "-e", standarderr_decompress, "-v", f"THREADS={threads}",
                        decompress_script, complete_file_path, str(threads)] 

            logger.info(f"Running petasuite with args: {peta_args}")
            subprocess.call(peta_args)
            logger.info(f"Done with petasuite for file {filename}")
        return complete_decompressed_file_path
    else:
        logger.error(f"Unknown compression type {compression_type} for file {filename}")
        return None


def link_fastqs_to_outputdir(fastq_dict, outputdir, logger):
    """
    Link the fastq files in the dictionary to the outputdir/fastq/ directory.
    """
    if outputdir is None:
        raise ValueError("outputdir is not defined")

    fastq_dir = os.path.join(outputdir, 'fastq')
    os.makedirs(fastq_dir, exist_ok=True)

    for sample_tag, fastq_paths in fastq_dict.items():
        for fq_path in fastq_paths:
            link_name = os.path.join(fastq_dir, os.path.basename(fq_path))
            if not os.path.exists(link_name):
                logger.info(f'Linking {fq_path} to {link_name}')
                os.symlink(fq_path, link_name)
            else:
                logger.info(f'Link {link_name} already exists')

    return fastq_dir


def download_and_decompress(remote_key, logger, hcp_runtag):
    downloaded_fq = download_hcp_fq(remote_key, logger, hcp_runtag)
    decompressed_fq = decompress_downloaded_fastq(downloaded_fq, logger)
    return decompressed_fq


def find_or_download_fastqs(sample_name, logger):
    """
    If a sample name has fastqs from additional sequencing runs - fetch those fastq objects and link them to Demultiplexdir of current run. 
    """
    fq_objs = slims_connection.fetch('Content', conjunction()
                              .add(equals('cntn_id', sample_name))
                              .add(equals('cntn_fk_contentType', 22)))
                              # Removed the condition to exclude the current run_tag
    fastq_dict = {}
    if fq_objs:
        runtags = []
        for fq_obj in fq_objs:
            fqs_runtag = fq_obj.cntn_cstm_runTag.value
            runtags.append(fqs_runtag)
        with ThreadPoolExecutor() as executor:
            future_to_fq = {}
            for tag in runtags:
                fqSSample = SlimsSample(sample_name, tag)
                json_info = json.loads(fqSSample.fastq.cntn_cstm_demuxerSampleResult.value)
                fq_paths = json_info['fastq_paths']
                fq_matched = False
                for fq_path in fq_paths:
                    if os.path.exists(fq_path):
                        fq_matched = True
                        logger.info(f'Found fastq {fq_path}')
                        if f'{sample_name}_{tag}' in fastq_dict:
                            fastq_dict[f'{sample_name}_{tag}'].append(fq_path)
                        else:
                            fastq_dict[f'{sample_name}_{tag}'] = [fq_path]
                    else:
                        logger.info(f'Fastq {fq_path} does not exist. Need to download from HCP')
                        json_backup = json.loads(fqSSample.fastq.cntn_cstm_demuxerBackupSampleResult.value)
                        remote_keys = json_backup['remote_keys']
                        fq_basename_fasterq = os.path.basename(fq_path).replace('.fastq.gz', '.fasterq')
                        fq_basename_spring = os.path.basename(fq_path).replace('.fastq.gz', '.spring')
                        matching_key = [key for key in remote_keys if fq_basename_fasterq in key or fq_basename_spring in key]
                        if matching_key:
                            fq_matched = True
                            future = executor.submit(download_and_decompress, matching_key[0], logger, tag)
                            future_to_fq[future] = f'{sample_name}_{tag}'
                        else:
                            logger.info(f'No matching remote keys found for {fq_basename_fasterq} or {fq_basename_spring}')
                if not fq_matched:
                    logger.info(f"None of the remote fastqs for {sample_name}_{tag} were matched")
                    logger.info(f"Downloading all remote fastqs for {sample_name}_{tag}")
                    for remote_key in remote_keys:
                        future = executor.submit(download_and_decompress, remote_key, logger, tag)
                        future_to_fq[future] = f'{sample_name}_{tag}'
            for future in as_completed(future_to_fq):
                samplename_tag = future_to_fq[future]
                try:
                    decompressed_fq = future.result()
                    if samplename_tag in fastq_dict:
                        fastq_dict[samplename_tag].append(decompressed_fq)
                    else:
                        fastq_dict[samplename_tag] = [decompressed_fq]
                except Exception as exc:
                    logger.error(f'{tag} generated an exception: {exc}')
            logger.info(f'Found fastqs for {sample_name}_{tag}')
    return fastq_dict


def get_pair_dict(Sctx, Rctx, logger):
    """
    If tumor and normal are sequenced in different runs - find the pairs. 
    Then use the find_more_fastqs function to find paths of fastqs that are sequenced in different runs and link fastqs.
    Returns a dict of T/N info 
    """

    run_tag = Rctx.run_tag
    pair_dict = {}

    # FIXME: using equals tumorNormalID here won't work when we change it to pairID...
    Sctx.slims_info = get_sample_slims_info(Sctx, run_tag)
    pairs = slims_connection.fetch('Content', conjunction()
                              .add(equals('cntn_fk_contentType', 6))
                              .add(equals('cntn_cstm_tumorNormalID', 
                              Sctx.slims_info['tumorNormalID'])))

    # If they don't have the same tumorNormalID
    pairs2 = slims_connection.fetch('Content', conjunction()
                              .add(equals('cntn_fk_contentType', 6))
                              .add(equals('cntn_id','DNA'+Sctx.slims_info['tumorNormalID'])))
    
    # We want to make sure that we complement the tumorNormalType correctly, i.e. a tumor goes with normal and vice versa
    pair_type_d = {'tumor':'normal', 'normal':'tumor'}
    pair_type = pair_type_d.get(Sctx.slims_info['tumorNormalType'],None)
    if not pair_type:
        logger.warning(f'The sample {Sctx.slims_info["content_id"]} does not have any assigned tumorNormalType')

    for pair in pairs:
        pair_slims_sample = translate_slims_info(pair)
        # Check if the sample we have found is either our newly sequenced sample (including the same sample previously sequenced) OR a complementing tumorNormalType to our newly sequenced sample
        if pair_slims_sample['content_id'] == Sctx.slims_info['content_id'] or\
                pair_slims_sample['tumorNormalType'] == pair_type:
            pair_dict[pair_slims_sample['content_id']] = [pair_slims_sample['tumorNormalType'], pair_slims_sample['tumorNormalID'], pair_slims_sample['department'], pair_slims_sample['is_priority']]
            # Check if there are additional fastqs in other runs and symlink fastqs
    for p in pairs2:
        pair_slims_sample = translate_slims_info(p)
        if not pair_slims_sample['tumorNormalType'] == pair_type:
            continue
        if not pair_slims_sample['content_id'] in pair_dict:
            logger.info(f'Using sample {pair_slims_sample["content_id"]} as a pair to {Sctx.slims_info["content_id"]} (tumorNormalID {Sctx.slims_info["tumorNormalID"]}), but it does not have a matching tumorNormalID')
        if not pair_slims_sample['tumorNormalType']:
            logger.warning(f'The sample {pair_slims_sample["content_id"]} does not have any assigned tumorNormalType, will not be used in pairing')
        pair_dict[pair_slims_sample['content_id']] = [pair_slims_sample['tumorNormalType'], pair_slims_sample['tumorNormalID'], pair_slims_sample['department'], pair_slims_sample['is_priority']]
    return pair_dict
