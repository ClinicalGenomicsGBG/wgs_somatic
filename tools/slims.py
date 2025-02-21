import os
import re
import json
import subprocess
import glob
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from slims.slims import Slims
from slims.criteria import is_one_of, equals, conjunction, not_equals
from slims.content import Status

from tools.helpers import read_config
from definitions import WRAPPER_CONFIG_PATH, ROOT_DIR

# TODO: Change to point to a credentials file in the config
# These are currently set in env/wgs_somatic_env/etc/conda/activate.d/set_slims_env.sh
class slims_credentials:
    url = os.environ.get('SLIMS_URL')
    user = os.environ.get('SLIMS_USER')
    password = os.environ.get('SLIMS_PASSWORD')


instance = 'wgs-somatic_query'
url = slims_credentials.url
user = slims_credentials.user
password = slims_credentials.password
slims_connection = Slims(instance, url, user, password)



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

def download_hcp_fq(bucket, remote_key, logger, hcp_runtag):
    '''Find and download fqs from HCP to fastqdir on seqstore for run'''
    config = read_config(WRAPPER_CONFIG_PATH)

    queue = config["hcp"]["queue"]
    qsub_script = config["hcp"]["qsub_script"]
    credentials = config["hcp"]["credentials"]
    hcp_downloads = config["hcp_download_dir"]
    wrapper_log_path = config["wrapper_log_path"]

    hcp_download_runpath = f'{hcp_downloads}/{hcp_runtag}' # This is the directory where the downloaded files will be stored

    hcp_path = f'{hcp_download_runpath}/{os.path.basename(remote_key)}' # This is the complete path of the downloaded file
    if not os.path.exists(hcp_path):

        standardout = os.path.join(wrapper_log_path, f"hcp_download_{os.path.basename(remote_key)}.stdout")
        standarderr = os.path.join(wrapper_log_path, f"hcp_download_{os.path.basename(remote_key)}.stderr")

        os.makedirs(hcp_download_runpath, exist_ok=True)

        qsub_args = ["qsub", "-N", f"hcp_download_{os.path.basename(remote_key)}", "-q", queue, "-sync", "y", "-o", standardout, "-e", standarderr, qsub_script, credentials, bucket, remote_key, hcp_path] 
        logger.info(f'Downloading {os.path.basename(remote_key)} from HCP')
        subprocess.call(qsub_args)

        # Poll for the existence of the downloaded file
        while not os.path.exists(hcp_path):
            logger.info(f'Waiting for {hcp_path} to be downloaded...')
            time.sleep(10)  # Wait for 10 seconds before checking again

    else:
        logger.info(f'{os.path.basename(remote_key)} already exists in {hcp_download_runpath}')
    return hcp_path


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
            dir_of_downloaded_file = os.path.dirname(complete_file_path)
            os.chdir(f'{dir_of_downloaded_file}') # Change to the directory of the downloaded file

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


def download_and_decompress(bucket, remote_key, logger, hcp_runtag):
    downloaded_fq = download_hcp_fq(bucket, remote_key, logger, hcp_runtag)
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
                for fq_path in fq_paths:
                    if os.path.exists(fq_path):
                        logger.info(f'Found fastq {fq_path}')
                        if f'{sample_name}_{tag}' in fastq_dict:
                            fastq_dict[f'{sample_name}_{tag}'].append(fq_path)
                        else:
                            fastq_dict[f'{sample_name}_{tag}'] = [fq_path]
                    else:
                        logger.info(f'Fastq {fq_path} does not exist. Need to download from HCP')
                        json_backup = json.loads(fqSSample.fastq.cntn_cstm_demuxerBackupSampleResult.value)
                        bucket = json_backup['bucket']
                        remote_keys = json_backup['remote_keys']
                        fq_basename_fasterq = os.path.basename(fq_path).replace('.fastq.gz', '.fasterq')
                        fq_basename_spring = os.path.basename(fq_path).replace('.fastq.gz', '.spring')
                        matching_key = [key for key in remote_keys if fq_basename_fasterq in key or fq_basename_spring in key]
                        if matching_key:
                            future = executor.submit(download_and_decompress, bucket, matching_key[0], logger, tag)
                            future_to_fq[future] = f'{sample_name}_{tag}'
                        else:
                            logger.info(f'No matching remote keys found for {fq_basename_fasterq} or {fq_basename_spring}')
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


# Need to do:
# Paired T+N comes from different runs. Do nothing with the first sample - wait for its pair to run pipeline.
# Tumor only. Flag in samplesheet? Info about this in slims?
# Normal only. Flag in samplesheet? Info about this in slims?
# Tumor has run as tumor only in a previous run. Normal comes in a later run and needs to be paired with its tumor and run as paired.
