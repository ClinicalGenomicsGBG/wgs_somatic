import os
import re
import json
import yaml
import subprocess
import glob

from slims.slims import Slims
from slims.criteria import is_one_of, equals, conjunction, not_equals
from slims.content import Status

from definitions import CONFIG_PATH, ROOT_DIR, ROOT_LOGGING_PATH

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

def download_hcp_fqs(fqSSample, run_path, logger, hcp_runtag):
    '''Find and download fqs from HCP to fastqdir on seqstore for run'''
    with open(CONFIG_PATH, 'r') as conf:
        config = yaml.safe_load(conf)

    json_info = json.loads(fqSSample.fastq.cntn_cstm_demuxerBackupSampleResult.value)
    bucket = json_info['bucket']
    remote_keys = json_info['remote_keys']

    queue = config["hcp"]["queue"]
    qsub_script = config["hcp"]["qsub_script"]
    credentials = config["hcp"]["credentials"]
    hcp_downloads = config["hcp_download_dir"]
    
    hcp_download_runpath = f'{hcp_downloads}/{hcp_runtag}' # This is the directory where the downloaded files will be stored
    
    for key in remote_keys:
        local_path = f'{run_path}/fastq/{os.path.basename(key)}' # This is were the file
        hcp_path = f'{hcp_download_runpath}/{os.path.basename(key)}' # This is the complete path of the downloaded file
        if not os.path.exists(local_path) or not os.path.exists(hcp_path):
            
            standardout = os.path.join(ROOT_LOGGING_PATH, f"hcp_download_{os.path.basename(key)}.stdout")
            standarderr = os.path.join(ROOT_LOGGING_PATH, f"hcp_download_{os.path.basename(key)}.stderr")

            try:
                if not os.path.exists(hcp_download_runpath):
                    os.makedirs(f'{hcp_download_runpath}')
                
                qsub_args = ["qsub", "-N", f"hcp_download_{os.path.basename(key)}", "-q", queue, "-sync", "y", "-o", standardout, "-e", standarderr, qsub_script, credentials, bucket, key, hcp_path] 
                logger.info(f'Downloading {os.path.basename(key)} from HCP')
                subprocess.call(qsub_args)
                
                decompress_downloaded_fastq(hcp_path, logger)
                
            except FileExistsError:
                pass

def decompress_downloaded_fastq(complete_file_path, logger):
    with open(CONFIG_PATH, 'r') as conf:
        config = yaml.safe_load(conf)
    
    filename = os.path.basename(complete_file_path) # This is the filename of the downloaded file
    standardout_decompress = os.path.join(ROOT_LOGGING_PATH, f"decompress_{filename}.stdout")
    standarderr_decompress = os.path.join(ROOT_LOGGING_PATH, f"decompress_{filename}.stderr")
    
    queue = config["hcp"]["queue"]
    threads = config["hcp"]["threads"]
    compression_type = filename.split('.')[-1] # This is the compression type of the downloaded file, could be either 'spring' or 'fasterq'
    
    if compression_type == 'spring':
        # Decompress the file using spring
        
        decompress_script = os.path.join(ROOT_DIR, config["hcp"]["spring_script"])
        complete_decompressed_file_path = complete_file_path.replace('.spring', '.fastq.gz')
        
        qsub_args = ["qsub", "-N", f"decompressing_{filename}", "-q", queue, "-sync", "y", "-pe", "mpi", f"{threads}",
                     "-o", standardout_decompress, "-e", standarderr_decompress, "-v", f"THREADS={threads}",
                     decompress_script, complete_file_path, complete_decompressed_file_path, str(threads)]
        
        logger.info(f"Decompressing {filename} using spring with args: {qsub_args}")
        subprocess.call(qsub_args)
        logger.info(f"Done decompressing {filename}")
    
    elif compression_type == 'fasterq':
        # Decompress the file using petasuite
        
        decompress_script = os.path.join(ROOT_DIR, config["hcp"]["peta_script"])
        dir_of_downloaded_file = os.path.dirname(complete_file_path)
        cwd = os.getcwd() # Save current working directory
        os.chdir(f'{dir_of_downloaded_file}') # Change to the directory of the downloaded file
        
        peta_args = ["qsub", "-N", f"decompressing_file_{filename}", "-q", queue, "-sync", "y", 
                     "-pe", "mpi", f"{threads}", "-o", standardout_decompress, "-e", standarderr_decompress, "-v",f"THREADS={threads}",
                     decompress_script, str(threads)] 
        
        logger.info(f"Running petasuite with args: {peta_args}")
        subprocess.call(peta_args)
        logger.info(f"Done with petasuite for file {filename}")
        os.chdir(cwd) # Change back to the original working directory
    
    

def link_fastqs(list_of_fq_paths, run_path, fqSSample, logger):
    '''Link fastqs to fastq-folder in demultiplexdir of current run.'''
    with open(CONFIG_PATH, 'r') as conf:
        config = yaml.safe_load(conf)
    hcp_downloads = config["hcp_download_dir"]
    # TODO: additional fastqs need to still be in demultiplexdir. not considering downloading from hcp right now. need to consider this later...
    for fq_path in list_of_fq_paths:
        fq_link = os.path.join(run_path, "fastq", os.path.basename(fq_path))
#        hcp_runtag = fq_path.split("/")[-3] # This could cause problems if we change the file structure
        hcp_runtag = fqSSample.fastq.cntn_cstm_runTag.value
        hcp_path =  f'{hcp_downloads}/{hcp_runtag}'
        if os.path.exists(fq_path): # If fq still on seqstore
        # Only links if link doesn't already exist
            if not os.path.islink(fq_link):
            # Now symlinks all additional paths to fastqs for tumor and normal in other runs.
                os.symlink(fq_path, fq_link)
        #If fq not in seqstore
        elif not os.path.exists(fq_path):
            # Check /medstore/tmp/hcp_downloads 
            # If one DNA object has multiple fastq objects linked we want them to be in the same folder
            previously_downloaded_samples_same_DNA = sorted(glob.glob(os.path.join(hcp_downloads,"*",os.path.basename(fq_path).split('_')[0]+"*.gz")))
            if previously_downloaded_samples_same_DNA:
                hcp_runtag = os.path.basename(os.path.dirname(previously_downloaded_samples_same_DNA[0]))
                hcp_path = f'{hcp_downloads}/{hcp_runtag}'
                logger.info(f"Found previously downloaded file {os.path.basename(fq_path).split('_')[0]} in {hcp_path}")
                logger.info(f"Setting download directory to {hcp_path}")
            if os.path.isdir(hcp_path):
                if os.path.basename(fq_path) in os.listdir(hcp_path):
                    logger.info(f'The file {os.path.basename(fq_path)} exists in {hcp_path}')
                else:
                    logger.info(f'{fq_path} does not exist. Need to download from hcp')
                    download_hcp_fqs(fqSSample, run_path, logger, hcp_runtag)
            else:
                download_hcp_fqs(fqSSample, run_path, logger, hcp_runtag)
            
            # Link the downloaded sample to the current rundir
            downloaded_fq_path = os.path.join(hcp_path, os.path.basename(fq_path))
            logger.info(f"The downloaded file is located {downloaded_fq_path}")
            fq_link = os.path.join(run_path, "fastq", os.path.basename(fq_path))
            if os.path.exists(downloaded_fq_path):
                if not os.path.islink(fq_link):
                    logger.info(f"Linking {downloaded_fq_path} to {fq_link}")
                    os.symlink(downloaded_fq_path, fq_link)

def find_more_fastqs(sample_name, Rctx, logger):
    """
    If a sample name has fastqs from additional sequencing runs - fetch those fastq objects and link them to Demultiplexdir of current run. 
    """
    run_tag = Rctx.run_tag
    more_fastqs = slims_connection.fetch('Content', conjunction()
                              .add(equals('cntn_id', sample_name))
                              .add(equals('cntn_fk_contentType', 22))
                              .add(not_equals('cntn_cstm_runTag', run_tag)))
    if more_fastqs:
        logger.info('There are more fastqs in other sequencing runs')
        runtags = []
        for fq in more_fastqs:
            fqs_runtag = fq.cntn_cstm_runTag.value
            runtags.append(fqs_runtag)
        for tag in runtags:
            fqSSample = SlimsSample(sample_name, tag)
            json_info = json.loads(fqSSample.fastq.cntn_cstm_demuxerSampleResult.value)
            fq_paths = json_info['fastq_paths']
            logger.info(f'linking fastqs for {sample_name}_{tag}')
            link_fastqs(fq_paths, Rctx.run_path, fqSSample, logger)

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
            find_more_fastqs(pair.cntn_id.value, Rctx, logger)
    for p in pairs2:
        pair_slims_sample = translate_slims_info(p)
        if not pair_slims_sample['content_id'] in pair_dict:
            logger.info(f'Using sample {pair_slims_sample["content_id"]} as a pair to {Sctx.slims_info["content_id"]} (tumorNormalID {Sctx.slims_info["tumorNormalID"]}), but it does not have a matching tumorNormalID')
        if not pair_slims_sample['tumorNormalType']:
            logger.warning(f'The sample {pair_slims_sample["content_id"]} does not have any assigned tumorNormalType, will not be used in pairing')
        pair_dict[pair_slims_sample['content_id']] = [pair_slims_sample['tumorNormalType'], pair_slims_sample['tumorNormalID'], pair_slims_sample['department'], pair_slims_sample['is_priority']]
        #FIND MORE FQS
        find_more_fastqs(p.cntn_id.value, Rctx, logger)
    return pair_dict


# Need to do:
# Paired T+N comes from different runs. Do nothing with the first sample - wait for its pair to run pipeline.
# Tumor only. Flag in samplesheet? Info about this in slims?
# Normal only. Flag in samplesheet? Info about this in slims?
# Tumor has run as tumor only in a previous run. Normal comes in a later run and needs to be paired with its tumor and run as paired.
