import os
import re
import json
import yaml
import subprocess

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

def download_hcp_fqs(fqSSample, run_path, logger):
    '''Find and download fqs from HCP to tmp folder on medstore and then link to fastqdir on seqstore for run'''
    with open(CONFIG_PATH, 'r') as conf:
        config = yaml.safe_load(conf)

    json_info = json.loads(fqSSample.fastq.cntn_cstm_demuxerBackupSampleResult.value)
    bucket = json_info['bucket']
    remote_keys = json_info['remote_keys']

    queue = config["hcp"]["queue"]
    threads = config["hcp"]["threads"]
    qsub_script = config["hcp"]["qsub_script"]
    credentials = config["hcp"]["credentials"]

    for key in remote_keys:
        local_path = f'{run_path}/fastq/{os.path.basename(key)}'
        if not os.path.exists(local_path):
            standardout = os.path.join(ROOT_LOGGING_PATH, f"hcp_download_{os.path.basename(key)}.stdout")
            standarderr = os.path.join(ROOT_LOGGING_PATH, f"hcp_download_{os.path.basename(key)}.stderr")
            qsub_args = ["qsub", "-N", f"hcp_download_{os.path.basename(key)}", "-q", queue, "-sync", "y", "-o", standardout, "-e", standarderr, qsub_script, credentials, bucket, key, local_path]
            logger.info(f'Downloading {os.path.basename(key)} from HCP')
            subprocess.call(qsub_args)



def link_fastqs(list_of_fq_paths, run_path, fqSSample, logger):
    '''Link fastqs to fastq-folder in demultiplexdir of current run.'''
    # TODO: additional fastqs need to still be in demultiplexdir. not considering downloading from hcp right now. need to consider this later...
    for fq_path in list_of_fq_paths:
        fq_link = os.path.join(run_path, "fastq", os.path.basename(fq_path))
        if os.path.exists(fq_path): # If fq still on seqstore
        # Only links if link doesn't already exist
            if not os.path.islink(fq_link):
            # Now symlinks all additional paths to fastqs for tumor and normal in other runs. 
                os.symlink(fq_path, fq_link)
        else:
            logger.info(f'{fq_path} does not exist. Need to download from hcp')
            download_hcp_fqs(fqSSample, run_path, logger)

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
    for pair in pairs:
        pair_slims_sample = translate_slims_info(pair)
        pair_dict[pair_slims_sample['content_id']] = [pair_slims_sample['tumorNormalType'], pair_slims_sample['tumorNormalID'], pair_slims_sample['department'], pair_slims_sample['is_priority']]
        # Check if there are additional fastqs in other runs and symlink fastqs
        find_more_fastqs(pair.cntn_id.value, Rctx, logger)
    for p in pairs2:
        pair_slims_sample = translate_slims_info(p)
        pair_dict[pair_slims_sample['content_id']] = [pair_slims_sample['tumorNormalType'], pair_slims_sample['tumorNormalID'], pair_slims_sample['department'], pair_slims_sample['is_priority']]
        #FIND MORE FQS
        find_more_fastqs(p.cntn_id.value, Rctx, logger)
    return pair_dict


# Need to do:
# Paired T+N comes from different runs. Do nothing with the first sample - wait for its pair to run pipeline.
# Tumor only. Flag in samplesheet? Info about this in slims?
# Normal only. Flag in samplesheet? Info about this in slims?
# Tumor has run as tumor only in a previous run. Normal comes in a later run and needs to be paired with its tumor and run as paired.
