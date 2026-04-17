import os
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Any
from pathlib import Path
from datetime import datetime

from slims.slims import Slims
from slims.criteria import equals, conjunction, not_equals

from tools.helpers import read_config
from definitions import WRAPPER_CONFIG_PATH, ROOT_DIR

RUN_TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")


def latest_result(results):
    return max(results, key=lambda r: r.column("rslt_createdOn").value, default=None)


class Sample:
    def __init__(self, slims_record):
        self.pk = slims_record.pk()
        self.date_created = self._value(slims_record, "rslt_createdOn")
        self.subject_barcode = self._value(slims_record, "rslt_cf_subjectBarcode")
        self.id = self._value(slims_record, "rslt_cf_sampleId")
        self.fastq_file_paths = self._value(slims_record, "rslt_cf_fastqFilePaths")
        self.long_term_storage_info = self._value(slims_record, "rslt_cf_longTermStorageInfo")
        self.family_id = self._value(slims_record, "rslt_cf_familyId")
        self.type_somatic = self._normalize_type(self._value(slims_record, "rslt_cf_sampleTypeSomatic"))
        self.sex = self._value(slims_record,"rslt_cf_sex")
        self.priority = self._value(slims_record,"rslt_cf_priority")
        self.status = self._display(slims_record,"rslt_fk_status")
        self.fastq_merge = self._value(slims_record,"rslt_cf_fqMerge")

        self.r1_path = Path(self.fastq_file_paths[0]) if self.fastq_file_paths else None
        self.r2_path = Path(self.fastq_file_paths[1]) if self.fastq_file_paths and len(self.fastq_file_paths) > 1 else None
        self.r1_linked_path = None
        self.r2_linked_path = None
        self.fastq_local = self.r1_path and self.r1_path.exists() and self.r2_path and self.r2_path.exists()

        self.r1_remote = None  # FIXME: Some info from the remote key field in slims should go here to be used for downloading if the files are not local
        self.r2_remote = None

    def missing_required_fields(self) -> list[str]:
        missing: list[str] = []
        if not self.id:
            missing.append("id")
        if not self.subject_barcode:
            missing.append("subject_barcode")
        if not self.date_created:
            missing.append("date_created")
        if not self.type_somatic:
            missing.append("type_somatic")
        return missing

    def validate_required_fields(self) -> None:
        missing = self.missing_required_fields()
        if missing:
            raise ValueError(f"Sample is missing required fields: {', '.join(missing)}")

    def _existing_fastq_paths(self) -> list[Path]:
        if not self.fastq_file_paths:
            return []
        return [Path(path) for path in self.fastq_file_paths if path and Path(path).exists()]

    def _refresh_fastq_state(self) -> None:
        if self.fastq_file_paths and len(self.fastq_file_paths) > 1:
            r1_candidate = Path(self.fastq_file_paths[0])
            r2_candidate = Path(self.fastq_file_paths[1])
            self.r1_path = r1_candidate if r1_candidate.exists() else None
            self.r2_path = r2_candidate if r2_candidate.exists() else None
        else:
            self.r1_path = None
            self.r2_path = None
        self.fastq_local = bool(self.r1_path and self.r2_path and self.r1_path.exists() and self.r2_path.exists())

    def has_local_fastqs(self) -> bool:
        self._refresh_fastq_state()
        return self.fastq_local

    def _parse_long_term_storage_info(self) -> tuple[list[str], str | None]:
        info: Any = self.long_term_storage_info
        if not info:
            return [], None

        if isinstance(info, str):
            try:
                info = json.loads(info)
            except json.JSONDecodeError:
                return [], None

        remote_keys: list[str] = []
        bucket: str | None = None

        def _collect(obj: Any) -> None:
            nonlocal bucket
            if isinstance(obj, dict):
                keys = obj.get("remote_keys")
                if isinstance(keys, list):
                    remote_keys.extend([key for key in keys if isinstance(key, str)])
                if bucket is None and isinstance(obj.get("bucket"), str):
                    bucket = obj["bucket"]
                for value in obj.values():
                    _collect(value)
            elif isinstance(obj, list):
                for item in obj:
                    _collect(item)

        _collect(info)
        deduped_remote_keys = list(dict.fromkeys(remote_keys))
        return deduped_remote_keys, bucket

    def has_remote_keys(self) -> bool:
        remote_keys, _ = self._parse_long_term_storage_info()
        return bool(remote_keys)

    def download(self, config, logger, hcp_runtag: str | None = None) -> list[Path]:
        if self.has_local_fastqs():
            return self._existing_fastq_paths()

        remote_keys, bucket = self._parse_long_term_storage_info()
        downloaded_paths: list[Path] = []
        runtag = hcp_runtag or self.subject_barcode or self.id

        if not remote_keys:
            raise ValueError(
                f"Sample {self.id} has no local FASTQs and no remote_keys in long_term_storage_info"
            )

        logger.info(f"Downloading missing FASTQs for sample {self.id} from {len(remote_keys)} remote keys")
        with ThreadPoolExecutor() as executor:
            for downloaded in executor.map(
                lambda remote_key: download_and_decompress(bucket, remote_key, logger, runtag), remote_keys
            ):
                if downloaded:
                    downloaded_paths.append(Path(downloaded))

        existing = self._existing_fastq_paths()
        merged = existing + [path for path in downloaded_paths if path.exists()]
        deduped = [str(path) for path in dict.fromkeys(str(path) for path in merged)]
        self.fastq_file_paths = deduped
        self._refresh_fastq_state()
        return [Path(path) for path in self.fastq_file_paths]


    def _value(self, record, name, default=None):
        try:
            return record.column(name).value
        except Exception:
            return default

    def _display(self, record, name, default=None):
        try:
            return getattr(record.column(name), "displayValue", default)
        except Exception:
            return default

    def _normalize_type(self, sample_type: str | None) -> str | None:
        if not sample_type:
            return None
        sample_type = sample_type.lower()
        if sample_type in {"tumor", "tumour"}:
            return "tumor"
        if sample_type == "normal":
            return "normal"
        return None


class Run:
    def __init__(
        self,
        logger,
        samples: list[Sample],
        run_root_dir: Optional[Path] = None,
        run_work_dir: Optional[Path] = None,
        main_id: Optional[str] = None,
        est_tumor_cov: Optional[float] = None,
        est_normal_cov: Optional[float] = None,
    ):
        self.samples = samples
        self.tumor_samples = [s for s in samples if s.type_somatic == "tumor"]
        self.normal_samples = [s for s in samples if s.type_somatic == "normal"]
        self.run_root_dir = run_root_dir
        self.main_id = main_id
        self.est_tumor_cov = est_tumor_cov
        self.est_normal_cov = est_normal_cov
        self.ready_for_pipeline = False
        self.prepared_fastq_dir: Optional[Path] = None
        self.prepared_tumor_r1: Optional[Path] = None
        self.prepared_tumor_r2: Optional[Path] = None
        self.prepared_normal_r1: Optional[Path] = None
        self.prepared_normal_r2: Optional[Path] = None
        self.tumor_name: Optional[str] = None
        self.normal_name: Optional[str] = None

        if run_work_dir is not None:
            self.run_work_dir = run_work_dir
        elif run_root_dir is not None:
            if not self.main_id:
                self.main_id = self._determine_main_id(self.tumor_samples, self.normal_samples)
            self.run_work_dir = run_root_dir / f"{self.main_id}_{RUN_TIMESTAMP}"
        else:
            logger.error("No run_root_dir or run_work_dir provided, cannot determine run_work_dir")
            raise ValueError("No run_root_dir or run_work_dir provided, cannot determine run_work_dir")

    def _determine_main_id(
        self,
        tumor_samples: list[Sample],
        normal_samples: list[Sample],
    ) -> str:
        """Return sample_id of most recent tumor, otherwise most recent normal sample."""
        if tumor_samples:
            return max(tumor_samples, key=lambda s: s.date_created).id
        if normal_samples:
            return max(normal_samples, key=lambda s: s.date_created).id

        raise ValueError("No tumor or normal samples found")

    def latest_sample_id(self, sample_type: str) -> str | None:
        if sample_type == "tumor":
            samples = self.tumor_samples
        elif sample_type == "normal":
            samples = self.normal_samples
        else:
            raise ValueError(f"Unknown sample_type: {sample_type}")

        if not samples:
            return None
        return max(samples, key=lambda sample: sample.date_created).id
        

def return_pending_samples(config, logger) -> list[Sample]:
    # Query slims for pending samples based on filters in config
    query = conjunction()
    query.add(equals("test_name", "test_pipeline_somatic"))
    query.add(equals("rslt_value", "Pending"))
    slims_records = slims_connection.fetch('Result', query)

    pending_samples = []
    for slims_record in slims_records:
        sample = Sample(slims_record)
        try:
            sample.validate_required_fields()
        except ValueError as exc:
            logger.warning(f"Skipping sample due to missing required fields: {exc}")
            continue
        # if not sample.fastq_merge:  # FIXME: Change to "StartPipeline" or similar
        #     logger.info(f"Skipping sample {sample.id} because it is not marked with StartPipeline")
        #     logger.info(f"Setting {sample.id} status to Successfull")
        #     r.update({"rslt_value": "Successfull"})
        #     continue
        logger.info(f"Setting {sample.id} status to In progress")
        #r.update({"rslt_value": "In progress"})
        pending_samples.append(sample)

    return pending_samples


def add_matched_samples(samples: list[Sample], config, logger) -> list[Sample]:
    # For all samples, check if there is another sample with oppposite tumorNormalType but same subject barcode
    # If not, we query for samples with the same subject barcode but opposite tumorNormalType and add those to the list of samples to process.
    barcodes = {"tumor": set(), "normal": set()}
    
    for sample in samples:
        if sample.type_somatic is None:
            logger.warning(f"Sample {sample.id} has unrecognized sample_type_somatic {sample.type_somatic}, skipping")
            continue
        barcodes[sample.type_somatic].add(sample.subject_barcode)

    new_samples = []

    for sample in samples:
        if sample.type_somatic is None:
            continue

        opposite_type = "normal" if sample.type_somatic == "tumor" else "tumor"

        if sample.subject_barcode in barcodes[opposite_type]:
            logger.info(f"Sample {sample.id} already has matching {opposite_type} for barcode {sample.subject_barcode}")
            continue

        logger.info(f"Missing {opposite_type} for barcode {sample.subject_barcode}, querying SLIMS")

        query = conjunction()
        query.add(equals("rslt_cf_subjectBarcode", sample.subject_barcode))
        query.add(equals("rslt_cf_sampleTypeSomatic", opposite_type))
        query.add(not_equals("pk", sample.pk))

        slims_record = latest_result(slims_connection.fetch("Result", query))

        if not slims_record:
            logger.info(f"No matching {opposite_type} sample found for barcode {sample.subject_barcode}")
            continue

        matched_sample = Sample(slims_record)
        new_samples.append(matched_sample)

        logger.info(f"Added matched {opposite_type} sample {matched_sample.id} for barcode {matched_sample.subject_barcode}")

    samples.extend(new_samples)
    return samples


def add_merge_samples(samples: list[Sample], config, logger) -> list[Sample]:
    merge_samples: list[Sample] = []
    existing_pks = {sample.pk for sample in samples}

    for sample in samples:
        if not sample.fastq_merge:
            logger.debug(f"Sample {sample.id} is not marked for fastq merging, skipping")
            continue

        query = conjunction()
        query.add(equals("rslt_cf_subjectBarcode", sample.subject_barcode))
        query.add(equals("rslt_cf_sampleTypeSomatic", sample.type_somatic))

        merge_results = slims_connection.fetch("Result", query)

        for slims_record in merge_results:
            merge_sample = Sample(slims_record)
            if merge_sample.pk in existing_pks:
                continue
            logger.info(f"Adding merge sample {merge_sample.pk} for sample {sample.id}")
            merge_samples.append(merge_sample)

    samples.extend(merge_samples)
    return samples


class SlimsCredentials:
    def __init__(self, slims_credentials_path):
        config = read_config(slims_credentials_path)
        self.url = config['slims']['url']
        self.user = config['slims']['user']
        self.password = config['slims']['password']


config = read_config(WRAPPER_CONFIG_PATH)
slims_credentials_path = config['slims']['credentials_path']
credentials = SlimsCredentials(slims_credentials_path)

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


def download_hcp_fq(bucket, remote_key, logger, hcp_runtag):
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
            # Remove getting the bucket from the config if we want to get it from slims
            bucket = location_details["bucket"]  # Bucket is specified in the config

            logger.info(f"Trying to download from {location_name}")

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
                # In rare cases, there may be a delay post-download before the file appears in the directory
                start_time = time.time()
                while not os.path.exists(hcp_path):
                    logger.info(f'Waiting for {hcp_path} to be downloaded...')
                    time.sleep(10)  # Wait for 10 seconds before checking again

                    # The delay should not be more than a minute
                    elapsed_time = time.time() - start_time
                    if elapsed_time > 60:
                        logger.error(f"The hcp_download finished successfully, but no file was found at {hcp_path}")
                        raise RuntimeError(f"The hcp_download finished successfully, but no file was found at {hcp_path}")

                logger.info(f"Successfully downloaded {os.path.basename(remote_key)} from {location_name}")
                return hcp_path

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
                        # Get bucket from slims. May or may not be correct depending on how the data is migrated
                        # bucket = json_backup['bucket']
                        # For now we will use the bucket from the config
                        bucket = None
                        remote_keys = json_backup['remote_keys']
                        fq_basename_fasterq = os.path.basename(fq_path).replace('.fastq.gz', '.fasterq')
                        fq_basename_spring = os.path.basename(fq_path).replace('.fastq.gz', '.spring')
                        matching_key = [key for key in remote_keys if fq_basename_fasterq in key or fq_basename_spring in key]
                        if matching_key:
                            fq_matched = True
                            future = executor.submit(download_and_decompress, bucket, matching_key[0], logger, tag)
                            future_to_fq[future] = f'{sample_name}_{tag}'
                        else:
                            logger.info(f'No matching remote keys found for {fq_basename_fasterq} or {fq_basename_spring}')
                if not fq_matched:
                    logger.info(f"None of the remote fastqs for {sample_name}_{tag} were matched")
                    logger.info(f"Downloading all remote fastqs for {sample_name}_{tag}")
                    for remote_key in remote_keys:
                        future = executor.submit(download_and_decompress, bucket, remote_key, logger, tag)
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
    if Sctx.slims_info['tumorNormalID'] is None:
        logger.error(f"Sample {Sctx.slims_info['content_id']} does not have a tumorNormalID assigned in SLIMS, stopping execution.")
        logger.info(f"Slims info: {Sctx.slims_info}")
        raise ValueError(f"Sample {Sctx.slims_info['content_id']} does not have a tumorNormalID assigned in SLIMS.")
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
