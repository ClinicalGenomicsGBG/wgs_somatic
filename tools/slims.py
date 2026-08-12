import os
import json
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Any
from pathlib import Path
import shutil

from slims.slims import Slims
from slims.criteria import equals, conjunction, not_equals, is_one_of

from tools.helpers import read_config
from definitions import WRAPPER_CONFIG_PATH, ROOT_DIR

def latest_result(results):
    return max(results, key=lambda r: r.column("rslt_createdOn").value, default=None)

class Patient:
    def __init__(self, subject_barcode: str):
        self.subject_barcode = subject_barcode
        self.tumor_samples: list["Sample"] = []
        self.normal_samples: list["Sample"] = []

    def __repr__(self):
        return (f"{self.subject_barcode!r}")

    def add_sample(self, sample: "Sample") -> None:
        if sample.type_somatic == "tumor":
            self.tumor_samples.append(sample)
        elif sample.type_somatic == "normal":
            self.normal_samples.append(sample)
        else:
            raise ValueError(
                f"Unknown sample_type_somatic '{sample.type_somatic}' "
                f"for sample {sample.id}"
            )
        if self.sequencing_id is None:
            self.sequencing_id = sample.sequencing_id
        elif self.sequencing_id != sample.sequencing_id:
            raise ValueError(
                f"Patient {self.subject_barcode} has samples from different sequencing runs: "
                f"{self.sequencing_id} and {sample.sequencing_id}"
            )


    @property
    def tumor_name(self) -> str | None:
        return self.tumor_samples[0].id if self.tumor_samples else None

    @property
    def normal_name(self) -> str | None:
        return self.normal_samples[0].id if self.normal_samples else None

    @property
    def has_sample_pk(self, pk: int) -> bool:
        return any(sample.pk == pk for sample in self.samples)

    @property
    def samples(self) -> list["Sample"]:
        return self.tumor_samples + self.normal_samples

    @property
    def has_tumor(self) -> bool:
        return bool(self.tumor_samples)

    @property
    def has_normal(self) -> bool:
        return bool(self.normal_samples)

    @property
    def missing_sample_types(self) -> list[str]:
        missing = []
        if not self.has_tumor:
            missing.append("tumor")
        if not self.has_normal:
            missing.append("normal")
        return missing

    def validate_sample_setup(self) -> None:
        if len(self.tumor_samples) != 1:
            raise ValueError(
                f"Patient {self.subject_barcode} has {len(self.tumor_samples)} tumor samples; "
                "exactly one is required."
            )
        if len(self.normal_samples) > 1:
            raise ValueError(
                f"Patient {self.subject_barcode} has {len(self.normal_samples)} normal samples; "
                "at most one is supported."
            )

class Sample:
    def __init__(self, slims_record):
        self.pk = slims_record.pk()
        self.date_created = self._value(slims_record, "rslt_createdOn")
        self.subject_barcode = self._value(slims_record, "rslt_cf_subjectBarcode")
        self.id = self._value(slims_record, "rslt_cf_sampleId")
        self.sequencing_id = self._value(slims_records, "rslt_cf_sequencingRunId")
        self.local_fastq_file_paths = self._parse_fastq_file_paths(
            self._value(slims_record, "rslt_cf_pipelineFilePaths")
        )

        #self.long_term_storage_info = self._value(slims_record, "rslt_cf_longTermStorageInfo")
        self.long_term_storage_info = self._value(slims_record, "rslt_cf_fastqRemotePaths")
        self.remote_keys, self.bucket = self._parse_long_term_storage_info()
        self.family_id = self._value(slims_record, "rslt_cf_familyId")
        self.type_somatic = self._normalize_type(self._value(slims_record, "rslt_cf_sampleTypeSomatic"))
        self.sex = self._value(slims_record,"rslt_cf_sex")
        self.priority = self._value(slims_record,"rslt_cf_priority")
        self.status = self._display(slims_record,"rslt_fk_status")
        self.fastq_merge = self._value(slims_record,"rslt_cf_fqMerge")

        self.r1_local = Path(self.local_fastq_file_paths[0]) if self.local_fastq_file_paths else None
        self.r2_local = Path(self.local_fastq_file_paths[1]) if self.local_fastq_file_paths and len(self.local_fastq_file_paths) > 1 else None
        self.r1_paths: list[Path] = [self.r1_local] if self.has_local_fastq else []
        self.r2_paths: list[Path] = [self.r2_local] if self.has_local_fastq else []
        self.r1_remote: Path | None = None
        self.r2_remote: Path | None = None


    def __repr__(self):
        return (f'''
        Sample: {self.id}
        pk: {self.pk}
        Local fastqs: {self.local_fastq_file_paths}
        Has local fastq: {self.has_local_fastq}
        Has remote fastq: {self.has_remote_fastq}
        Has final fastq: {self.has_final_fastq}
        Paths:
        R1:{self.r1_path}
        R2:{self.r2_path}
        ''')

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

    @staticmethod
    def _parse_fastq_file_paths(value):
        if value is None:
            return []

        if isinstance(value, list):
            return value

        if isinstance(value, str):
            value = value.strip()

            if not value:
                return []

            try:
                parsed = json.loads(value)
            except json.JSONDecodeError:
                # Fallback: treat a plain string as one path
                return [value]

            if isinstance(parsed, list):
                return parsed

            if isinstance(parsed, str):
                return [parsed]

        return []

    @property
    def has_local_fastq(self) -> bool:
        return (
            self.r1_local is not None
            and self.r2_local is not None
            and self.r1_local.exists()
            and self.r2_local.exists()
            )

    @property
    def has_remote_fastq(self) -> bool:
        return (
            self.r1_remote is not None
            and self.r2_remote is not None
            and self.r1_remote.exists()
            and self.r2_remote.exists()
            )

    @property
    def has_final_fastq(self) -> bool:
        return (
            bool(self.r1_paths)
            and bool(self.r2_paths)
            and all(path.exists() for path in self.r1_paths)
            and all(path.exists() for path in self.r2_paths)
            )
    
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

    def check_previous_merge(self, config):
        runs = sorted(
            "_".join(Path(key).stem.split("_")[1:3])
            for key in self.remote_keys
        )
        run_string = "+".join(runs)

        hcp_download_runpath = f'{config["hcp_download_dir"]}/{self.id}'

        r1_merged = (Path(hcp_download_runpath)/ f"{self.id}_{run_string}_R1_001.fastq.gz")
        r2_merged = (Path(hcp_download_runpath) / f"{self.id}_{run_string}_R2_001.fastq.gz")

        if r1_merged.exists():
            self.r1_remote = r1_merged
        if r2_merged.exists():
            self.r2_remote = r2_merged

    def download(self, config, logger) -> list[Path]:
        
        downloaded_paths: list[Path] = []

        if not self.remote_keys:
            raise ValueError(
                f"Sample {self.id} has no remote_keys in long_term_storage_info"
            )
        
        logger.info(f"Preapring missing FASTQs for sample {self.id} from {len(self.remote_keys)} remote keys")
        
        with ThreadPoolExecutor() as executor:
            for downloaded in executor.map(
                lambda remote_key: download_and_decompress(self.bucket, remote_key, logger, self.id), self.remote_keys
            ):
                if downloaded:
                    downloaded_paths.append(Path(downloaded))

        downloaded = [path for path in downloaded_paths if path.exists()]

        for read in ("R1", "R2"):
            paths = [
                path for path in downloaded
                if f"_{read}_" in path.name
            ]

            setattr(self, f"{read.lower()}_remote", paths)
            setattr(self, f"{read.lower()}_paths", paths)

            # This is kept outcommented for now. Use if a merge of samples are wanted 
            #setattr(
            #    self,
            #    f"{read.lower()}_remote",
            #    self.merge_fastqs(paths, read, logger),
            #)


    def merge_fastqs(self, source_paths: list[Path], read, logger,) -> list[Path]:

        runs = sorted("_".join(path.name.split("_")[1:3])for path in source_paths)
        run_string = "+".join(runs)
        merged_name = (f"{self.id}_{run_string}_{read}_001.fastq.gz")
        merged_path = source_paths[0].with_name(f"{self.id}_{run_string}_{read}_001.fastq.gz")
        if not merged_path.exists():
            logger.info(f"Merging {len(source_paths)} FASTQs into {merged_path}")

            with open(merged_path, "wb") as out_handle:
                for source_path in source_paths:
                    with open(source_path, "rb") as in_handle:
                        shutil.copyfileobj(
                            in_handle,
                            out_handle,
                            length=1024 * 1024,
                        )

        return merged_path


    def resolve_fastq_pair(self, run, config, logger) -> None:
        """Determine where to find fastq files to use"""

        if self.has_local_fastq:
            pass
        else:
            logger.info(f"sample {self.id} has no local FASTQ-files")

            if not self.has_remote_fastq:
                self.download(config=config, logger=logger)


        if not self.has_final_fastq:
            raise FileNotFoundError(f"could not resolve paired FASTQ files for sample {self.id}")


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

        

def return_pending_patients(config, logger) -> list[Patient]:
    # Query slims for pending samples based on filters in config
    # This is not a current functionality, remove config?
    query = conjunction()
    query.add(equals("test_name", "test_pipeline_somatic"))
    query.add(equals("rslt_value", "Pending"))

    slims_records = slims_connection.fetch('Result', query)

    patients: dict[str, Patient] = {}
    
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
        #ogger.info(f"Setting {sample.id} status to In progress")
        #slims_record.update({"rslt_value": "In progress"})
        #inding_samples.append(sample)
        
        patient_barcode = sample.subject_barcode

        patient = patients.setdefault(
            patient_barcode,
            Patient(patient_barcode)
        )

        patient.add_sample(sample)
    return list(patients.values())


def add_matched_samples(patients: list[Patient], logger) -> None:
    # For all patients, check if there are both tumor and normal samples
    # If not, we query for samples for that patient 
    for patient in patients:
        for sample_type in patient.missing_sample_types:

            logger.info(f"Missing {sample_type} sample for patient {patient.barcode}, querying SLIMS")

            query = conjunction()
            query.add(equals("rslt_cf_subjectBarcode", patient.barcode))

            if sample_type == "tumor":
                query.add(is_one_of("rslt_cf_sampleTypeSomatic", ["tumor", "tumour", "Tumor", "Tumour"]))
            else:
                query.add(is_one_of("rslt_cf_sampleTypeSomatic", ["normal", "Normal"]))

            slims_record = latest_result(slims_connection.fetch("Result", query))

            if slims_record is None:
                logger.info(f"No matching {sample_type} sample found for patient {patient.barcode}")
                continue

            patient.add_sample(Sample(slims_record))

def add_merge_samples(patients: list[Patient], logger) -> list[Patient]:
    for patient in patients:
        for sample in patient.samples:
            if not sample.fastq_merge:
                logger.debug(f"Sample {sample.id} is not marked for fastq merging, skipping")
                continue

            query = conjunction()
            query.add(equals("rslt_cf_subjectBarcode", patiente.barcode))
            query.add(equals("rslt_cf_sampleTypeSomatic", sample.type_somatic))

            merge_results = slims_connection.fetch("Result", query)

            for slims_record in merge_results:
                merge_sample = Sample(slims_record)
                
                if patient.has_sample_pk(merge_sample.pk):
                    continue

                #if merge_sample.pk in existing_pks:
                #    continue
                logger.info(f"Adding merge sample {merge_sample.pk} for patient {patient.barcode}")
                patient.add_sample(merge_sample)
                #merge_samples.append(merge_sample)


class SlimsCredentials:
    def __init__(self, slims_credentials_path):
        config = read_config(slims_credentials_path)
        self.url = config['slims']['url']
        self.user = config['slims']['user']
        self.password = config['slims']['password']


config = read_config(WRAPPER_CONFIG_PATH)
slims_credentials_path = config['slims_credentials_path']
credentials = SlimsCredentials(slims_credentials_path)

slims_connection = Slims('wgs-somatic_query',
                         credentials.url,
                         credentials.user,
                         credentials.password)


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

