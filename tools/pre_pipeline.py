from pathlib import Path
from typing import Optional
from tools.slims import Patient, Sample
from datetime import datetime

class Run:
    def __init__(
        self,
        logger,
        patient: Patient,
        run_root_dir: Optional[Path] = None,
        run_work_dir: Optional[Path] = None,
        main_id: Optional[str] = None,
        est_tumor_cov: Optional[float] = None,
        est_normal_cov: Optional[float] = None,
    ):
        self.patient = patient
        self.tumor_samples = self.patient.tumor_samples
        self.normal_samples = self.patient.normal_samples
        self.run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
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
        #self.tumor_name: Optional[str] = None
        #self.normal_name: Optional[str] = None

        if run_work_dir is not None:
            self.run_work_dir = run_work_dir
        elif run_root_dir is not None:
            if not self.main_id:
                self.main_id = self._determine_main_id()
            self.run_work_dir = run_root_dir / f"{self.main_id}_{self.run_timestamp}"
        else:
            logger.error("No run_root_dir or run_work_dir provided, cannot determine run_work_dir")
            raise ValueError("No run_root_dir or run_work_dir provided, cannot determine run_work_dir")

    def _determine_main_id(self) -> str:
        """Return sample_id of most recent tumor, otherwise most recent normal sample."""
        if self.patient.tumor_samples:
            return max(self.patient.tumor_samples, key=lambda s: s.date_created).id
        if seld.patient.normal_samples:
            return max(self.patient.normal_samples, key=lambda s: s.date_created).id

    #def latest_sample_id(self, sample_type: str) -> str | None:
    #    if sample_type == "tumor":
    #        samples = self.tumor_samples
    #    elif sample_type == "normal":
    #        samples = self.normal_samples
    #    else:
    #        raise ValueError(f"Unknown sample_type: {sample_type}")
#
 #       if not samples:
 #           return None
 #       return max(samples, key=lambda sample: sample.date_created).id



def pre_pipeline(runs: list[Run], config, logger) -> list[Run]:
    """
    Build and prepare run objects for pipeline submission.

    For each run:
    - ensure sample FASTQs are available locally (or downloaded)
    - merge multiple tumor/normal FASTQ sets into one R1/R2 per role
    - place final files in run_work_dir/fastq
    """
    prepared_runs: list[Run] = []
    


    for run in runs:
        
        logger.info(f"Preparing run in {run.run_work_dir}")

        run.prepared_fastq_dir = run.run_work_dir / "fastq"
        run.prepared_fastq_dir.mkdir(parents=True, exist_ok=True)
        #run.tumor_name = run.latest_sample_id("tumor")
        #run.normal_name = run.latest_sample_id("normal")
        tumor_r1_sources: list[Path] = []
        tumor_r2_sources: list[Path] = []
        normal_r1_sources: list[Path] = []
        normal_r2_sources: list[Path] = []


        for sample in run.patient.samples:
            print("Before resolving")
            print(sample)
            sample.resolve_fastq_pair(run.prepared_fastq_dir, logger)
            print("After resolving")
            print(sample)
            #if sample.type_somatic == "tumor":
            #    tumor_r1_sources.append(sample_r1)
            #    tumor_r2_sources.append(sample_r2)
            #elif sample.type_somatic == "normal":
            #    normal_r1_sources.append(sample_r1)
            #    normal_r2_sources.append(sample_r2)
            #else:
            #    logger.warning(f"Skipping sample {sample.id} with unknown somatic type: {sample.type_somatic}")

        has_tumor = run.patient.has_tumor 
        has_normal = run.patient.has_normal 
        if has_tumor and len(tumor_r1_sources) != len(tumor_r2_sources):
            raise ValueError(
                f"Run {run.run_work_dir} has unbalanced tumor FASTQs: {len(tumor_r1_sources)} R1 and {len(tumor_r2_sources)} R2"
            )

        if has_normal and len(normal_r1_sources) != len(normal_r2_sources):
            raise ValueError(
                f"Run {run.run_work_dir} has unbalanced normal FASTQs: {len(normal_r1_sources)} R1 and {len(normal_r2_sources)} R2"
            )

        if not has_tumor and not has_normal:
            logger.warning(f"Run {run.run_work_dir} has no resolved FASTQs, will not be marked ready.")
            #prepared_runs.append(run)
            continue

        if has_tumor:
            run.prepared_tumor_r1 = run.prepared_fastq_dir / f"{run.patient.tumor_name}_R1_001.fastq.gz"
            run.prepared_tumor_r2 = run.prepared_fastq_dir / f"{run.patient.tumor_name}_R2_001.fastq.gz"
            _materialize_fastq_pair(sample.r1_sources, run.prepared_tumor_r1, logger)
            _materialize_fastq_pair(tumor_r2_sources, run.prepared_tumor_r2, logger)

        if has_normal:
            run.prepared_normal_r1 = run.prepared_fastq_dir / f"{run.patient.normal_name}_R1_001.fastq.gz"
            run.prepared_normal_r2 = run.prepared_fastq_dir / f"{run.patient.normal_name}_R2_001.fastq.gz"
            _materialize_fastq_pair(normal_r1_sources, run.prepared_normal_r1, logger)
            _materialize_fastq_pair(normal_r2_sources, run.prepared_normal_r2, logger)

        #run.ready_for_pipeline = True
        outputdir = str(run.run_work_dir)

        pipeline_args = {
            "outputdir": outputdir,
        }

        tumor = run.patient.tumor_name
        normal = run.patient.normal_name

        if tumor:
            pipeline_args["tumorname"] = tumor
            pipeline_args["tumorfastqs"] = str(run.prepared_fastq_dir)

        if normal:
            pipeline_args["normalname"] = normal
            pipeline_args["normalfastqs"] = str(run.prepared_fastq_dir)

        run.pipeline_args = pipeline_args

        if tumor and normal:
            run.display_name = f"{tumor} (T) {normal} (N)"
            run.analysis_end_args = (outputdir, tumor, normal)
        elif tumor:
            run.display_name = f"{tumor} (T)"
            run.analysis_end_args = (outputdir, tumor, None)
        elif normal:
            run.display_name = f"{normal} (N)"
            run.analysis_end_args = (outputdir, None, normal)
        else:
            continue

        prepared_runs.append(run)

    return prepared_runs
