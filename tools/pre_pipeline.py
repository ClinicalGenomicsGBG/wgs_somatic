from pathlib import Path
from typing import Optional
from tools.slims import Patient, Sample
from datetime import datetime
from launch_snakemake import get_timestamp

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
        self.tumor_sample = self.patient.tumor_samples[0] if self.patient.tumor_samples else None
        self.normal_sample = self.patient.normal_samples[0] if self.patient.normal_samples else None
        self.run_timestamp = datetime.now().strftime("%y%m%d-%H%M%S")
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

    def materialize_fastq(self, sample, logger) -> None:

        for source_path in sample.r1_paths + sample.r2_paths:
            target_path = self.prepared_fastq_dir / source_path.name

            if target_path.exists() or target_path.is_symlink():
                logger.info(f"Link already exists: {target_path}")
                continue

            logger.info(f"Linking {source_path} -> {target_path}")
            target_path.symlink_to(source_path)


def pre_pipeline(runs: list[Run], config, logger) -> None:
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

        for sample in run.patient.samples:
            sample.resolve_fastq_pair(run.prepared_fastq_dir, config, logger)
            run.materialize_fastq(sample, logger)

        #run.patient.validate_sample_setup()

        if not run.tumor_sample.has_final_fastq:
            raise FileNotFoundError(
                f"Tumor sample {tumor.id} is missing paired FASTQs."
            )

        if run.normal_sample and not run.normal_sample.has_final_fastq:
            raise FileNotFoundError(
                f"Normal sample {normal.id} is missing paired FASTQs."
            )

        run.pipeline_args = {
                "outputdir": str(run.run_work_dir),
                "tumorname": run.tumor_sample.id,
                "tumorfastqs": str(run.prepared_fastq_dir),
                #"tumorfastq1": str(run.tumor_sample.r1_path), This is probably a better solution
                #"tumorfastq2": str(run.tumor_sample.r2_path), This is probably a better solution
                }
        if run.normal_sample:
            run.pipeline_args["normalname"] = run.normal_sample.id
            run.pipeline_args["normalfastqs"] = str(run.prepared_fastq_dir)
            #run.pipeline_args["normalfastq1"] = str(run.normal_sample.r1_path) This is probably a better solution
            #run.pipeline_args["normalfastq2"] = str(run.normal_sample.r2_path) This is probably a better solution

    if run.patient.has_normal:
        run.display_name = f"{run.tumor_sample.id} (T) {run.normal_sample.id} (N)"
        run.analysis_end_args = (
            run.tumor_sample.id,
            run.normal_sample.id,
        )
    else:
        run.display_name = f"{run.tumor_sample.id} (T)"
        run.analysis_end_args = (
            run.tumor_sample.id,
            None,
        )

    run.ready_for_pipeline = True
