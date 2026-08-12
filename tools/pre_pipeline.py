from pathlib import Path
from typing import Optional
from tools.slims import Patient, Sample
from collections import Counter
from datetime import datetime
from launch_snakemake import get_timestamp

class Batch:
    def __init__(self):
        self.runs = []
        self.batch_name = None
        self.threads = []
        self.end_threads = []
        self.run_summary = []
        self.ok_samples = []
        self.bad_samples = []
    def determine_batch_name(self) -> None:
        sequencing_ids = [
            sample.sequencing_id
            for run in self.runs
            for sample in run.patient.samples
            if sample.sequencing_id
        ]

        most_common_id, _ = Counter(sequencing_ids).most_common(1)[0]
         
        today = datetime.now().strftime("%y%m%d")

        if most_common_id.startswith(today):
            self.batch_name = most_common_id
        else:
            sample_count = sum(len(run.patient.samples) for run in self.runs)
            self.batch_name = f'{today}_{sample_count}_samples'
    

class Run:
    def __init__(
        self,
        logger,
        patient: Patient,
        run_root_dir: OptionaloPath] = None,
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
        self.run_work_dir = Path(run_work_dir) if run_work_dir else None
        self.est_tumor_cov = est_tumor_cov
        self.est_normal_cov = est_normal_cov
        self.ready_for_pipeline = False
        self.prepared_fastq_dir: Optional[Path] = None
        self.prepared_tumor_r1: Optional[Path] = None
        self.prepared_tumor_r2: Optional[Path] = None
        self.prepared_normal_r1: Optional[Path] = None
        self.prepared_normal_r2: Optional[Path] = None
        self.setup_folders()

    def setup_folders(self):
        if self.run_work_dir:
            return
        elif self.run_root_dir is not None:
            self.main_id = self.main_id or self._determine_main_id()
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

    @property
    def check_ok(outputdir) -> bool:
    '''Function to check if analysis has finished correctly'''
        if os.path.isfile(f"{self.run_work_dir}/workflow_finished.txt"):
            return True
        else:
            return False
    
    @property
    def has_tumor(self) -> bool:
        return bool(self.tumor_sample)

    @property
    def has_normal(self) -> bool:
        return bool(self.normal_sample)

def pre_pipeline(runs: list[Run], config, logger) -> None:
    """
    Build and prepare run objects for pipeline submission.

    For each run:
    - ensure sample FASTQs are available locally (or downloaded)
    - merge multiple tumor/normal FASTQ sets into one R1/R2 per role
    - place final files in run_work_dir/fastq
    """
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
                }

        if run.normal_sample:
            run.pipeline_args["normalname"] = run.normal_sample.id
            run.pipeline_args["normalfastqs"] = str(run.prepared_fastq_dir)

            run.display_name = (f"{run.tumor_sample.id} (T) {run.normal_sample.id} (N)")

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
