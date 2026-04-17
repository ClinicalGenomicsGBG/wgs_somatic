import shutil
from pathlib import Path

from tools.slims import Sample, Run

def _resolve_sample_fastq_pair(sample: Sample, config, logger) -> tuple[Path, Path]:
    if sample.has_local_fastqs():
        logger.info(f"Using local FASTQs from Sample record for {sample.id}")
        return sample.r1_path, sample.r2_path

    logger.info(f"FASTQs not fully local for {sample.id}, downloading from sample remote keys")
    sample.download(config=config, logger=logger)
    if not sample.has_local_fastqs():
        raise FileNotFoundError(f"Could not resolve paired FASTQs for sample {sample.id}")
    return sample.r1_path, sample.r2_path


def _materialize_fastq_pair(source_paths: list[Path], target_path: Path, logger) -> None:
    target_path.parent.mkdir(parents=True, exist_ok=True)
    if target_path.exists() or target_path.is_symlink():
        target_path.unlink()

    if len(source_paths) == 1:
        logger.info(f"Linking {source_paths[0]} -> {target_path}")
        target_path.symlink_to(source_paths[0])
        return

    logger.info(f"Merging {len(source_paths)} FASTQs into {target_path}")
    with open(target_path, "wb") as out_handle:
        for source_path in source_paths:
            with open(source_path, "rb") as in_handle:
                shutil.copyfileobj(in_handle, out_handle, length=1024 * 1024)

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
        run.prepared_fastq_dir = run.run_work_dir / "fastq"
        run.prepared_fastq_dir.mkdir(parents=True, exist_ok=True)
        run.tumor_name = run.latest_sample_id("tumor")
        run.normal_name = run.latest_sample_id("normal")

        tumor_r1_sources: list[Path] = []
        tumor_r2_sources: list[Path] = []
        normal_r1_sources: list[Path] = []
        normal_r2_sources: list[Path] = []

        logger.info(f"Preparing run in {run.run_work_dir}")

        for sample in run.samples:
            sample_r1, sample_r2 = _resolve_sample_fastq_pair(sample, config, logger)

            if sample.type_somatic == "tumor":
                tumor_r1_sources.append(sample_r1)
                tumor_r2_sources.append(sample_r2)
            elif sample.type_somatic == "normal":
                normal_r1_sources.append(sample_r1)
                normal_r2_sources.append(sample_r2)
            else:
                logger.warning(f"Skipping sample {sample.id} with unknown somatic type: {sample.type_somatic}")

        has_tumor = bool(tumor_r1_sources or tumor_r2_sources)
        has_normal = bool(normal_r1_sources or normal_r2_sources)

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
            prepared_runs.append(run)
            continue

        if has_tumor:
            run.prepared_tumor_r1 = run.prepared_fastq_dir / f"{run.tumor_name}_R1_001.fastq.gz"
            run.prepared_tumor_r2 = run.prepared_fastq_dir / f"{run.tumor_name}_R2_001.fastq.gz"
            _materialize_fastq_pair(tumor_r1_sources, run.prepared_tumor_r1, logger)
            _materialize_fastq_pair(tumor_r2_sources, run.prepared_tumor_r2, logger)

        if has_normal:
            run.prepared_normal_r1 = run.prepared_fastq_dir / f"{run.normal_name}_R1_001.fastq.gz"
            run.prepared_normal_r2 = run.prepared_fastq_dir / f"{run.normal_name}_R2_001.fastq.gz"
            _materialize_fastq_pair(normal_r1_sources, run.prepared_normal_r1, logger)
            _materialize_fastq_pair(normal_r2_sources, run.prepared_normal_r2, logger)

        run.ready_for_pipeline = True
        prepared_runs.append(run)

    return prepared_runs
