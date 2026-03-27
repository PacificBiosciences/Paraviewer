#!/usr/bin/env python

import logging
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from os import path
from typing import Dict, List, Optional, Tuple

from paraviewer.orographer_builder import generate_orographer_plots
from paraviewer.page_builder import build_review_page
from paraviewer.process_paraphase import (
    get_paraphase_results,
    make_table_entries,
    make_trio_table_entries,
)
from paraviewer.process_puretarget import get_puretarget_results
from paraviewer.utils import (
    OROGRAPHER_OUTPUT_PATH,
    ParaphaseResults,
    PedigreeEntry,
    RegionEntry,
    find_vcf_file,
    make_output_dirs,
    split_bam,
)
from paraviewer.workers import (
    SinglePlotTask,
    SplitBamArgs,
    TrioPlotTask,
    generate_single_plot_worker,
    generate_trio_plot_worker,
    split_bam_worker,
)

logger = logging.getLogger(__name__)
logging.getLogger("urllib3").setLevel(logging.CRITICAL)

TASK_TIMEOUT_SECONDS = 600


def get_trio_samples(
    pedigree_dict: Dict[str, PedigreeEntry],
    all_paraphase_results: Dict[str, ParaphaseResults],
) -> Tuple[Dict[str, PedigreeEntry], set]:
    """
    Get the trio samples from the pedigree dictionary and the paraphase results.

    Returns:
        Dict mapping trio name ("{proband_id}-trio") to PedigreeEntry;
        set of all sample IDs in trios.
    """
    trio_samples = {}
    trio_ids = set()
    for sample in pedigree_dict:
        if sample in all_paraphase_results:
            pedigree_entry = pedigree_dict[sample]
            if (
                pedigree_entry.PaternalID in all_paraphase_results
                and pedigree_entry.MaternalID in all_paraphase_results
            ):
                trio_samples[f"{sample}-trio"] = pedigree_entry
                trio_ids.add(sample)
                trio_ids.add(pedigree_entry.PaternalID)
                trio_ids.add(pedigree_entry.MaternalID)
    return trio_samples, trio_ids


def read_pedigree_file(
    ped_file: Optional[str],
    include_only_samples: list,
    exclude_samples: list,
) -> Dict[str, PedigreeEntry]:
    """Read a GATK-format PED file and return a dictionary of PedigreeEntry objects."""
    if not ped_file or not path.exists(ped_file):
        return {}

    pedigree_dict = {}
    try:
        with open(ped_file, "r") as ped_handle:
            for line in ped_handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                fields = line.split()
                if len(fields) != 6:
                    logger.warning("Skipping malformed PED line: %s", line)
                    continue

                (
                    family_id,
                    individual_id,
                    paternal_id,
                    maternal_id,
                    sex,
                    phenotype,
                ) = fields
                if (
                    len(include_only_samples) > 0
                    and individual_id.lower() not in include_only_samples
                ):
                    continue
                if individual_id.lower() in exclude_samples:
                    continue

                paternal_id = paternal_id if paternal_id != "0" else ""
                maternal_id = maternal_id if maternal_id != "0" else ""

                if sex == "1":
                    sex = "Male"
                elif sex == "2":
                    sex = "Female"
                else:
                    sex = "Unknown"

                entry = PedigreeEntry(
                    FamilyID=family_id,
                    IndividualID=individual_id,
                    PaternalID=paternal_id,
                    MaternalID=maternal_id,
                    Sex=sex,
                    Phenotype=phenotype,
                )
                pedigree_dict[individual_id] = entry

    except (OSError, ValueError) as read_error:
        logger.error("Error reading PED file %s: %s", ped_file, read_error)
        sys.exit(1)

    return pedigree_dict


def process_individual_sample(
    sample_paraphase_results: ParaphaseResults,
    pedigree_dict: Dict[str, PedigreeEntry],
    outdir: str,
    clobber: bool,
    include_only_regions: list,
    exclude_regions: list,
    max_reads_per_haplotype: int,
    reference_path: str,
    gtf_file: Optional[str],
    input_dir: str,
    is_puretarget: bool,
    is_trio_sample: bool,
) -> Tuple[List[RegionEntry], dict]:
    """
    Process one sample: split BAM, build table entries, generate orographer plots.
    Used by tests and kept for compatibility; main pipeline uses workers.
    """
    make_output_dirs(outdir, sample_paraphase_results.Sample, clobber)
    split_bams = split_bam(
        sample_paraphase_results.BAM,
        sample_paraphase_results.BAI,
        outdir,
        sample_paraphase_results.Sample,
        include_only_regions,
        exclude_regions,
        max_reads_per_haplotype,
    )
    if not split_bams:
        raise ValueError(
            "No specified regions found in %s." % sample_paraphase_results.BAM
        )

    sample_region_entries = make_table_entries(
        sample_paraphase_results,
        pedigree_dict.get(sample_paraphase_results.Sample),
        split_bams,
        is_trio_sample,
    )
    successful_entries = generate_orographer_plots(
        sample_region_entries,
        outdir,
        reference_path,
        gtf_file,
        input_dir,
        is_puretarget,
    )
    return successful_entries, split_bams


def process_trio(
    trio: PedigreeEntry,
    all_paraphase_results: Dict[str, ParaphaseResults],
    all_split_bams: Dict[str, dict],
    outdir: str,
    clobber: bool,
    reference_path: str,
    gtf_file: Optional[str],
    input_dir: str,
    is_puretarget: bool,
) -> List[RegionEntry]:
    """
    Process one trio: build trio table entries, generate orographer plots.
    Used by tests and kept for compatibility; main pipeline uses workers.
    """
    make_output_dirs(outdir, trio.IndividualID + "-trio", clobber)
    trio_region_entries = make_trio_table_entries(
        trio,
        all_paraphase_results[trio.IndividualID],
        all_paraphase_results[trio.PaternalID],
        all_paraphase_results[trio.MaternalID],
        all_split_bams,
        outdir,
    )
    return generate_orographer_plots(
        trio_region_entries,
        outdir,
        reference_path,
        gtf_file,
        input_dir,
        is_puretarget,
    )


def _run_stage_splits(
    split_args_list: List[Tuple[str, SplitBamArgs]],
    executor: Optional[ProcessPoolExecutor],
    timeout: int,
) -> Dict[str, dict]:
    """Run split_bam for all samples. Returns all_split_bams[sample][region]."""
    all_split_bams = {}
    if executor is None:
        for sample, args in split_args_list:
            all_split_bams[sample] = split_bam_worker(args)
        return all_split_bams

    futures = {
        executor.submit(split_bam_worker, args): sample
        for sample, args in split_args_list
    }
    try:
        for future in as_completed(futures.keys(), timeout=timeout):
            sample = futures[future]
            all_split_bams[sample] = future.result(timeout=timeout)
    except TimeoutError:
        pending_samples = [s for f, s in futures.items() if not f.done()]
        logger.error(
            "BAM split timed out (%s s). Samples that did not complete in time: %s",
            timeout,
            pending_samples,
        )
        for future in futures:
            future.cancel()
        sys.exit(1)
    except Exception:
        for future in futures:
            future.cancel()
        raise
    return all_split_bams


def _run_stage_plots(
    single_tasks: List[SinglePlotTask],
    trio_tasks: List[TrioPlotTask],
    executor: Optional[ProcessPoolExecutor],
    timeout: int,
) -> Tuple[Dict[int, str], Dict[int, str]]:
    """
    Run plot workers. Returns (single_results, trio_results) index -> html_path.
    Raises on first failure.
    """
    single_results = {}
    trio_results = {}

    if executor is None:
        for task in single_tasks:
            index, html_path = generate_single_plot_worker(task)
            single_results[index] = html_path
        for task in trio_tasks:
            index, html_path = generate_trio_plot_worker(task)
            trio_results[index] = html_path
        return single_results, trio_results

    futures = []
    task_by_future = {}
    for task in single_tasks:
        future = executor.submit(generate_single_plot_worker, task)
        futures.append(future)
        task_by_future[future] = ("single", task.index)
    for task in trio_tasks:
        future = executor.submit(generate_trio_plot_worker, task)
        futures.append(future)
        task_by_future[future] = ("trio", task.index)

    try:
        for future in as_completed(futures, timeout=timeout):
            kind, _ = task_by_future[future]
            result_index, html_path = future.result(timeout=timeout)
            if kind == "single":
                single_results[result_index] = html_path
            else:
                trio_results[result_index] = html_path
    except TimeoutError:
        pending = []
        for future in futures:
            if not future.done():
                kind, idx = task_by_future[future]
                if kind == "single":
                    t = single_tasks[idx]
                    pending.append(f"{t.sample}/{t.region}")
                else:
                    t = trio_tasks[idx]
                    pending.append(f"{t.proband_id}-trio/{t.region}")
        logger.error(
            "Plot generation timed out (%s s). Tasks that did not complete in time: %s",
            timeout,
            pending,
        )
        for future in futures:
            future.cancel()
        raise
    except Exception:
        for future in futures:
            future.cancel()
        raise

    return single_results, trio_results


def _cancel_futures_and_shutdown(futures: list, executor: ProcessPoolExecutor) -> None:
    """Cancel all pending futures and shut down the executor without waiting."""
    for future in futures:
        future.cancel()
    executor.shutdown(wait=False)


def paraviewer(
    paraphase_dir: str,
    ptcp_dir: str,
    outdir: str,
    pedigree: str,
    include_only_samples: list,
    exclude_samples: list,
    include_only_regions: list,
    exclude_regions: list,
    max_reads_per_haplotype: int,
    ref: str,
    gtf: str,
    clobber: bool,
    task_timeout_seconds: int,
    threads: int,
) -> int:
    """
    Run the paraviewer pipeline. Returns 0=success, 1=failure, 130=interrupt.
    """
    cpu_count = os.cpu_count()
    if threads > 1 and cpu_count is not None and threads > cpu_count:
        effective_workers = cpu_count
        logger.warning(
            "--threads %s exceeds cpu_count %s; using %s workers",
            threads,
            cpu_count,
            effective_workers,
        )
    else:
        effective_workers = min(threads, cpu_count or 1) if threads > 1 else 1

    base_timeout = task_timeout_seconds

    include_only_regions = include_only_regions or []
    exclude_regions = exclude_regions or []
    include_tuple = tuple(include_only_regions)
    exclude_tuple = tuple(exclude_regions)

    pedigree_dict = read_pedigree_file(pedigree, include_only_samples, exclude_samples)

    all_paraphase_results = {}
    is_puretarget = False
    input_dir = None
    if paraphase_dir:
        all_paraphase_results = get_paraphase_results(
            paraphase_dir,
            include_only_samples,
            exclude_samples,
            pedigree_dict,
        )
        input_dir = paraphase_dir
    elif ptcp_dir:
        all_paraphase_results = get_puretarget_results(
            ptcp_dir,
            include_only_samples,
            exclude_samples,
            pedigree_dict,
        )
        is_puretarget = True
        input_dir = ptcp_dir

    if len(all_paraphase_results) == 0:
        logger.error("No results found in input directory")
        return 1

    trio_samples, trio_member_ids = get_trio_samples(
        pedigree_dict, all_paraphase_results
    )

    for sample in all_paraphase_results:
        make_output_dirs(outdir, sample, clobber)
    for sample in trio_samples:
        make_output_dirs(outdir, sample, clobber)

    split_args_list = []
    for sample in all_paraphase_results:
        results = all_paraphase_results[sample]
        split_args_list.append(
            (
                sample,
                SplitBamArgs(
                    bam_path=results.BAM,
                    bai_path=results.BAI,
                    outdir=outdir,
                    sample=sample,
                    include_only_regions=include_tuple,
                    exclude_regions=exclude_tuple,
                    max_reads_per_hap=max_reads_per_haplotype,
                ),
            )
        )

    executor = None
    if effective_workers > 1:
        executor = ProcessPoolExecutor(max_workers=effective_workers)

    try:
        if executor is None:
            all_split_bams = _run_stage_splits(split_args_list, None, base_timeout)
        else:
            all_split_bams = _run_stage_splits(split_args_list, executor, base_timeout)

        if not all_split_bams:
            logger.error("No split BAMs produced for any sample.")
            return 1

        single_entries = []
        for sample in all_paraphase_results:
            results = all_paraphase_results[sample]
            if sample not in all_split_bams:
                logger.error("No regions found for sample %s.", sample)
                return 1
            sample_entries = make_table_entries(
                results,
                pedigree_dict.get(sample),
                all_split_bams[sample],
                sample in trio_member_ids,
            )
            single_entries.extend(sample_entries)

        output_dir = path.join(outdir, OROGRAPHER_OUTPUT_PATH)
        os.makedirs(output_dir, exist_ok=True)

        single_tasks = []
        for idx, entry in enumerate(single_entries):
            bam_path = (
                path.join(outdir, entry.BAM) if isinstance(entry.BAM, str) else None
            )
            if bam_path is None:
                continue
            vcf_path = find_vcf_file(
                input_dir, entry.Sample, entry.Region, is_puretarget
            )
            prefix = f"{entry.Sample}_{entry.Region}"
            single_tasks.append(
                SinglePlotTask(
                    index=idx,
                    sample=entry.Sample,
                    region=entry.Region,
                    bam_path=bam_path,
                    chrom=entry.Chrom,
                    start=entry.Start,
                    end=entry.End,
                    ref=ref,
                    gtf=gtf,
                    vcf=vcf_path,
                    output_dir=output_dir,
                    prefix=prefix,
                )
            )

        trio_entries = []
        for sample in trio_samples:
            trio = trio_samples[sample]
            trio_region_entries = make_trio_table_entries(
                trio,
                all_paraphase_results[trio.IndividualID],
                all_paraphase_results[trio.PaternalID],
                all_paraphase_results[trio.MaternalID],
                all_split_bams,
                outdir,
            )
            trio_entries.extend(trio_region_entries)

        trio_tasks = []
        for idx, entry in enumerate(trio_entries):
            if not isinstance(entry.BAM, list) or len(entry.BAM) != 3:
                continue
            bam_paths = [
                path.join(outdir, entry.BAM[0]),
                path.join(outdir, entry.BAM[1]),
                path.join(outdir, entry.BAM[2]),
            ]
            proband_id = entry.Sample.replace("-trio", "")
            vcf_paternal = find_vcf_file(
                input_dir, entry.PaternalID, entry.Region, is_puretarget
            )
            vcf_maternal = find_vcf_file(
                input_dir, entry.MaternalID, entry.Region, is_puretarget
            )
            vcf_proband = find_vcf_file(
                input_dir, proband_id, entry.Region, is_puretarget
            )
            vcf_paths = [vcf_paternal, vcf_maternal, vcf_proband]
            prefix = f"{entry.Sample}_{entry.Region}"
            trio_tasks.append(
                TrioPlotTask(
                    index=idx,
                    proband_id=proband_id,
                    paternal_id=entry.PaternalID,
                    maternal_id=entry.MaternalID,
                    region=entry.Region,
                    bam_paths=bam_paths,
                    vcf_paths=vcf_paths,
                    chrom=entry.Chrom,
                    start=entry.Start,
                    end=entry.End,
                    ref=ref,
                    gtf=gtf,
                    output_dir=output_dir,
                    prefix=prefix,
                )
            )

        plot_timeout = TASK_TIMEOUT_SECONDS * 3 if trio_tasks else base_timeout
        single_results, trio_results = _run_stage_plots(
            single_tasks, trio_tasks, executor, plot_timeout
        )

        for idx, html_path in single_results.items():
            single_entries[idx] = single_entries[idx]._replace(OrographerHTML=html_path)
        for idx, html_path in trio_results.items():
            trio_entries[idx] = trio_entries[idx]._replace(OrographerHTML=html_path)

        all_entries = single_entries + trio_entries
        all_entries.sort(key=lambda entry: (entry.Sample, entry.Region))

        build_review_page(outdir, all_entries)

        data_dir = path.join(outdir, "data")
        if path.exists(data_dir):
            shutil.rmtree(data_dir)

        return 0

    except KeyboardInterrupt:
        if executor is not None:
            executor.shutdown(wait=False)
        logger.warning(
            "Exited early (keyboard interrupt); %s partial. Delete or use --clobber.",
            outdir,
        )
        return 130
    except (TimeoutError, Exception) as pipeline_error:
        if executor is not None:
            executor.shutdown(wait=False)
        logger.exception("Pipeline failed: %s; exiting.", pipeline_error)
        sys.exit(1)
    finally:
        if executor is not None:
            executor.shutdown(wait=False)
