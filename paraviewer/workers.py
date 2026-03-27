"""
Top-level, stateless worker functions for multiprocessing.
All functions and task types are picklable (no closures; plain data).
"""

from __future__ import print_function

import logging
from collections import namedtuple
from os import path

from orographer.orographer import orographer
from orographer.utils import OutputConfig

from paraviewer.utils import (
    OROGRAPHER_OUTPUT_PATH,
    split_bam,
)

logger = logging.getLogger(__name__)

PARAPHASE_REGION_TYPE = "paraphase"

# Picklable task DTOs: plain strings/ints/bools or lists of strings
SplitBamArgs = namedtuple(
    "SplitBamArgs",
    [
        "bam_path",
        "bai_path",
        "outdir",
        "sample",
        "include_only_regions",
        "exclude_regions",
        "max_reads_per_hap",
    ],
)

SinglePlotTask = namedtuple(
    "SinglePlotTask",
    [
        "index",
        "sample",
        "region",
        "bam_path",
        "chrom",
        "start",
        "end",
        "ref",
        "gtf",
        "vcf",
        "output_dir",
        "prefix",
    ],
)

TrioPlotTask = namedtuple(
    "TrioPlotTask",
    [
        "index",
        "proband_id",
        "paternal_id",
        "maternal_id",
        "region",
        "bam_paths",
        "vcf_paths",
        "chrom",
        "start",
        "end",
        "ref",
        "gtf",
        "output_dir",
        "prefix",
    ],
)


def split_bam_worker(args: SplitBamArgs):
    """
    Run split_bam for one sample. Returns dict[region_name] -> RegionBamPaths(BAM, BAI).
    Raises on failure (all-or-nothing).
    """
    include_regions = (
        list(args.include_only_regions) if args.include_only_regions else None
    )
    exclude_regions = list(args.exclude_regions) if args.exclude_regions else None
    result = split_bam(
        args.bam_path,
        args.bai_path,
        args.outdir,
        args.sample,
        include_regions,
        exclude_regions,
        args.max_reads_per_hap,
    )
    if not result:
        raise ValueError(f"No regions produced for sample {args.sample}")
    return result


def generate_single_plot_worker(task: SinglePlotTask):
    """
    Run orographer for one single-sample region. Returns (task.index, html_path_rel).
    Raises on failure or if output file is missing.
    """
    coordinate_str = f"{task.chrom}:{task.start}-{task.end}"
    output_config = OutputConfig(task.output_dir, task.prefix)
    orographer(
        region_type=PARAPHASE_REGION_TYPE,
        bam_file=task.bam_path,
        coordinate_strs=[coordinate_str],
        reference_path=task.ref,
        output_config=output_config,
        gtf_file=task.gtf,
        vcf_file=task.vcf,
        other_bam_files=None,
        other_vcf_files=None,
        sample_label=task.sample,
        other_sample_labels=None,
    )
    html_filename = f"{task.prefix}_{task.chrom}_{task.start}_{task.end}_bokeh.html"
    html_path_rel = path.join(OROGRAPHER_OUTPUT_PATH, html_filename)
    html_full = path.join(path.dirname(task.output_dir), html_path_rel)
    if not path.exists(html_full):
        raise FileNotFoundError(f"Orographer did not create {html_full}")
    return (task.index, html_path_rel)


def generate_trio_plot_worker(task: TrioPlotTask):
    """
    Run orographer for one trio (3 BAMs). Returns (task.index, html_path_rel).
    Raises on failure or if output file is missing.
    """
    coordinate_str = f"{task.chrom}:{task.start}-{task.end}"
    proband_bam = task.bam_paths[2]
    other_bam_files = [task.bam_paths[0], task.bam_paths[1]]
    other_vcf_files = (
        [task.vcf_paths[0], task.vcf_paths[1]] if task.vcf_paths else [None, None]
    )
    output_config = OutputConfig(task.output_dir, task.prefix)
    sample_label = f"Proband ({task.proband_id})"
    other_sample_labels = [
        f"Paternal ({task.paternal_id})",
        f"Maternal ({task.maternal_id})",
    ]
    orographer(
        region_type=PARAPHASE_REGION_TYPE,
        bam_file=proband_bam,
        coordinate_strs=[coordinate_str],
        reference_path=task.ref,
        output_config=output_config,
        gtf_file=task.gtf,
        vcf_file=(
            task.vcf_paths[2] if task.vcf_paths and len(task.vcf_paths) > 2 else None
        ),
        other_bam_files=other_bam_files,
        other_vcf_files=other_vcf_files,
        sample_label=sample_label,
        other_sample_labels=other_sample_labels,
    )
    coord_part = coordinate_str.replace(":", "_").replace("-", "_")
    html_filename = f"{task.prefix}_{coord_part}_bokeh.html"
    html_path_rel = path.join(OROGRAPHER_OUTPUT_PATH, html_filename)
    outdir_root = path.dirname(task.output_dir)
    html_full = path.join(outdir_root, html_path_rel)
    if not path.exists(html_full):
        raise FileNotFoundError(f"Orographer did not create {html_full}")
    return (task.index, html_path_rel)
