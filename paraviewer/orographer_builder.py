#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generate orographer HTML/JSON plots for display in web viewer
"""

from __future__ import print_function

import logging
import os
from os import path
from typing import List, Optional

from orographer.orographer import orographer
from orographer.utils import OutputConfig

from paraviewer.utils import OROGRAPHER_OUTPUT_PATH, RegionEntry, find_vcf_file

logger = logging.getLogger(__name__)

PARAPHASE_REGION_TYPE = "paraphase"


def generate_trio_plots_for_entry(
    region_entry: RegionEntry,
    outdir: str,
    reference_path: str,
    gtf_file: Optional[str],
    input_dir: str,
    is_puretarget: bool,
) -> Optional[RegionEntry]:
    """
    Generate a single orographer plot with three BAMs (paternal, maternal, proband).
    Orographer plots them in one HTML (top to bottom: paternal, maternal, proband).
    """
    coordinate_str = f"{region_entry.Chrom}:{region_entry.Start}-{region_entry.End}"
    proband_id = region_entry.Sample.replace("-trio", "")
    paternal_id = region_entry.PaternalID
    maternal_id = region_entry.MaternalID

    # region_entry.BAM is [paternal_rel, maternal_rel, proband_rel]
    paternal_bam = path.join(outdir, region_entry.BAM[0])
    maternal_bam = path.join(outdir, region_entry.BAM[1])
    proband_bam = path.join(outdir, region_entry.BAM[2])

    for bam_path in (paternal_bam, maternal_bam, proband_bam):
        if not path.exists(bam_path):
            logger.warning(
                f"BAM file not found for trio: {bam_path}, "
                f"skipping {region_entry.Sample} {region_entry.Region}"
            )
            return None

    output_dir = path.join(outdir, OROGRAPHER_OUTPUT_PATH)
    prefix = f"{region_entry.Sample}_{region_entry.Region}"
    output_config = OutputConfig(output_dir, prefix)

    proband_vcf = find_vcf_file(
        input_dir, proband_id, region_entry.Region, is_puretarget
    )
    paternal_vcf = find_vcf_file(
        input_dir, paternal_id, region_entry.Region, is_puretarget
    )
    maternal_vcf = find_vcf_file(
        input_dir, maternal_id, region_entry.Region, is_puretarget
    )
    other_bam_files = [paternal_bam, maternal_bam]
    other_vcf_files = [paternal_vcf, maternal_vcf]
    # Primary BAM = proband; others = paternal, maternal (top, middle, bottom)
    # Labels shown in plot: "Paternal (id)", "Maternal (id)", "Proband (id)"
    sample_label = f"Proband ({proband_id})"
    other_sample_labels = [f"Paternal ({paternal_id})", f"Maternal ({maternal_id})"]

    orographer(
        region_type=PARAPHASE_REGION_TYPE,
        bam_file=proband_bam,
        coordinate_strs=[coordinate_str],
        reference_path=reference_path,
        output_config=output_config,
        gtf_file=gtf_file,
        vcf_file=proband_vcf,
        other_bam_files=other_bam_files,
        other_vcf_files=other_vcf_files,
        sample_label=sample_label,
        other_sample_labels=other_sample_labels,
    )

    # Filename: prefix_coordstr_bokeh.html (":" and "-" → "_")
    coord_part = coordinate_str.replace(":", "_").replace("-", "_")
    html_filename = f"{prefix}_{coord_part}_bokeh.html"
    html_path_rel = path.join(OROGRAPHER_OUTPUT_PATH, html_filename)
    html_full_path = path.join(outdir, html_path_rel)
    if not path.exists(html_full_path):
        logger.warning(
            "Orographer HTML not created for trio: %s, skipping %s %s",
            html_full_path,
            region_entry.Sample,
            region_entry.Region,
        )
        return None

    return region_entry._replace(OrographerHTML=html_path_rel)


def generate_single_plot_for_entry(
    region_entry: RegionEntry,
    outdir: str,
    reference_path: str,
    gtf_file: Optional[str],
    input_dir: str,
    is_puretarget: bool,
) -> Optional[RegionEntry]:
    coordinate_str = f"{region_entry.Chrom}:{region_entry.Start}-{region_entry.End}"
    bam_file = path.join(outdir, region_entry.BAM)
    if not path.exists(bam_file):
        logger.warning(
            "BAM not found: %s, skipping %s %s",
            bam_file,
            region_entry.Sample,
            region_entry.Region,
        )
        return None

    prefix = f"{region_entry.Sample}_{region_entry.Region}"
    output_dir = path.join(outdir, OROGRAPHER_OUTPUT_PATH)
    output_config = OutputConfig(output_dir, prefix)
    vcf_file = find_vcf_file(
        input_dir, region_entry.Sample, region_entry.Region, is_puretarget
    )

    orographer(
        region_type=PARAPHASE_REGION_TYPE,
        bam_file=bam_file,
        coordinate_strs=[coordinate_str],
        reference_path=reference_path,
        output_config=output_config,
        gtf_file=gtf_file,
        vcf_file=vcf_file,
        other_bam_files=None,
        other_vcf_files=None,
        sample_label=region_entry.Sample,
        other_sample_labels=None,
    )

    html_filename = (
        f"{prefix}_{region_entry.Chrom}_{region_entry.Start}_"
        f"{region_entry.End}_bokeh.html"
    )
    html_path_rel = path.join(OROGRAPHER_OUTPUT_PATH, html_filename)
    html_full_path = path.join(outdir, html_path_rel)
    if not path.exists(html_full_path):
        logger.warning(
            "Orographer HTML not created: %s, skipping %s %s",
            html_full_path,
            region_entry.Sample,
            region_entry.Region,
        )
        return None

    return region_entry._replace(OrographerHTML=html_path_rel)


def generate_orographer_plots(
    sample_region_entries: List[RegionEntry],
    outdir: str,
    reference_path: str,
    gtf_file: Optional[str],
    input_dir: str,
    is_puretarget: bool,
) -> List[RegionEntry]:
    """
    Generate orographer HTML/JSON plots for each region entry.

    Args:
        sample_region_entries: List of RegionEntry namedtuples
        outdir: Output directory path
        reference_path: Path to reference FASTA file (required)
        gtf_file: Path to GTF/GFF3 file (optional, can be None)
        input_dir: Paraphase/puretarget input dir (for VCF lookup)
        is_puretarget: True if puretarget (affects VCF path)

    Returns:
        List of RegionEntry namedtuples with OrographerHTML paths populated.
        Failed regions are excluded from the returned list.
    """
    if len(sample_region_entries) == 0:
        logger.warning("No valid regions for orographer plot generation")
        return []

    orographer_output_dir = path.join(outdir, OROGRAPHER_OUTPUT_PATH)
    os.makedirs(orographer_output_dir, exist_ok=True)

    successful_entries: List[RegionEntry] = []

    for region_entry in sample_region_entries:
        try:
            if isinstance(region_entry.BAM, list):
                if len(region_entry.BAM) != 3:
                    logger.warning(
                        "Expected 3 BAMs for trio, got %s, skipping %s %s",
                        len(region_entry.BAM),
                        region_entry.Sample,
                        region_entry.Region,
                    )
                    continue
                updated = generate_trio_plots_for_entry(
                    region_entry,
                    outdir,
                    reference_path,
                    gtf_file,
                    input_dir,
                    is_puretarget,
                )
                if updated:
                    successful_entries.append(updated)
            else:
                updated = generate_single_plot_for_entry(
                    region_entry,
                    outdir,
                    reference_path,
                    gtf_file,
                    input_dir,
                    is_puretarget,
                )
                if updated:
                    successful_entries.append(updated)

        except ValueError as value_error:
            logger.warning(
                "ValueError for %s %s: %s",
                region_entry.Sample,
                region_entry.Region,
                value_error,
            )
            continue
        except FileNotFoundError as file_error:
            logger.warning(
                "FileNotFoundError for %s %s: %s",
                region_entry.Sample,
                region_entry.Region,
                file_error,
            )
            continue
        except OSError as os_error:
            logger.warning(
                "OSError for %s %s: %s",
                region_entry.Sample,
                region_entry.Region,
                os_error,
            )
            continue

    return successful_entries
