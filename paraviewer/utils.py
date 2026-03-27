#!/usr/bin/env python

import gzip
import json
import logging
import os
import pathlib
import sys
from collections import namedtuple
from os import makedirs, path
from typing import Optional, Set

import pysam
from cachetools import LRUCache

logger = logging.getLogger(__name__)

# Fractional padding applied to realign regions (e.g., 0.10 == 10% of the target region)
REGION_PADDING = 0.10

# Define all namedtuples at module level
RegionEntry = namedtuple(
    "RegionEntry",
    [
        "Chrom",
        "Start",
        "End",
        "Region",
        "Sample",
        "BAM",
        "BAI",
        "CopyNumber",
        "SpecialInfo",
        "OrographerHTML",
        "FamilyID",
        "PaternalID",
        "MaternalID",
        "Sex",
        "Phenotype",
        "IsTrioSample",
    ],
)
PedigreeEntry = namedtuple(
    "PedigreeEntry",
    ["FamilyID", "IndividualID", "PaternalID", "MaternalID", "Sex", "Phenotype"],
)
ParaphaseResults = namedtuple(
    "ParaphaseResults", ["Sample", "BAM", "BAI", "JSON", "F8_INV", "HAVANNO"]
)
HAVANNO_INFO = namedtuple(
    "HAVANNO_INFO", ["Haplotype", "PathogenicVariants", "Insertion", "Deletion"]
)
GenomicInterval = namedtuple("GenomicInterval", ["Chrom", "Start", "End"])
RegionBamPaths = namedtuple("RegionBamPaths", ["BAM", "BAI"])

OROGRAPHER_OUTPUT_PATH = "orographer_output"
BAMS_PATH = "data/{sample}/bams"


def genomic_interval_from_str(region_str):
    """
    Given a string with a coordinate, return the coordinate as a GenomicInterval
    """
    try:
        chrom_part, pos_part = region_str.strip().split(":")
        start_str, end_str = pos_part.split("-")
        start = int(start_str.replace(",", ""))
        end = int(end_str.replace(",", ""))

        if start < 0 or end < 0:
            raise ValueError("Coordinates must be non-negative.")

        return GenomicInterval(chrom_part, start, end)

    except ValueError as e:
        raise ValueError(f"Invalid region format or values: {e}") from e
    except Exception as e:
        raise ValueError(
            "Input must be in the format 'chrN:start-end' with numeric coordinates."
        ) from e


def parse_phase_region(phase_region: str) -> "GenomicInterval":
    """
    Parse Paraphase phase_region strings, e.g. ``38:chr6:32013300-32046200``
    (genome build prefix, chromosome, start-end).
    """
    s = phase_region.strip()
    parts = s.split(":")
    if len(parts) < 3:
        raise ValueError(
            f"phase_region must look like '38:chr6:start-end', got: {phase_region!r}"
        )
    pos_part = parts[-1]
    chrom = parts[-2]
    try:
        start_str, end_str = pos_part.split("-", 1)
        start = int(start_str.replace(",", ""))
        end = int(end_str.replace(",", ""))
    except ValueError as e:
        raise ValueError(f"Invalid phase_region coordinates: {phase_region!r}") from e
    if start < 0 or end < 0:
        raise ValueError("Coordinates must be non-negative.")
    return GenomicInterval(chrom, start, end)


OLD_PARAPHASE_NO_PHASE_REGION = (
    "Paraphase JSON for region %r has no 'phase_region' field. "
    "This output is from an older Paraphase version; rerun with paraphase >=v3.3.0"
)


def is_gzipped(putative_zipfile):
    """
    Check if file is zipped
    """
    with open(putative_zipfile, "rb") as filehandle:
        id_bytes = filehandle.read(2)
        return id_bytes == b"\x1f\x8b"


def unpack_json(json_filename):
    """
    unpacks a json or json.gz file into a dict and returns it
    """
    if not (path.exists(json_filename) and path.isfile(json_filename)):
        logger.warning(" {} does not exist".format(json_filename))
        return
    if json_filename.endswith(".gz"):
        if not is_gzipped(json_filename):
            logger.warning(
                "{} is identified as gzipped but is not".format(json_filename)
            )
            return
        with gzip.open(json_filename, "rt", encoding="UTF-8") as json_fh:
            try:
                return json.load(json_fh)
            except json.decoder.JSONDecodeError:
                logger.warning(" {} is empty or misformatted".format(json_filename))
                return
    else:
        with open(json_filename, "r") as json_fh:
            try:
                return json.load(json_fh)
            except json.decoder.JSONDecodeError:
                logger.warning(" {} is empty or misformatted".format(json_filename))
                return


def is_mac():
    return sys.platform == "darwin"


def is_linux():
    return sys.platform.startswith("linux")


def parse_sample_name_from_paraphase_output(file_path):
    """
    Parse a sample name from a paraphase output file
    """
    base_name = path.basename(file_path)
    file_name, ext = path.splitext(base_name)
    while ext and ext != ".paraphase":
        file_name, ext = path.splitext(file_name)
    return file_name


def make_output_dirs(
    outdir: str,
    sample: str,
    clobber,
):
    """
    Create the expected output directory structure.
    If already exists and clobber is not set, exit.
    """
    data_dir = path.join(outdir, "data")
    if path == pathlib.Path("/") or path == pathlib.Path.home():
        logger.error("Paraviewer cannot write to root or home (safety).")
        sys.exit()

    if path.exists(data_dir) and not clobber:
        logger.error("Output data dir %s exists; use --clobber to overwrite.", data_dir)
        sys.exit(1)

    # Create the data directories for this sample
    for new_path in (path.join(outdir, BAMS_PATH.format(sample=sample)),):
        if not path.exists(new_path):
            makedirs(new_path, exist_ok=True)

    # Create orographer output directory (shared across all samples)
    orographer_output_dir = path.join(outdir, OROGRAPHER_OUTPUT_PATH)
    if not path.exists(orographer_output_dir):
        makedirs(orographer_output_dir, exist_ok=True)


def split_bam(
    bam_path: str,
    bai_path: str,
    outdir: str,
    sample: str,
    include_only_regions: Optional[str],
    exclude_regions: Optional[str],
    max_reads_per_hap: int,
) -> dict:
    """
    Split BAM into chunks by RN tags for Orographer.
    Max max_reads_per_hap per haplotype.

    Args:
        bam_path: Path to the input BAM file
        bai_path: Path to the input BAI file
        outdir: Output directory path
        sample: Sample ID

    Returns:
        Dictionary mapping region names to namedtuples containing BAM and BAI paths
    """
    # Create an LRU cache for file handles with max size of 5
    region_files_cache = LRUCache(maxsize=5)
    result_abs_paths = {}
    result_relative_paths = {}
    new_bam_pattern = path.join(BAMS_PATH.format(sample=sample), "{}_{}.bam")

    try:
        with pysam.AlignmentFile(bam_path, "rb") as in_bam:
            # Create a new simplified header
            new_header = pysam.AlignmentHeader.from_dict(
                {"HD": in_bam.header["HD"], "SQ": in_bam.header["SQ"]}
            )
            haplotype_read_counts = {}
            for read in in_bam:
                if not read.has_tag("RN"):
                    continue

                region_name = read.get_tag("RN")
                if (
                    include_only_regions
                    and region_name.lower() not in include_only_regions
                ):
                    continue
                if exclude_regions and region_name.lower() in exclude_regions:
                    continue

                hp = "unknown"
                if read.has_tag("HP"):
                    hp = read.get_tag("HP")

                if hp and hp not in haplotype_read_counts:
                    haplotype_read_counts[hp] = 0
                if hp and haplotype_read_counts[hp] >= max_reads_per_hap:
                    continue

                # If we haven't created a file for this region yet, create it
                if region_name not in region_files_cache:
                    bam_out_path = path.join(
                        outdir,
                        new_bam_pattern.format(sample, region_name),
                    )
                    region_files_cache[region_name] = pysam.AlignmentFile(
                        bam_out_path, "wb", header=new_header, index_filename=bai_path
                    )
                    result_abs_paths[region_name] = bam_out_path

                region_files_cache[region_name].write(read)
                if hp:
                    haplotype_read_counts[hp] += 1
            for hap in haplotype_read_counts:
                logger.debug(
                    "%s haplotype %s: %s reads",
                    sample,
                    hap,
                    haplotype_read_counts[hap],
                )

    finally:
        for file_handle in region_files_cache.values():
            if file_handle is not None:
                file_handle.close()

        for region_name, bam_path in result_abs_paths.items():
            pysam.index(bam_path)
            result_relative_paths[region_name] = RegionBamPaths(
                BAM=new_bam_pattern.format(sample, region_name),
                BAI=new_bam_pattern.format(sample, region_name) + ".bai",
            )
    return result_relative_paths


def strip_suffix_from_path(path_str: str, suffix: str) -> str:
    """
    Given a path as a string, remove a suffix from it and send it back
    """

    path_obj = pathlib.Path(path_str).expanduser()
    if path_obj.name.endswith(suffix):
        new_name = path_obj.name[: -len(suffix)]
        return str(path_obj.with_name(new_name))
    return path_str


def find_vcf_file(
    paraphase_dir: str,
    sample: str,
    region: str,
    is_puretarget: bool = False,
) -> Optional[str]:
    """
    Find VCF file for a given sample and region in paraphase directory structure.

    Args:
        paraphase_dir: Path to paraphase directory (or ptcp directory for puretarget)
        sample: Sample name
        region: Region name
        is_puretarget: If True, look in {sample}_paraphase subdirectory

    Returns:
        Path to VCF file if found, None otherwise
    """
    if is_puretarget:
        # Puretarget: ptcp_dir/{sample}_paraphase/{sample}_paraphase_vcfs/
        vcf_dir = path.join(
            paraphase_dir, f"{sample}_paraphase", f"{sample}_paraphase_vcfs"
        )
    else:
        # For paraphase: {paraphase_dir}/{sample}_paraphase_vcfs/{sample}_{region}.vcf
        vcf_dir = path.join(paraphase_dir, f"{sample}_paraphase_vcfs")

    # Try both .vcf and .vcf.gz extensions
    vcf_candidates = [
        path.join(vcf_dir, f"{sample}_{region}.vcf"),
        path.join(vcf_dir, f"{sample}_{region}.vcf.gz"),
    ]

    for vcf_path in vcf_candidates:
        if path.exists(vcf_path) and path.isfile(vcf_path):
            return vcf_path

    return None


def delete_bam(bam_path: str):
    """
    Delete BAM and its paired BAI index (if present) to reduce disk size.

    Args:
        bam_path: Absolute BAM or BAI file path to delete (pair inferred)
    """
    candidates: Set[str] = {bam_path}
    if bam_path.endswith(".bam"):
        candidates.add(bam_path + ".bai")
    elif bam_path.endswith(".bai"):
        # delete paired BAM if looks like *.bam.bai or generic .bai
        if bam_path.endswith(".bam.bai"):
            candidates.add(bam_path[:-4])
        else:
            candidates.add(bam_path[:-4])

    for filepath in candidates:
        try:
            if path.exists(filepath) and path.isfile(filepath):
                os.remove(filepath)
        except OSError as os_error:
            logger.warning(f"Failed to delete {filepath}: {os_error}")
