#!/usr/bin/env python

from os import path, makedirs
import pathlib
import shutil
import sys
import json
from typing import Optional
from cachetools import LRUCache
import pysam
import yaml
import gzip
import logging
from collections import namedtuple

logger = logging.getLogger(__name__)

REGION_PADDING = 1000

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
        "Image",
        "IGVSession",
        "FamilyID",
        "PaternalID",
        "MaternalID",
        "Sex",
        "Phenotype",
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

IMAGES_PATH = "data/{sample}/images"
IGV_SESSIONS_PATH = "data/{sample}/igv_sessions"
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
        raise ValueError(f"Invalid region format or values: {e}")
    except Exception:
        raise ValueError(
            "Input must be in the format 'chrN:start-end' with numeric coordinates."
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


def get_config(genome_build, source_pipeline):
    """
    Read region config files and return as a dictionary.

    Args:
        genome_build (str): The genome build identifier (e.g., 'hg38', 'hg19')

    Returns:
        dict: The configuration data from the YAML file keyed by region name

    Raises:
        FileNotFoundError: If the config file doesn't exist
        yaml.YAMLError: If the YAML file is malformed
    """
    assert source_pipeline in ["paraphase", "puretarget"]
    data_path = path.join(path.dirname(__file__), "data", genome_build)
    config_path = path.join(data_path, f"config_{source_pipeline}.yaml")

    if not path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    try:
        with open(config_path, "r") as f:
            config_data = yaml.safe_load(f)
        return config_data
    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML file {config_path}: {e}")
        raise


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
        logger.error(
            "For safety reasons, Paraviewer cannot output to root or home directories (`/` or $HOME)"
        )
        sys.exit()

    if path.exists(data_dir) and not clobber:
        logger.error(
            f"Output data directory {data_dir} already exists and --clobber is not set"
        )
        sys.exit(1)

    # Create the data directories for this sample
    for new_path in (
        path.join(outdir, IMAGES_PATH.format(sample=sample)),
        path.join(outdir, IGV_SESSIONS_PATH.format(sample=sample)),
        path.join(outdir, BAMS_PATH.format(sample=sample)),
    ):
        if not path.exists(new_path):
            makedirs(new_path, exist_ok=True)


def copy_trio_bams(
    trio,
    region,
    outdir,
    sample_bam,
    paternal_bam,
    maternal_bam,
):
    bams = []
    new_bam_pattern = new_target_bam = path.join(
        BAMS_PATH.format(sample=trio.IndividualID + "-trio"),
        "{}_{}.bam",
    )
    for sample, bam in (
        (trio.PaternalID, paternal_bam),
        (trio.MaternalID, maternal_bam),
        (trio.IndividualID, sample_bam),
    ):
        bam = path.join(outdir, bam)
        bai = bam + ".bai"
        new_target_bam = path.join(new_bam_pattern.format(sample, region))
        new_target_bai = new_target_bam + ".bai"

        shutil.copy(bam, path.join(outdir, new_target_bam))
        shutil.copy(bai, path.join(outdir, new_target_bai))

        bams.append(new_target_bam)

    return bams


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
    Split a BAM file into smaller chunks based on RN tags for improved IGV visualization.
    Each haplotype is allowed at most max_reads_per_hap reads.

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
    region_bam_paths = namedtuple("RegionBamPaths", ["BAM", "BAI"])
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
                    "{} haplotype {}: {} reads".format(
                        sample, hap, haplotype_read_counts[hap]
                    )
                )

    finally:
        for file_handle in region_files_cache.values():
            if file_handle is not None:
                file_handle.close()

        for region_name, bam_path in result_abs_paths.items():
            pysam.index(bam_path)
            result_relative_paths[region_name] = region_bam_paths(
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
