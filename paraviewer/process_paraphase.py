import logging
import sys
from glob import glob
from os import path

from paraviewer.special_info import get_special_info
from paraviewer.utils import (
    OLD_PARAPHASE_NO_PHASE_REGION,
    OROGRAPHER_OUTPUT_PATH,
    REGION_PADDING,
    ParaphaseResults,
    PedigreeEntry,
    RegionBamPaths,
    RegionEntry,
    normalize_region_data,
    parse_phase_region,
    parse_sample_name_from_paraphase_output,
    unpack_json,
    warn_if_unsupported_paraphase_version,
)

logger = logging.getLogger(__name__)


def get_paraphase_results(
    paraphase_dir: str,
    include_only_samples: list[str],
    exclude_samples: list[str],
    pedigree_dict: dict[str, PedigreeEntry],
) -> dict[str, ParaphaseResults] | None:
    """
    Validates that expected result files are where
    they should be and returns their paths. For each included sample, should find:
    * BAM
    * BAI
    * JSON(.GZ)

    Returns list[ParaphaseResults]
    """
    all_results = {}
    # check JSON file
    json_matches = glob(path.join(paraphase_dir, "*paraphase.json")) + glob(
        path.join(paraphase_dir, "*paraphase.json.gz")
    )
    if not json_matches or len(json_matches) < 1:
        logger.warning(f"No JSON result file found in {paraphase_dir}")
        return
    for json_filename in json_matches:
        sample = parse_sample_name_from_paraphase_output(json_filename)
        if len(include_only_samples) > 0 and sample.lower() not in include_only_samples:
            continue
        if sample.lower() in exclude_samples:
            continue
        if len(pedigree_dict) > 0 and sample not in pedigree_dict:
            continue

        # check BAM file
        bam_name = path.join(paraphase_dir, f"{sample}.paraphase.bam")
        if not path.isfile(bam_name):
            logger.warning(f"No BAM result file found in {paraphase_dir}")
            continue

        # check BAI file
        bai_name = path.join(paraphase_dir, f"{sample}.paraphase.bam.bai")
        if not path.isfile(bai_name):
            logger.warning(f"No BAM index file found in {paraphase_dir}")
            continue

        warn_if_unsupported_paraphase_version(bam_name, sample)

        all_results[sample] = ParaphaseResults(
            Sample=sample,
            BAI=bai_name,
            BAM=bam_name,
            JSON=json_filename,
            F8_INV="",
            HAVANNO="",
        )
    if len(all_results) == 0:
        logger.debug(f"No samples found for {paraphase_dir}.")

    return all_results


def make_trio_table_entries(
    trio: PedigreeEntry,
    proband_paraphase_results: ParaphaseResults,
    paternal_paraphase_results: ParaphaseResults,
    maternal_paraphase_results: ParaphaseResults,
    all_split_bams: dict[str, dict[str, RegionBamPaths]],
    outdir: str,
) -> list[RegionEntry]:
    """
    Reads the info that will be used for page building from trio paraphase jsons
    and creates namedtuples for the table rows.

    Returns:
        List of RegionEntry namedtuples
    """
    if (
        trio.IndividualID not in all_split_bams
        or trio.PaternalID not in all_split_bams
        or trio.MaternalID not in all_split_bams
    ):
        return []
    proband_paraphase_json_calls = unpack_json(proband_paraphase_results.JSON)

    trio_entries = []
    for region in proband_paraphase_json_calls:
        proband_region_data = normalize_region_data(proband_paraphase_json_calls[region])

        if (
            region not in all_split_bams[trio.IndividualID]
            or region not in all_split_bams[trio.MaternalID]
            or region not in all_split_bams[trio.PaternalID]
        ):
            continue

        # Use split BAM paths (paternal, maternal, proband) for trio plot
        paternal_paths = all_split_bams[trio.PaternalID][region]
        maternal_paths = all_split_bams[trio.MaternalID][region]
        proband_paths = all_split_bams[trio.IndividualID][region]
        bam_paths = [paternal_paths.BAM, maternal_paths.BAM, proband_paths.BAM]
        bai_paths = [paternal_paths.BAI, maternal_paths.BAI, proband_paths.BAI]

        phase_region = (
            proband_region_data.get("phase_region")
            if isinstance(proband_region_data, dict)
            else None
        )
        if not phase_region:
            logger.error(OLD_PARAPHASE_NO_PHASE_REGION, region)
            sys.exit(1)
        try:
            realign_region = parse_phase_region(str(phase_region))
        except ValueError as e:
            logger.error("Invalid phase_region for %r: %s", region, e)
            sys.exit(1)

        total_cn, special_info = get_special_info(
            region, proband_region_data, proband_paraphase_results
        )

        # Compute dynamic padding as 5% of region length
        region_len = max(0, realign_region.End - realign_region.Start)
        pad = max(0, int(region_len * REGION_PADDING))
        padded_start = max(0, realign_region.Start - pad)
        padded_end = max(0, realign_region.End + pad)

        # Single HTML path for native trio plot (populated after plot generation)
        sample_name = trio.IndividualID + "-trio"
        prefix = f"{sample_name}_{region}"
        orographer_html_path = path.join(
            OROGRAPHER_OUTPUT_PATH,
            f"{prefix}_{realign_region.Chrom}_{padded_start}_{padded_end}_bokeh.html",
        )

        trio_entries.append(
            RegionEntry(
                realign_region.Chrom,
                padded_start,
                padded_end,
                region,
                trio.IndividualID + "-trio",
                bam_paths,
                bai_paths,
                total_cn,
                special_info,
                orographer_html_path,
                trio.FamilyID,
                trio.PaternalID,
                trio.MaternalID,
                trio.Sex,
                trio.Phenotype,
                True,
            )
        )
    return trio_entries


def make_table_entries(
    paraphase_results: ParaphaseResults,
    pedigree_entry: PedigreeEntry | None,
    split_bams: dict[str, RegionBamPaths],
    is_trio_sample: bool,
) -> list[RegionEntry]:
    """
    Reads the info that will be used for page building from a json file
    and creates namedtuples for the table rows.

    Args:
        ParaphaseResults namedtuple
        pedigree_entry: Optional PedigreeEntry for this sample

    Returns:
        List of RegionEntry namedtuples
    """
    paraphase_json_calls = unpack_json(paraphase_results.JSON)

    sample_entries = []
    for region in paraphase_json_calls:
        region_data = normalize_region_data(paraphase_json_calls[region])
        if region not in split_bams:
            continue
        bam_path = split_bams[region].BAM
        bai_path = split_bams[region].BAI

        pr = region_data.get("phase_region") if isinstance(region_data, dict) else None
        if not pr:
            logger.error(OLD_PARAPHASE_NO_PHASE_REGION, region)
            sys.exit(1)
        try:
            realign_region = parse_phase_region(str(pr))
        except ValueError as e:
            logger.error("Invalid phase_region for %r: %s", region, e)
            sys.exit(1)

        total_cn, special_info = get_special_info(region, region_data, paraphase_results)

        # Compute dynamic padding as 5% of region length
        region_len = max(0, realign_region.End - realign_region.Start)
        pad = max(0, int(region_len * REGION_PADDING))
        padded_start = max(0, realign_region.Start - pad)
        padded_end = max(0, realign_region.End + pad)

        # Generate orographer HTML path (will be populated after plot generation)
        prefix = f"{paraphase_results.Sample}_{region}"
        orographer_html_path = path.join(
            OROGRAPHER_OUTPUT_PATH,
            f"{prefix}_{realign_region.Chrom}_{padded_start}_{padded_end}_bokeh.html",
        )

        sample_entries.append(
            RegionEntry(
                realign_region.Chrom,
                padded_start,
                padded_end,
                region,
                paraphase_results.Sample,
                bam_path,
                bai_path,
                total_cn,
                special_info,
                orographer_html_path,
                pedigree_entry.FamilyID if pedigree_entry else "",
                pedigree_entry.PaternalID if pedigree_entry else "",
                pedigree_entry.MaternalID if pedigree_entry else "",
                pedigree_entry.Sex if pedigree_entry else "",
                pedigree_entry.Phenotype if pedigree_entry else "",
                is_trio_sample,
            )
        )
    return sample_entries
