from glob import glob
import logging
from os import path
from typing import Dict, Optional
from paraviewer.special_info import get_special_info
from paraviewer.utils import (
    IGV_SESSIONS_PATH,
    IMAGES_PATH,
    REGION_PADDING,
    ParaphaseResults,
    PedigreeEntry,
    RegionEntry,
    copy_trio_bams,
    genomic_interval_from_str,
    parse_sample_name_from_paraphase_output,
    unpack_json,
)

logger = logging.getLogger(__name__)


def get_paraphase_results(
    paraphase_dir: str,
    include_only_samples: list[str],
    exclude_samples: list[str],
    pedigree_dict: Dict[str, PedigreeEntry],
):
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
        logger.warning("No JSON result file found in {}".format(paraphase_dir))
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
        bam_name = path.join(paraphase_dir, "{}.paraphase.bam".format(sample))
        if not path.isfile(bam_name):
            logger.warning("No BAM result file found in {}".format(paraphase_dir))
            continue

        # check BAI file
        bai_name = path.join(paraphase_dir, "{}.paraphase.bam.bai".format(sample))
        if not path.isfile(bai_name):
            logger.warning("No BAM index file found in {}".format(paraphase_dir))
            continue

        all_results[sample] = ParaphaseResults(
            Sample=sample,
            BAI=bai_name,
            BAM=bam_name,
            JSON=json_filename,
            F8_INV="",
            HAVANNO="",
        )
    if len(all_results) == 0:
        logger.warning(f"No samples found for {paraphase_dir}.")

    return all_results


def make_trio_table_entries(
    trio: PedigreeEntry,
    proband_paraphase_results: ParaphaseResults,
    paternal_paraphase_results: ParaphaseResults,
    maternal_paraphase_results: ParaphaseResults,
    all_split_bams: Dict[str, str],
    paraphase_config: Dict[str, Dict],
    outdir: str,
):
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
    paternal_paraphase_json_calls = unpack_json(paternal_paraphase_results.JSON)
    maternal_paraphase_json_calls = unpack_json(maternal_paraphase_results.JSON)

    trio_entries = []
    for region in proband_paraphase_json_calls:
        proband_region_data = proband_paraphase_json_calls[region]
        paternal_region_data = paternal_paraphase_json_calls[region]
        maternal_region_data = maternal_paraphase_json_calls[region]

        if (
            region not in all_split_bams[trio.IndividualID]
            or region not in all_split_bams[trio.MaternalID]
            or region not in all_split_bams[trio.PaternalID]
        ):
            continue

        bam_paths = copy_trio_bams(
            trio,
            region,
            outdir,
            all_split_bams[trio.IndividualID][region].BAM,
            all_split_bams[trio.PaternalID][region].BAM,
            all_split_bams[trio.MaternalID][region].BAM,
        )

        try:
            realign_region = genomic_interval_from_str(
                paraphase_config[region]["realign_region"]
            )
        except Exception as e:
            logger.info(
                f"Failed to find config info for region f{region}, skipping. \n{e}"
            )
            continue

        total_cn, special_info = get_special_info(
            region, proband_region_data, proband_paraphase_results
        )

        igv_session_path = path.join(
            IGV_SESSIONS_PATH.format(sample=trio.IndividualID + "-trio"),
            f"{region}_igv.xml",
        )
        image_path = path.join(
            IMAGES_PATH.format(sample=trio.IndividualID + "-trio"),
            f"{trio.IndividualID + "-trio"}_{region}.png",
        )

        trio_entries.append(
            RegionEntry(
                realign_region.Chrom,
                max(0, realign_region.Start - REGION_PADDING),
                max(0, realign_region.End + REGION_PADDING),
                region,
                trio.IndividualID + "-trio",
                bam_paths,
                [x + ".bai" for x in bam_paths],
                total_cn,
                special_info,
                image_path,
                igv_session_path,
                trio.FamilyID,
                trio.PaternalID,
                trio.MaternalID,
                trio.Sex,
                trio.Phenotype,
            )
        )
    return trio_entries


def make_table_entries(
    paraphase_results: ParaphaseResults,
    pedigree_entry: Optional[PedigreeEntry],
    split_bams: Dict[str, str],
    paraphase_config: Dict[str, Dict],
):
    """
    Reads the info that will be used for page building from a json file
    and creates namedtuples for the table rows.

    Args:
        ParaphaseResults namedtuple
        genome_build: Genome build to use
        pedigree_entry: Optional PedigreeEntry for this sample

    Returns:
        List of RegionEntry namedtuples
    """
    paraphase_json_calls = unpack_json(paraphase_results.JSON)

    sample_entries = []
    for region in paraphase_json_calls:
        region_data = paraphase_json_calls[region]
        if region not in split_bams:
            continue
        bam_path = split_bams[region].BAM
        bai_path = split_bams[region].BAI

        try:
            realign_region = genomic_interval_from_str(
                paraphase_config[region]["realign_region"]
            )
        except Exception as e:
            logger.info(
                f"Failed to find config info for region f{region}, skipping. \n{e}"
            )
            continue

        total_cn, special_info = get_special_info(
            region, region_data, paraphase_results
        )

        igv_session_path = path.join(
            IGV_SESSIONS_PATH.format(sample=paraphase_results.Sample),
            f"{region}_igv.xml",
        )
        image_path = path.join(
            IMAGES_PATH.format(sample=paraphase_results.Sample),
            f"{paraphase_results.Sample}_{region}.png",
        )

        sample_entries.append(
            RegionEntry(
                realign_region.Chrom,
                max(0, realign_region.Start - REGION_PADDING),
                max(0, realign_region.End + REGION_PADDING),
                region,
                paraphase_results.Sample,
                bam_path,
                bai_path,
                total_cn,
                special_info,
                image_path,
                igv_session_path,
                pedigree_entry.FamilyID if pedigree_entry else "",
                pedigree_entry.PaternalID if pedigree_entry else "",
                pedigree_entry.MaternalID if pedigree_entry else "",
                pedigree_entry.Sex if pedigree_entry else "",
                pedigree_entry.Phenotype if pedigree_entry else "",
            )
        )
    return sample_entries
