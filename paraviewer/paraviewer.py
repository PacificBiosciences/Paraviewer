#!/usr/bin/env python

from os import path
import logging
from paraviewer.igv_builder import generate_igv_images, generate_igv_sessions
from paraviewer.process_paraphase import (
    get_paraphase_results,
    make_table_entries,
    make_trio_table_entries,
)
from paraviewer.page_builder import build_review_page
from paraviewer.process_puretarget import get_puretarget_results
from paraviewer.utils import (
    ParaphaseResults,
    PedigreeEntry,
    get_config,
    make_output_dirs,
    split_bam,
)
from . import __version__
import sys
from typing import Dict, List, Optional
from tqdm import tqdm

logger = logging.getLogger(__name__)
logging.getLogger("urllib3").setLevel(logging.CRITICAL)


def get_trio_samples(
    pedigree_dict: Dict[str, PedigreeEntry],
    all_paraphase_results: List[ParaphaseResults],
):
    """
    Get the trio samples from the pedigree dictionary and the paraphase results
    """
    trio_samples = {}
    for sample in pedigree_dict:
        if sample in all_paraphase_results:
            pedigree_entry = pedigree_dict[sample]
            if (
                pedigree_entry.PaternalID in all_paraphase_results
                and pedigree_entry.MaternalID in all_paraphase_results
            ):
                trio_samples[f"{sample}-trio"] = pedigree_entry
    return trio_samples


def read_pedigree_file(
    ped_file: Optional[str],
    include_only_samples: list[str],
    exclude_samples: list[str],
) -> Dict[str, PedigreeEntry]:
    """
    Read a GATK-format PED file and return a dictionary of PedigreeEntry objects.

    Args:
        ped_file: Path to the PED file, or None if no file provided

    Returns:
        Dictionary mapping sample IDs to PedigreeEntry objects
    """
    if not ped_file or not path.exists(ped_file):
        return {}

    pedigree_dict = {}
    try:
        with open(ped_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                fields = line.split()
                if len(fields) != 6:
                    logger.warning(f"Warning: Skipping malformed PED line: {line}")
                    continue

                family_id, individual_id, PaternalID, MaternalID, sex, phenotype = (
                    fields
                )
                if (
                    len(include_only_samples) > 0
                    and individual_id.lower() not in include_only_samples
                ):
                    continue
                if individual_id.lower() in exclude_samples:
                    continue

                # Convert empty values to empty strings
                PaternalID = PaternalID if PaternalID != "0" else ""
                MaternalID = MaternalID if MaternalID != "0" else ""

                # Recode sex values
                if sex == "1":
                    sex = "Male"
                elif sex == "2":
                    sex = "Female"
                else:
                    sex = "Unknown"

                entry = PedigreeEntry(
                    FamilyID=family_id,
                    IndividualID=individual_id,
                    PaternalID=PaternalID,
                    MaternalID=MaternalID,
                    Sex=sex,
                    Phenotype=phenotype,
                )
                pedigree_dict[individual_id] = entry

    except Exception as e:
        logger.error(f"Error reading PED file {ped_file}: {e}")
        sys.exit(1)

    return pedigree_dict


def process_individual_sample(
    sample_paraphase_results,
    pedigree_dict,
    paraphase_config,
    outdir,
    clobber,
    include_only_regions,
    exclude_regions,
    genome,
    max_reads_per_hap,
    no_igv_rerun,
):
    """
    Process an individual sample and return the region entries

    Args:
        sample_paraphase_results: ParaphaseResults namedtuple
        pedigree_dict: Dictionary mapping sample IDs to PedigreeEntry objects
        paraphase_config: Dictionary mapping region names to config information

    Returns:
        List of RegionEntry namedtuples
    """
    make_output_dirs(outdir, sample_paraphase_results.Sample, clobber)
    split_bams = split_bam(
        sample_paraphase_results.BAM,
        sample_paraphase_results.BAI,
        outdir,
        sample_paraphase_results.Sample,
        include_only_regions,
        exclude_regions,
        max_reads_per_hap,
    )
    if len(split_bams) == 0:
        logger.error(f"No specified regions found in {sample_paraphase_results.BAM}.")
        sys.exit(1)

    sample_region_entries = make_table_entries(
        sample_paraphase_results,
        pedigree_dict.get(sample_paraphase_results.Sample),
        split_bams,
        paraphase_config,
    )

    igv_batch_entries = generate_igv_sessions(sample_region_entries, outdir, genome)
    if not no_igv_rerun:
        generate_igv_images(igv_batch_entries, outdir, genome, False)
    return sample_region_entries, split_bams


def process_trio(
    trio,
    all_paraphase_results,
    all_split_bams,
    paraphase_config,
    outdir,
    clobber,
    genome,
    no_igv_rerun,
):
    """
    Process a trio sample and return the region entries
    """
    make_output_dirs(outdir, trio.IndividualID + "-trio", clobber)

    trio_region_entries = make_trio_table_entries(
        trio,
        all_paraphase_results[trio.IndividualID],
        all_paraphase_results[trio.PaternalID],
        all_paraphase_results[trio.MaternalID],
        all_split_bams,
        paraphase_config,
        outdir,
    )

    igv_batch_entries = generate_igv_sessions(trio_region_entries, outdir, genome)
    if not no_igv_rerun:
        generate_igv_images(igv_batch_entries, outdir, genome, True)
    return trio_region_entries


def paraviewer(args):
    all_region_entries = []
    all_split_bams = {}
    pedigree_dict = read_pedigree_file(
        args.pedigree, args.include_only_samples, args.exclude_samples
    )

    all_paraphase_results = {}
    if args.paraphase_dir:
        paraphase_config = get_config(args.genome, "paraphase")
        all_paraphase_results = get_paraphase_results(
            args.paraphase_dir,
            args.include_only_samples,
            args.exclude_samples,
            pedigree_dict,
        )
    elif args.ptcp_dir:
        paraphase_config = get_config(args.genome, "puretarget")
        all_paraphase_results = get_puretarget_results(
            args.ptcp_dir,
            args.include_only_samples,
            args.exclude_samples,
            pedigree_dict,
        )
    if len(all_paraphase_results) == 0:
        logger.error("No results found in input directory")
        sys.exit(1)

    trio_samples = get_trio_samples(pedigree_dict, all_paraphase_results)

    for i, sample in enumerate(
        tqdm(
            list(all_paraphase_results.keys()) + list(trio_samples.keys()),
            desc="Samples processed",
        )
    ):
        if i < len(all_paraphase_results):
            logger.info(f"Processing sample {sample}")
            sample_paraphase_results = all_paraphase_results[sample]
            sample_region_entries, sample_split_bams = process_individual_sample(
                sample_paraphase_results,
                pedigree_dict,
                paraphase_config,
                args.outdir,
                args.clobber,
                args.include_only_regions,
                args.exclude_regions,
                args.genome,
                args.max_reads_per_haplotype,
                args.no_igv_rerun,
            )
            all_region_entries += sample_region_entries
            all_split_bams[sample] = sample_split_bams
        else:
            logger.info(f"Processing trio {sample}")
            all_region_entries += process_trio(
                trio_samples[sample],
                all_paraphase_results,
                all_split_bams,
                paraphase_config,
                args.outdir,
                args.clobber,
                args.genome,
                args.no_igv_rerun,
            )
    build_review_page(args.outdir, all_region_entries)
