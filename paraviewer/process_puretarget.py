from glob import glob
import logging
from os import path
from typing import Dict
from paraviewer.process_paraphase import get_paraphase_results
from paraviewer.utils import (
    ParaphaseResults,
    PedigreeEntry,
    parse_sample_name_from_paraphase_output,
    strip_suffix_from_path,
    unpack_json,
)

logger = logging.getLogger(__name__)


def get_puretarget_results(
    puretarget_dir: str,
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
    for paraphase_dir in glob(path.join(puretarget_dir, "*_paraphase")):
        output_file = glob(path.join(paraphase_dir, "*"))[0]
        sample = parse_sample_name_from_paraphase_output(output_file)
        sample_results = get_paraphase_results(
            paraphase_dir, include_only_samples, exclude_samples, pedigree_dict
        )
        if len(sample_results) == 0:
            continue
        base_dirname = strip_suffix_from_path(paraphase_dir, "_paraphase")
        f8_json_path = base_dirname + ".f8inversion.json"
        f8_special_info = get_f8_inv_annotation(f8_json_path)
        sample_results[sample] = sample_results[sample]._replace(F8_INV=f8_special_info)

        havanno_json_path = base_dirname + ".havanno.json"
        havanno_special_info = get_havanno_annotations(havanno_json_path)
        sample_results[sample] = sample_results[sample]._replace(
            HAVANNO=havanno_special_info
        )

        all_results.update(sample_results)

    return all_results


def get_f8_inv_annotation(f8_json_path):
    """
    read the sample-associated file from the f8inversion.json
    and add annotation for the inversion
    """
    f8_inv_info = unpack_json(f8_json_path)
    f8_inv_annotation = ""
    if f8_inv_info is None:
        return f8_inv_annotation
    if f8_inv_info["f8inv1"]["has_inversion"]:
        f8_inv_annotation += "{},".format(f8_inv_info["f8inv1"]["inversion_genotype"])
    if f8_inv_info["f8inv22"]["has_inversion"]:
        f8_inv_annotation += "{},".format(f8_inv_info["f8inv22"]["inversion_genotype"])
    return f8_inv_annotation


def get_havanno_annotations(havanno_json_path):
    """
    read the sample-associated file from the havanno.json
    and return dict of special info strings by region name
    """
    results = {}
    havanno_info = unpack_json(havanno_json_path)

    if havanno_info is None:
        return results

    for region in havanno_info["annotations"]:
        region_havanno_annotations = []
        for haplotype in havanno_info["annotations"][region]:
            if "hba_homology" in haplotype:
                continue
            annotation = havanno_info["annotations"][region][haplotype]
            pathogenic_variant_count = annotation["num_pathogenic_variants"]
            insertion_size = annotation["total_insertion_size"]
            deletion_size = annotation["total_deletion_size"]

            haplotype_annotations = []

            if pathogenic_variant_count and pathogenic_variant_count > 0:
                haplotype_annotations.append(
                    f"{pathogenic_variant_count} possible pathogenic vars"
                )

            if insertion_size and insertion_size > 0:
                haplotype_annotations.append(f"{insertion_size}bp INS")

            if deletion_size and deletion_size > 0:
                haplotype_annotations.append(f"{deletion_size}bp DEL")
            if len(haplotype_annotations) > 0:
                haplotype_annotations = ", ".join(haplotype_annotations)
                region_havanno_annotations.append(
                    f"{haplotype},{haplotype_annotations}"
                )
        if len(region_havanno_annotations) > 0:
            results[region] = ";".join(region_havanno_annotations)

    return results
