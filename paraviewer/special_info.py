from itertools import chain

from paraviewer.utils import ParaphaseResults


def get_copy_number(region_data: dict) -> str | dict | int:
    """
    Figure out the correct copy number from the gene_cn,
    total_cn, and highest_total_cn fields in paraphase output
    """
    total_cn = ""
    if "gene_cn" in region_data and region_data["gene_cn"] is not None:
        total_cn = region_data["gene_cn"]
    elif "total_cn" in region_data and region_data["total_cn"] is not None:
        total_cn = region_data["total_cn"]
    elif "highest_total_cn" in region_data and region_data["highest_total_cn"] is not None:
        total_cn = region_data["highest_total_cn"]
    return total_cn


def define_special_fields() -> dict[str, list[str]]:
    """
    Special info fields for known complex medically relevant gene regions.
    Gene region names are given for organization only.
    """
    cmrg_special_fields = {
        "CFH": ["fusions_called"],
        "F8": ["sv_called"],
        "GBA": ["fusions_called"],
        "HBA": ["genotype", "sv_called", "alleles_final"],
        "IKBKG": ["deletion_haplotypes"],
        "NEB": ["alleles_final"],
        "OPN1LW": ["annotated_haplotypes", "alleles_final", "annotated_alleles"],
        "RCCX": ["alleles_final", "ending_hap", "annotated_alleles"],
        "SMN1": [
            "smn1_cn",
            "smn2_cn",
            "smn2_del78_cn",
        ],
        "all_regions": [],
    }
    return cmrg_special_fields


def get_special_info(
    region_name: str, region_data: dict, paraphase_results: ParaphaseResults
) -> tuple[str | dict | int, str]:
    """
    Get total_cn field and special info field for complex medically-relevant genes
    """
    total_cn = get_copy_number(region_data)
    cmrg_special_fields = define_special_fields()
    special_info = []

    for special_field_name in sorted({str(x) for x in chain(*cmrg_special_fields.values())}):
        if special_field_name in region_data:
            special_info_field_data = region_data[special_field_name]
            if special_info_field_data is None or special_info_field_data == "NA":
                continue

            if isinstance(special_info_field_data, list):
                if len(special_info_field_data) == 0:
                    continue
                field_joined = ""
                for field in special_info_field_data:
                    if len(field_joined) > 0:
                        field_joined += ", "
                    if isinstance(field, list):
                        field_joined += " | ".join(str(x) for x in field)
                    else:
                        field_joined += field

                special_info.append(f"{special_field_name}: {field_joined}")
            elif isinstance(special_info_field_data, dict):
                if len(special_info_field_data) == 0:
                    continue
                special_info_items = [f"{key}" for key in special_info_field_data]
                special_info.append(
                    f"{special_field_name},{', '.join(special_info_items)}"
                )
            else:
                # append string data
                special_info.append(f"{special_field_name},{special_info_field_data}")

    # if this is f8, drop everything other than the inversion call
    is_f8 = "f8" in region_name.lower()
    if is_f8:
        special_info = []
    if paraphase_results.F8_INV:
        special_info.append(paraphase_results.F8_INV)
    elif not is_f8 and paraphase_results.HAVANNO and region_name in paraphase_results.HAVANNO:
        special_info.append(paraphase_results.HAVANNO[region_name])

    special_info = ";".join(special_info)
    return total_cn, special_info.strip(", ")
