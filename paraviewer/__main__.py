#!/usr/bin/env python
from __future__ import print_function

import argparse
import difflib
import json
import logging
import pathlib
import shutil
import subprocess
import sys
from glob import glob
from os import path

from orographer.deploy import run_deploy

from paraviewer.utils import parse_sample_name_from_paraphase_output, unpack_json

from .__init__ import __version__
from .paraviewer import paraviewer


def _closest_match(value, valid_set, cutoff=0.5):
    """Return the single closest match for value in valid_set, or None."""
    if not value or not valid_set:
        return None
    matches = difflib.get_close_matches(value, valid_set, n=1, cutoff=cutoff)
    return matches[0] if matches else None


def _get_valid_sample_names(input_dir, is_puretarget):
    """
    Discover sample names from paraphase or puretarget input directory.
    Returns a list of lowercase sample names, or None if discovery fails.
    """
    if not input_dir or not path.isdir(input_dir):
        return None
    try:
        if is_puretarget:
            names = []
            for subdir in glob(path.join(input_dir, "*_paraphase")):
                files = glob(path.join(subdir, "*"))
                if files:
                    names.append(
                        parse_sample_name_from_paraphase_output(files[0]).lower()
                    )
            return list(set(names)) if names else None
        else:
            json_matches = glob(path.join(input_dir, "*paraphase.json")) + glob(
                path.join(input_dir, "*paraphase.json.gz")
            )
            if not json_matches:
                return None
            return list(
                set(
                    parse_sample_name_from_paraphase_output(json_path).lower()
                    for json_path in json_matches
                )
            )
    except OSError:
        return None


logger = logging.getLogger(__name__)


def is_tool_installed(tool_name: str) -> bool:
    if shutil.which(tool_name) is None:
        return False
    else:
        return True


def is_tool_installed_via_conda(tool_name: str) -> bool:
    is_tool_installed(tool_name)
    for conda_option in ["mamba", "conda", "micromamba"]:
        try:
            result = subprocess.run(
                [conda_option, "list", "--json"],
                capture_output=True,
                text=True,
                check=True,
            )
            packages = json.loads(result.stdout)
            return any(pkg["name"] == tool_name for pkg in packages)
        except Exception:
            continue
    return False


def valid_parent_dir(dirpath):
    parent = pathlib.Path(dirpath).parent
    if parent.exists():
        return dirpath
    logger.error(f"Parent directory {parent} does not exist")
    sys.exit(1)


def valid_dir(dirpath):
    if not path.exists(dirpath):
        logger.error(f"Directory {dirpath} does not exist")
        sys.exit(1)
    if not path.isdir(dirpath):
        logger.error(f"{dirpath} is not a directory")
        sys.exit(1)
    return dirpath


def valid_file(filepath):
    if not path.exists(filepath):
        logger.error(f"File {filepath} does not exist")
        sys.exit(1)
    return filepath


def setup_args():
    parser = argparse.ArgumentParser(
        prog="paraviewer", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-v",
        "--version",
        help="Installed version ({})".format(__version__),
        action="version",
        version="%(prog)s " + str(__version__),
    )

    # Create subparsers for commands
    subparsers = parser.add_subparsers(
        dest="command", help="Command to run", metavar="COMMAND", required=True
    )

    # Create command
    create_parser = subparsers.add_parser(
        "create",
        help="Generate paraviewer HTML page with orographer plots",
        description=(
            "Process paraphase or puretarget results and generate "
            "interactive HTML viewer with orographer plots."
        ),
    )

    # Argument groups for nicer --help organization
    required_args_group = create_parser.add_argument_group("Required")
    filtering_args_group = create_parser.add_argument_group("Filtering")
    annotation_args_group = create_parser.add_argument_group("Annotation")
    other_args_group = create_parser.add_argument_group("Other")

    # Required (inputs/outputs and reference)
    required_args_group.add_argument(
        "--outdir",
        help="Path to output directory - should not already exist",
        required=True,
        type=valid_parent_dir,
    )
    required_args_group.add_argument(
        "--paraphase-dir",
        help="EITHER path to paraphase result directory.",
        required=False,
        type=valid_dir,
    )
    required_args_group.add_argument(
        "--ptcp-dir",
        help="OR path to PureTarget Carrier Panel result directory.",
        required=False,
        type=valid_dir,
    )

    required_args_group.add_argument(
        "--ref",
        help="Path to reference FASTA file",
        type=valid_file,
        required=True,
    )
    # Annotation (pedigree and gene models)
    annotation_args_group.add_argument(
        "--gtf",
        help="Optional path to bgzip+tabix GTF/GFF3 for gene track.",
        type=valid_file,
        required=False,
    )

    # Filtering (sample/region lists)
    filtering_args_group.add_argument(
        "--include-only-regions",
        help="Region names to include; others excluded.",
        type=str,
        nargs="+",
        required=False,
    )
    filtering_args_group.add_argument(
        "--exclude-regions",
        help="Space-delimited list of region names to exclude.",
        type=str,
        nargs="+",
        required=False,
    )
    annotation_args_group.add_argument(
        "--pedigree",
        help="Optional GATK-format PED; unrepresented samples excluded.",
        type=valid_file,
        required=False,
    )
    filtering_args_group.add_argument(
        "--include-only-samples",
        help="Sample IDs to include; others excluded.",
        type=str,
        nargs="+",
        required=False,
    )
    filtering_args_group.add_argument(
        "--exclude-samples",
        help="Space-delimited list of sample IDs to exclude.",
        type=str,
        nargs="+",
        required=False,
    )

    # Oother
    other_args_group.add_argument(
        "--max-reads-per-haplotype",
        help="Maximum number of reads to show per haplotype.",
        default=500,
        required=False,
    )
    other_args_group.add_argument(
        "--threads",
        help="Number of worker processes, up to CPU count (default 1)",
        type=int,
        default=1,
        required=False,
    )
    other_args_group.add_argument(
        "--clobber",
        help="Overwrite output directory if it already exists",
        action="store_true",
        required=False,
    )
    other_args_group.add_argument(
        "--verbose",
        help="Print verbose output for debugging purposes",
        action="store_true",
    )
    other_args_group.add_argument(
        "--task-timeout",
        type=int,
        default=600,
        help=argparse.SUPPRESS,
        metavar="SECONDS",
    )

    # Deploy command
    deploy_parser = subparsers.add_parser(
        "deploy",
        help="Start HTTP server to serve paraviewer output",
        description="Start HTTP server to serve HTML.",
    )

    deploy_parser.add_argument(
        "--outdir",
        help="Directory path containing HTML and JSON files to serve",
        required=True,
        type=valid_dir,
    )

    deploy_parser.add_argument(
        "--port",
        help="Port number to serve on (default: 8000)",
        type=int,
        default=8000,
    )

    return parser


def validate_include_exclude_lists(
    include_list, exclude_list, list_name, valid_regions=None
):
    if include_list is None:
        include_list = []
    if exclude_list is None:
        exclude_list = []
    for item in include_list:
        if exclude_list and item in exclude_list:
            logger.error(f"{item} is in both include and exclude {list_name} lists")
            sys.exit(1)
    include_list = [item.lower() for item in include_list]
    exclude_list = [item.lower() for item in exclude_list]

    if valid_regions:
        filtered_include_list = []
        filtered_exclude_list = []
        valid_lower = [x.lower() for x in valid_regions]
        for _i, include_item in enumerate(include_list):
            if include_item in valid_lower:
                filtered_include_list.append(include_item)
            else:
                suggestion = _closest_match(include_item, valid_lower)
                if suggestion:
                    logger.warning(
                        "Include list `%s` contains invalid entry %s, ignored "
                        "(did you mean: %s?)",
                        list_name,
                        include_item,
                        suggestion,
                    )
                else:
                    logger.warning(
                        "Include list `%s` contains invalid entry %s, ignored",
                        list_name,
                        include_item,
                    )

        for _i, exclude_item in enumerate(exclude_list):
            if exclude_item in valid_lower:
                filtered_exclude_list.append(exclude_item)
            else:
                suggestion = _closest_match(exclude_item, valid_lower)
                if suggestion:
                    logger.warning(
                        "Exclude list `%s` contains invalid entry %s, ignored "
                        "(did you mean: %s?)",
                        list_name,
                        exclude_item,
                        suggestion,
                    )
                else:
                    logger.warning(
                        "Exclude list `%s` contains invalid entry %s, ignored",
                        list_name,
                        exclude_item,
                    )
        return filtered_include_list, filtered_exclude_list
    return include_list, exclude_list


def _collect_union_region_keys(input_dir: str, is_puretarget: bool) -> set:
    """Lowercase region names present in any *paraphase.json(.gz) under input_dir."""
    keys: set = set()
    if is_puretarget:
        for paraphase_dir in glob(path.join(input_dir, "*_paraphase")):
            jm = glob(path.join(paraphase_dir, "*paraphase.json")) + glob(
                path.join(paraphase_dir, "*paraphase.json.gz")
            )
            for jp in jm:
                data = unpack_json(jp)
                if isinstance(data, dict):
                    keys.update(
                        str(k).lower() for k in data if isinstance(k, str)
                    )
    else:
        jm = glob(path.join(input_dir, "*paraphase.json")) + glob(
            path.join(input_dir, "*paraphase.json.gz")
        )
        for jp in jm:
            data = unpack_json(jp)
            if isinstance(data, dict):
                keys.update(str(k).lower() for k in data if isinstance(k, str))
    return keys


def validate_region_filters_strict(
    include_list, exclude_list, union_region_keys: set
) -> tuple:
    """
    Require every include/exclude region name to appear in union_region_keys.
    Exits on unknown name or overlap.
    """
    if include_list is None:
        include_list = []
    if exclude_list is None:
        exclude_list = []
    include_list = [str(x).lower() for x in include_list]
    exclude_list = [str(x).lower() for x in exclude_list]
    for item in include_list:
        if item in exclude_list:
            logger.error("%r is in both include and exclude region lists", item)
            sys.exit(1)
    if not union_region_keys and (include_list or exclude_list):
        logger.error(
            "No region keys found in input Paraphase JSON; cannot validate "
            "--include-only-regions / --exclude-regions."
        )
        sys.exit(1)
    for item in include_list:
        if item not in union_region_keys:
            sug = _closest_match(item, union_region_keys)
            if sug:
                logger.error(
                    "Unknown region %r (not in input JSON). Did you mean %r?",
                    item,
                    sug,
                )
            else:
                logger.error(
                    "Unknown region %r: not present in any sample's Paraphase JSON.",
                    item,
                )
            sys.exit(1)
    for item in exclude_list:
        if item not in union_region_keys:
            sug = _closest_match(item, union_region_keys)
            if sug:
                logger.error(
                    "Unknown region %r in --exclude-regions (not in input JSON). "
                    "Did you mean %r?",
                    item,
                    sug,
                )
            else:
                logger.error(
                    "Unknown region %r in --exclude-regions: not in input JSON.",
                    item,
                )
            sys.exit(1)
    return include_list, exclude_list


def main():
    print("\nParaViewer v{}".format(__version__), file=sys.stderr)
    parser = setup_args()
    args = parser.parse_args()

    if args.command == "deploy":
        run_deploy(args.outdir, args.port)
        return

    # Create command
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if not args.paraphase_dir and not args.ptcp_dir:
        logger.error("Either --paraphase-dir or --ptcp-dir must be specified")
        sys.exit(1)
    if args.paraphase_dir and args.ptcp_dir:
        logger.error(
            "--paraphase-dir and --ptcp-dir are mutually exclusive; use one only."
        )
        sys.exit(1)

    source_pipeline = "paraphase"
    input_dir = args.paraphase_dir
    if args.ptcp_dir:
        source_pipeline = "puretarget"
        input_dir = args.ptcp_dir
    valid_samples = None
    if args.include_only_samples or args.exclude_samples:
        valid_samples = _get_valid_sample_names(
            input_dir, source_pipeline == "puretarget"
        )
    had_include_only_samples = bool(args.include_only_samples)
    args.include_only_samples, args.exclude_samples = validate_include_exclude_lists(
        args.include_only_samples, args.exclude_samples, "sample", valid_samples
    )
    if had_include_only_samples and len(args.include_only_samples) == 0:
        logger.error(
            "--include-only-samples was set but none of the given sample names "
            "are valid; exiting."
        )
        sys.exit(1)

    union_keys = _collect_union_region_keys(
        input_dir, source_pipeline == "puretarget"
    )
    args.include_only_regions, args.exclude_regions = validate_region_filters_strict(
        args.include_only_regions,
        args.exclude_regions,
        union_keys,
    )

    exit_code = paraviewer(
        args.paraphase_dir,
        args.ptcp_dir,
        args.outdir,
        args.pedigree,
        args.include_only_samples,
        args.exclude_samples,
        args.include_only_regions,
        args.exclude_regions,
        args.max_reads_per_haplotype,
        args.ref,
        args.gtf,
        args.clobber,
        args.task_timeout,
        args.threads,
    )
    sys.exit(exit_code or 0)


if __name__ == "__main__":
    print(
        "You are running this module directly, which should only be done for debugging",
        file=sys.stderr,
    )
    sys.exit(main() or 0)
