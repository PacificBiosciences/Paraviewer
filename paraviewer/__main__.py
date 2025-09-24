#!/usr/bin/env python
from __future__ import print_function
import json
from re import A
import shutil
import argparse
import sys
from os import path
import logging
import subprocess
import pathlib
from paraviewer.utils import is_linux, get_config

logger = logging.getLogger(__name__)

from .__init__ import __version__
from .paraviewer import paraviewer


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
    parser.add_argument(
        "--outdir",
        help="Path to output directory - should not already exist",
        required=True,
        type=valid_parent_dir,
    )
    parser.add_argument(
        "--paraphase-dir",
        help="Path to paraphase result directory.",
        required=False,
        type=valid_dir,
    )
    parser.add_argument(
        "--ptcp-dir",
        help="Path to PureTarget Carrier Panel result directory.",
        required=False,
        type=valid_dir,
    )
    parser.add_argument(
        "--clobber",
        help="Overwrite output directory if it already exists",
        action="store_true",
    )
    parser.add_argument(
        "--genome",
        help="Desired genome build. Choose between GRCh37/HG19 (hg19) and GRCh38/HG38 (hg38)",
        type=str,
        required=True,
        choices=["hg19", "hg38"],
    )
    parser.add_argument(
        "--pedigree",
        help="Path to GATK-format PED file containing pedigree information - unrepresented samples will be excluded.",
        type=valid_file,
    )
    parser.add_argument(
        "--include-only-regions",
        help="Space-delimited list of region names to include. Regions not specified will be excluded.",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--exclude-regions",
        help="Space-delimited list of region names to exclude.",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--include-only-samples",
        help="Space-delimited list of sample IDs to include. Samples not specified will be excluded.",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--exclude-samples",
        help="Space-delimited list of sample IDs to exclude.",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--max-reads-per-haplotype",
        help="Maximum number of reads to show per haplotype.",
        default=500,
    )
    parser.add_argument(
        "--verbose",
        help="Print verbose output for debugging purposes",
        action="store_true",
    )
    parser.add_argument(
        "--no-igv-rerun",
        help=argparse.SUPPRESS,
        action="store_true",
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
        valid_regions = [x.lower() for x in valid_regions]
        for i, include_item in enumerate(include_list):
            if include_item in valid_regions:
                filtered_include_list.append(include_item)
            else:
                logger.warning(
                    f"Include list `{list_name}` contains invalid entry {include_item}, which will be ignored"
                )

        for i, exclude_item in enumerate(exclude_list):
            if exclude_item in valid_regions:
                filtered_exclude_list.append(exclude_item)
            else:
                logger.warning(
                    f"Exclude list `{list_name}` contains invalid entry {exclude_item}, which will be ignored"
                )
        return filtered_include_list, filtered_exclude_list
    return include_list, exclude_list


def main():
    print("\nParaViewer v{}".format(__version__), file=sys.stderr)
    parser = setup_args()
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if not is_tool_installed_via_conda("igv"):
        error_msg = "IGV not found. Paraviewer requires that you install IGV via conda/mamba, e.g."
        error_msg += "\n`mamba install -c bioconda igv`"
        logger.error(error_msg)
        sys.exit(1)
    if is_linux() and not is_tool_installed("Xvfb"):
        error_msg = "Xvfb not found. Paraviewer on linux requires Xvfb."
        logger.error(error_msg)
        sys.exit(1)

    if not args.paraphase_dir and not args.ptcp_dir:
        logger.error("Either --paraphase-dir or --ptcp-dir must be specified")
        sys.exit(1)
    if args.paraphase_dir and args.ptcp_dir:
        logger.error(
            "--paraphase-dir and --ptcp-dir are mutually exclusive: specify only one or the other"
        )
        sys.exit(1)

    args.include_only_samples, args.exclude_samples = validate_include_exclude_lists(
        args.include_only_samples, args.exclude_samples, "sample"
    )
    source_pipeline = "paraphase"
    if args.ptcp_dir:
        source_pipeline = "puretarget"
    valid_regions = get_config(args.genome, source_pipeline).keys()
    args.include_only_regions, args.exclude_regions = validate_include_exclude_lists(
        args.include_only_regions,
        args.exclude_regions,
        "region",
        valid_regions,
    )

    paraviewer(args)


if __name__ == "__main__":
    print(
        "You are running this module directly, which should only be done for debugging",
        file=sys.stderr,
    )
    sys.exit(main() or 0)
