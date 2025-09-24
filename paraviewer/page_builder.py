#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create HTML viewer for Paraphase genomic variant results.
"""
from __future__ import print_function

import logging
import os
import sys
from .__init__ import __version__

from jinja2 import Environment, FileSystemLoader, select_autoescape
import shutil
from .utils import RegionEntry

logger = logging.getLogger(__name__)


def write_site(table_data, out_dir):
    """
    Generate the HTML page for the viewer

    Args:
        table_data: Dictionary of sample data
        out_dir: Output directory path
    """
    # grab the template
    template_dir = os.path.join(os.path.dirname(__file__), "templates")
    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=select_autoescape(["html"]),
    )

    try:
        html_template = env.get_template("paraviewer.html")
    except Exception as e:
        logger.error(f"Failed to load template: {e}")
        sys.exit(1)

    # Copy static assets to output directory
    for filename in [
        "logo_Paraphase.svg",
        "github-mark-white.png",
        "js",
        "stylesheets",
    ]:
        src = os.path.join(template_dir, "static", filename)
        destination = os.path.join(out_dir, filename)

        if os.path.isfile(src):
            shutil.copy2(src, destination)
        elif os.path.isdir(src):
            shutil.copytree(src, destination, dirs_exist_ok=True)
        else:
            logger.error("Failed to copy static template files to outdir")
            sys.exit(1)

    # write index.html
    with open("{out_dir}/index.html".format(out_dir=out_dir), "w") as fh:
        print(
            html_template.render(
                data=table_data,
            ),
            file=fh,
        )


def generate_table(sample_entries: list[RegionEntry]):
    """
    Convert sample_entries list to a list of data entries for the HTML table

    Args:
        sample_entries: Dictionary of RegionEntry objects

    Returns:
        Dictionary containing table data and column visibility flags
    """
    table_data = []
    has_pedigree_columns = {
        "FamilyID": False,
        "PaternalID": False,
        "MaternalID": False,
        "Sex": False,
        "Phenotype": False,
    }

    for entry in sample_entries:
        # Check each pedigree field independently
        if entry.FamilyID:
            has_pedigree_columns["FamilyID"] = True
        if entry.PaternalID:
            has_pedigree_columns["PaternalID"] = True
        if entry.MaternalID:
            has_pedigree_columns["MaternalID"] = True
        if entry.Sex:
            has_pedigree_columns["Sex"] = True
        if entry.Phenotype:
            has_pedigree_columns["Phenotype"] = True

        table_data.append(entry._asdict())

    # Sort by chromosome, start position, and end position
    table_data.sort(key=lambda x: (x["Chrom"], x["Start"], x["End"]))

    return {"table_data": table_data, "has_pedigree_columns": has_pedigree_columns}


def build_review_page(outdir: str, sample_entries: list[RegionEntry]):
    """
    Generate a review html page for the paraphase results from one or more samples

    Args:
        outdir: Output directory path
        sample_entries: Dictionary of RegionEntry objects
    """
    if len(sample_entries) == 0:
        logger.error("No sample/region entries found, exiting")
        sys.exit(1)

    # Convert sample_entries to the format expected by the template
    table_data = generate_table(sample_entries)

    # Generate the HTML page
    write_site(table_data, outdir)

    resultpath = os.path.join(outdir, "index.html")
    logger.info(f"Open {resultpath} in your browser to view Paraphase results")
