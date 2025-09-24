#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create IGV images and session files for display in web viewer
"""
from __future__ import print_function
import subprocess
import tempfile

from collections import namedtuple
import logging
import sys
import os
import time
import atexit

from paraviewer.utils import BAMS_PATH, IGV_SESSIONS_PATH, RegionEntry, is_mac, is_linux

logger = logging.getLogger(__name__)

TIMEOUT_SECONDS = 1200
BATCH_SIZE = 100

"""
IGV automatically listens for localhost connections to port 60151.
These allow it to accept socket connections and use them to define commands.
"""
IGV_HOST = "127.0.0.1"
IGV_PORT = 60151
IGV_MAX_WAIT = 60
IGV_WAIT = 30
IGV_VIRTUAL_SCREEN = "1920x1080x24"
IGV_VIRTUAL_DISPLAY_NUMBER = 99
IGV_GENOME_URL = "genome https://raw.githubusercontent.com/igvteam/igv-data/refs/heads/main/genomes/json/{}.json"
VALID_GENOMES = ["hg38", "hg19"]

IGVSessionFields = namedtuple(
    "IGVSessionFields", ["Genome", "Chrom", "Start", "End", "BAM", "BAI", "Sample"]
)
TrioIGVSessionFields = namedtuple(
    "IGVSessionFields",
    [
        "Genome",
        "Chrom",
        "Start",
        "End",
        "BAM",
        "BAI",
        "PaternalBAM",
        "PaternalBAI",
        "MaternalBAM",
        "MaternalBAI",
        "Sample",
        "PaternalID",
        "MaternalID",
    ],
)

"""
Trio IGV session template. Must be used by calling format() with a TrioIGVSessionFields namedtuple
"""
TRIO_IGV_SESSION_TEMPLATE = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{0.Genome}" locus="{0.Chrom}:{0.Start}-{0.End}" version="8">
    <Resources>
        <Resource path="{0.PaternalBAM}" type="bam"/>
        <Resource path="{0.MaternalBAM}" type="bam"/>
        <Resource path="{0.BAM}" type="bam"/>
    </Resources>
    <Panel height="427" name="Panel1745015156681" width="2543">
        <Track attributeKey="{0.PaternalBAM} Coverage" autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" fontSize="12" id="{0.PaternalBAM}_coverage" name="Paternal Coverage" snpThreshold="0.2" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="87.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track attributeKey="{0.PaternalBAM}" clazz="org.broad.igv.sam.AlignmentTrack" color="185,185,185" displayMode="SQUISHED" experimentType="THIRD_GEN" fontSize="12" id="{0.PaternalBAM}" name="Paternal BAM" visible="true">
            <RenderOptions colorOption="YC_TAG" groupByOption="TAG" groupByTag="HP"/>
        </Track>
    </Panel>
    <Panel height="382" name="Panel1748548125549" width="2543">
        <Track attributeKey="{0.MaternalBAM} Coverage" autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" fontSize="12" id="{0.MaternalBAM}_coverage" name="Maternal Coverage" snpThreshold="0.2" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="131.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track attributeKey="{0.MaternalBAM}" clazz="org.broad.igv.sam.AlignmentTrack" color="185,185,185" displayMode="SQUISHED" experimentType="THIRD_GEN" fontSize="12" id="{0.MaternalBAM}" name="Maternal BAM" visible="true">
            <RenderOptions colorOption="YC_TAG" groupByOption="TAG" groupByTag="HP"/>
        </Track>
    </Panel>
    <Panel height="394" name="Panel1748548125582" width="2543">
        <Track attributeKey="{0.BAM} Coverage" autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" fontSize="12" id="{0.BAM}_coverage" name="{0.Sample} Coverage" snpThreshold="0.2" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="144.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track attributeKey="{0.BAM}" clazz="org.broad.igv.sam.AlignmentTrack" color="185,185,185" displayMode="SQUISHED" experimentType="THIRD_GEN" fontSize="12" id="{0.BAM}" name="{0.Sample}" visible="true">
            <RenderOptions colorOption="YC_TAG" groupByOption="TAG" groupByTag="HP"/>
        </Track>
    </Panel>
    <Panel height="40" name="FeaturePanel" width="2543">
        <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="12" id="Reference sequence" name="Reference sequence" sequenceTranslationStrandValue="+" shouldShowTranslation="false" visible="true"/>
        <Track attributeKey="Refseq Genes" clazz="org.broad.igv.track.FeatureTrack" colorScale="ContinuousColorScale;0.0;836.0;255,255,255;0,0,178" fontSize="10" groupByStrand="false" id="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz" name="Refseq Genes" visible="true"/>
    </Panel>
    <PanelLayout dividerFractions="0.005012531328320802,0.34502923976608185,0.647451963241437,0.9649122807017544,0.9724310776942355"/>
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>
"""


"""
IGV session template. Must be used by calling format() with an IGVSessionFields namedtuple
"""
IGV_SESSION_TEMPLATE = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{0.Genome}" locus="{0.Chrom}:{0.Start}-{0.End}" version="8">
    <Resources>
        <Resource path="{0.BAM}" type="bam"/>
    </Resources>

    <Panel height="3566" name="Panel1745015156681" width="2543">
        <Track attributeKey="{0.Sample}.paraphase.bam Coverage" autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" fontSize="12" id="{0.BAM}_coverage" name="{0.BAM} Coverage" snpThreshold="0.2" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="54.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track attributeKey="{0.BAM}" clazz="org.broad.igv.sam.AlignmentTrack" color="185,185,185" displayMode="SQUISHED" experimentType="THIRD_GEN" fontSize="12" id="{0.BAM}" name="{0.BAM}" visible="true">
            <RenderOptions colorOption="YC_TAG" groupByOption="TAG" groupByTag="HP"/>
        </Track>
    </Panel>
    
    <Panel height="163" name="FeaturePanel" width="2543">
        <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="12" id="Reference sequence" name="Reference sequence" sequenceTranslationStrandValue="+" shouldShowTranslation="false" visible="true"/>
        <Track attributeKey="Refseq Genes" clazz="org.broad.igv.track.FeatureTrack" colorScale="ContinuousColorScale;0.0;836.0;255,255,255;0,0,178" fontSize="10" groupByStrand="false" id="https://hgdownload.soe.ucsc.edu/goldenPath/{0.Genome}/database/ncbiRefSeq.txt.gz" name="Refseq Genes" visible="true"/>
    </Panel>
    
    <PanelLayout dividerFractions="0.005012531328320802,0.7794486215538847,0.9573934837092731"/>
    
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>"""


"""
IGV batch script entry template. Must be used by calling format() with:
* session path (str)
* output directory (str)
* sample (str)
* region (str)
"""
IGV_BATCH_TEMPLATE = """
new
load {}
snapshotDirectory {}
snapshot {}_{}.png
"""


def validate_IGV_session_fields(session_fields: IGVSessionFields):
    """
    Make sure that the IGV session fields are valid.
    Exit if not.
    """
    if session_fields.Genome not in VALID_GENOMES:
        logger.error("Invalid genome build {}".format(session_fields.Genome))
        sys.exit(1)
    valid_chromosomes = ["{}".format(x) for x in range(1, 23)] + ["X", "Y"]
    if session_fields.Chrom.strip("chr") not in valid_chromosomes:
        logger.error("Invalid chromosome {}".format(session_fields.Chrom))
        sys.exit(1)
    if type(session_fields.Start) is not int or type(session_fields.End) is not int:
        logger.error(
            "Non-integer coordinate in {}-{}".format(
                session_fields.Start, session_fields.End
            )
        )
        sys.exit(1)


def generate_igv_sessions(
    sample_region_entries: list[RegionEntry],
    outdir: str,
    genome: str,
):
    """
    * Create an IGV session file for each region in sample_region_entries.
    * Write the IGV session to outdir/data/sample/.
    * Write an IGV batch script section to generate the IGV snapshot from the session file.
    * Return list of IGV batch strings
    """
    igv_batch_entries = []
    for region_entry in sample_region_entries:
        if type(region_entry.BAM) != list:
            session_info = IGVSessionFields(
                genome,
                region_entry.Chrom,
                region_entry.Start,
                region_entry.End,
                os.path.basename(region_entry.BAM),
                os.path.basename(region_entry.BAI),
                region_entry.Sample,
            )
            validate_IGV_session_fields(session_info)
            igv_session = IGV_SESSION_TEMPLATE.format(session_info)
        elif len(region_entry.BAM) == 3:
            session_info = TrioIGVSessionFields(
                Genome=genome,
                Chrom=region_entry.Chrom,
                Start=region_entry.Start,
                End=region_entry.End,
                PaternalBAM=os.path.basename(region_entry.BAM[0]),
                PaternalBAI=os.path.basename(region_entry.BAI[0]),
                MaternalBAM=os.path.basename(region_entry.BAM[1]),
                MaternalBAI=os.path.basename(region_entry.BAI[1]),
                BAM=os.path.basename(region_entry.BAM[2]),
                BAI=os.path.basename(region_entry.BAI[2]),
                Sample=region_entry.Sample,
                PaternalID=region_entry.PaternalID,
                MaternalID=region_entry.MaternalID,
            )
            validate_IGV_session_fields(session_info)
            igv_session = TRIO_IGV_SESSION_TEMPLATE.format(session_info)

        # write one session file to the permanent location for later download
        igv_session_name = os.path.join(
            outdir,
            IGV_SESSIONS_PATH.format(sample=region_entry.Sample),
            f"{region_entry.Region}_igv.xml",
        )
        with open(igv_session_name, "wt") as igv_out:
            print(igv_session, file=igv_out)

        # write one session file to a temp location for immediate use
        tmp_igv_session_name = os.path.join(
            outdir,
            BAMS_PATH.format(sample=region_entry.Sample),
            f".{region_entry.Region}_igv.xml",
        )
        with open(tmp_igv_session_name, "wt") as igv_out:
            print(igv_session, file=igv_out)

        snapshot_directory = os.path.join(outdir, "data", region_entry.Sample, "images")
        igv_batch_entry = IGV_BATCH_TEMPLATE.format(
            tmp_igv_session_name,
            snapshot_directory,
            region_entry.Sample,
            region_entry.Region,
        )
        igv_batch_entries.append(igv_batch_entry)
    return igv_batch_entries


def write_prefs_file(outdir: str, is_trio: bool):
    """
    Writes a temporary file for image size preferences in IGV
    """
    left = 780
    top = 128
    width = 1546
    height = 300
    if is_trio:
        height *= 3

    prefs_file_path = os.path.join(outdir, "prefs.txt")
    with open(prefs_file_path, "w") as prefs_fh:
        print(f"IGV.Bounds={left},{top},{width},{height}", file=prefs_fh)
    return prefs_file_path


def run_batch_mac(batch_cmd, outdir):
    """
    Deploy the IGV GUI to generate images in batch
    mode using the given command.
    """
    # Run IGV in batch mode
    with open(os.path.join(outdir, "igv.log"), "a") as igv_log:
        complete_process = subprocess.run(
            batch_cmd,
            check=True,
            stderr=igv_log,
            stdout=igv_log,
            timeout=TIMEOUT_SECONDS,
        )
        complete_process.check_returncode()


def find_free_display(start, end):
    for display_num in range(start, end + 10):
        lock_file = f"/tmp/.X{display_num}-lock"
        if not os.path.exists(lock_file):
            return display_num
    logger.error(f"No free X displays found (range: {start}-{end+10}).")
    sys.exit(1)


def run_batch_linux(batch_cmd, outdir, batch_count):
    """
    Deploy IGV in headless mode with Xvfb to generate
    images using the given commands.
    """
    # Start Xvfb
    virtual_display_num = find_free_display(
        IGV_VIRTUAL_DISPLAY_NUMBER, IGV_VIRTUAL_DISPLAY_NUMBER + batch_count
    )
    xvfb = subprocess.Popen(
        ["Xvfb", f":{virtual_display_num}", "-screen", "0", IGV_VIRTUAL_SCREEN]
    )
    
    try:
        time.sleep(1)  # Give Xvfb time to initialize

        # Environment with virtual display
        env = os.environ.copy()
        env["DISPLAY"] = f":{virtual_display_num}"

        with open(os.path.join(outdir, "igv.log"), "a") as igv_log:
            complete_process = subprocess.run(
                batch_cmd,
                env=env,
                stdout=igv_log,
                stderr=igv_log,
                timeout=TIMEOUT_SECONDS,
            )
            complete_process.check_returncode()
    finally:
        # Ensure Xvfb is always terminated, regardless of success or failure
        xvfb.terminate()
        xvfb.wait()  # Wait for the process to actually terminate


def write_batch_scripts(outdir, genome, igv_batch_entries):
    """
    Write the batch scripts, given command entries.
    Each may have a maximum of BATCH_SIZE entries.
    Return the filenames.
    """
    batch_files = []
    for i, entry in enumerate(igv_batch_entries):
        # only put a bax of BATCH_SIZE entries per file
        file_idx = int(i / BATCH_SIZE)
        if len(batch_files) < file_idx + 1:
            batch_files.append(
                tempfile.NamedTemporaryFile(
                    mode="w", delete=False, suffix=".txt", dir=outdir
                )
            )
            batch_files[-1].write(IGV_GENOME_URL.format(genome))

        batch_file = batch_files[file_idx]

        batch_file.write(entry + "\n")

    batch_filenames = []
    for batch_file in batch_files:
        batch_filenames.append(batch_file.name)
        batch_file.write("exit\n")
        batch_file.close()

    return batch_filenames


def generate_igv_images(igv_batch_entries, outdir, genome, is_trio):
    """
    With a list of IGV batch entries, write them to file
    and then call IGV to generate the snapshots
    """
    if not is_mac() and not is_linux():
        logger.error(
            "Paraviewer supports only MacOs and Linus. Platform not recognized, exiting"
        )
        sys.exit(1)

    if len(igv_batch_entries) == 0:
        logger.error("No valid regions for IGV batch image generation")
        sys.exit(1)

    batch_filenames = write_batch_scripts(outdir, genome, igv_batch_entries)
    preferences_file = write_prefs_file(outdir, is_trio)

    for batch_filename in batch_filenames:
        cmd = [
            "igv",
            "-b",
            batch_filename,
            "--preferences",
            preferences_file,
        ]

        try:
            if is_mac():
                run_batch_mac(cmd, outdir)
            elif is_linux():
                run_batch_linux(cmd, outdir, len(batch_filename))
        except Exception as err:
            logger.error("IGV batch failed:\n{}\n{}".format(batch_filename, err))
        finally:
            os.remove(batch_filename)
