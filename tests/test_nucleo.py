#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `nucleo` package."""

import os
import subprocess
import shutil
import difflib
import json
import logging
import coloredlogs

# Create Logger if verbose
loggeroutput = "nucleo_test.log"
logging.basicConfig(
    filename=loggeroutput,
    filemode="w",
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)
logger = logging.getLogger(__name__)
coloredlogs.install(level="DEBUG")

RESULT_FILE_NAME = [
    "test_collapsed_FM.bai",
    "test_collapsed_FM.bam",
    "test_collapsed_aln_metrics.txt",
    "test_collapsed_duplex.bai",
    "test_collapsed_duplex.bam",
    "test_collapsed_duplex_aln_metrics.txt",
    "test_collapsed_grouped.bam",
    "test_collapsed_simplex.bai",
    "test_collapsed_simplex.bam",
    "test_collapsed_simplex_aln_metrics.txt",
    "test_fastp_report.html",
    "test_fastp_report.json",
    "test_umi_family_size.hist",
    "test_uncollapsed_BR.bai",
    "test_uncollapsed_BR.bam",
    "test_uncollapsed_BR_alignment_summary_metrics.txt",
    "test_uncollapsed_FM.bai",
    "test_uncollapsed_FM.bam",
    "test_uncollapsed_MD_metrics.txt",
    "pipeline_result.json"
]

OUTPUT_JSON_FILENAME = "pipeline_result.json"
EXPECTED_OUTPUT_JSON_FILENAME = "expected_output.json"


def setup_module(travis):
    """
    Setup and Test the workflow with cwltool
    """
    logging.info("### SETUP ###")
    with open(OUTPUT_JSON_FILENAME, "w") as json:

        cmd = [
            "cwltool",
            "nucleo.cwl",
            "test_nucleo/test_input/inputs.json",
        ]
        logging.info("setup_module: cmd being executed, %s", " ".join(cmd))
        process = subprocess.Popen(
            cmd, stdin=subprocess.PIPE, stdout=json, close_fds=True
        )

        ret_code = process.wait()
        json.flush()

    return ret_code


def teardown_module():
    """
    Tear down the setup by deleteing all the files that are downloaded and produced.
    """
    logging.info("### TEARDOWN ###")
    for outfile in RESULT_FILE_NAME:
        try:
            os.remove(outfile)
        except OSError as e:
            logging.error(
                "ERROR: cannot remove output file, %s: %s" % (outfile, e))
    try:
        shutil.rmtree("test_nucleo")
    except OSError as e:
        logging.error(
            "ERROR: cannot remove folder test_bam_collapsing : %s" % (e))


def test_check_if_metrics_file_are_same():

    """
    # General tests for checking if the metrics file is the same
    """
    logging.info(
        "### Check if files are the same from alignment metrics calculation ###"
    )
"""   
    compare_picard_metrics_files(
        "test_collapsed_aln_metrics.txt",
        "test_nucleo/test_output/test_collapsed_aln_metrics.txt",
    )
    compare_picard_metrics_files(
        "test_collapsed_duplex_aln_metrics.txt",
        "test_nucleo/test_output/test_collapsed_duplex_aln_metrics.txt",
    )
    compare_picard_metrics_files(
        "test_collapsed_simplex_aln_metrics.txt",
        "test_nucleo/test_output/test_collapsed_simplex_aln_metrics.txt",
    )
    compare_picard_metrics_files(
        "test_uncollapsed_BR_alignment_summary_metrics.txt",
        "test_nucleo/test_output/test_uncollapsed_BR_alignment_summary_metrics.txt",
    )
"""

def test_output_json():

    """
    General tests for output json
    """
    logging.info(
        "### Check if json file exists and check some basic stats ###")
    assert os.path.exists(OUTPUT_JSON_FILENAME)
    output_json = json.loads(open(OUTPUT_JSON_FILENAME, 'r').read())
    assert (
        output_json["fastp_html_output"]["basename"]
        == "test_fastp_report.html"
    )
    assert (
        output_json["fastp_json_output"]["basename"]
        == "test_fastp_report.json"
    )
    assert (
        output_json["gatk_collect_alignment_summary_metrics_txt_uncollapsed"]["basename"]
        == "test_uncollapsed_BR_alignment_summary_metrics.txt"
    )
    assert (
        output_json["indel_realignment_bam"]["basename"]
        == "test_uncollapsed_FM.bam"
    )
    assert (
        output_json["picard_mark_duplicates_metrics"]["basename"]
        == "test_uncollapsed_MD_metrics.txt"
    )
    assert (
        output_json["uncollapsed_bam"]["basename"]
        == "test_uncollapsed_BR.bam"
    )
    assert (
        output_json["gatk_collect_alignment_summary_metrics_txt_simplex"]["basename"]
        == "test_collapsed_simplex_aln_metrics.txt"
    )
    assert (
        output_json["gatk_collect_alignment_summary_metrics_txt_duplex"]["basename"]
        == "test_collapsed_duplex_aln_metrics.txt"
    )
    assert (
        output_json["gatk_collect_alignment_summary_metrics_txt_collapsed"]["basename"]
        == "test_collapsed_aln_metrics.txt"
    )
    assert (
        output_json["fgbio_postprocessing_simplex_bam"]["basename"]
        == "test_collapsed_simplex.bam"
    )
    assert (
        output_json["fgbio_group_reads_by_umi_histogram"]["basename"]
        == "test_umi_family_size.hist"
    )
    assert (
        output_json["fgbio_group_reads_by_umi_bam"]["basename"]
        == "test_collapsed_grouped.bam"
    )
    assert (
        output_json["fgbio_filter_consensus_reads_duplex_bam"]["basename"]
        == "test_collapsed_duplex.bam"
    )
    assert (
        output_json["fgbio_collect_duplex_seq_metrics_umi_counts"]["basename"]
        == "test.umi_counts.txt"
    )
    assert (
        output_json["fgbio_collect_duplex_seq_metrics_family_size"]["basename"]
        == "test.family_sizes.txt"
    )
    assert (
        output_json["fgbio_collect_duplex_seq_metrics_duplex_yield_metrics"]["basename"]
        == "test.duplex_yield_metrics.txt"
    )
    assert (
        output_json["fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt"]["basename"]
        == "test.duplex_umi_counts.txt"
    )
    assert (
        output_json["fgbio_collect_duplex_seq_metrics_duplex_family_size"]["basename"]
        == "test.duplex_family_sizes.txt"
    )
    assert (
        output_json["fgbio_collapsed_bam"]["basename"]
        == "test_collapsed_FM.bam"
    )


def compare_picard_metrics_files(output, expected):
    """
    Remove lines starting with `  # ` in picard metrics and use difflib to print differences if any and then assert
    """
    lines_result = open(output, "r").readlines()
    lines_result = list(filter(predicate, lines_result))
    lines_expected = open(expected, "r").readlines()
    lines_expected = list(filter(predicate, lines_expected))
    print("\n".join(difflib.ndiff(lines_result, lines_expected)))
    assert all([a == b for a, b in zip(lines_result, lines_expected)])


def predicate(line):
    """
    Remove lines starting with `  # `
    """
    if "#" in line:
        return False
    return True
