#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `nucleo` package."""

import os
import subprocess
import shutil
import difflib
import json

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
    print("\n### SETUP ###\n")
    with open(OUTPUT_JSON_FILENAME, "w") as json:

        cmd = [
            "cwltool",
            "--preserve-environment",
            "PATH",
            "nucleo.cwl",
            "test_nucleo/test_input/inputs.json",
        ]

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
    print("\n### TEARDOWN ###\n")
    for outfile in RESULT_FILE_NAME:
        try:
            os.remove(outfile)
        except OSError as e:
            print("ERROR: cannot remove output file, %s: %s" % (outfile, e))
    try:
        shutil.rmtree("test_nucleo")
    except OSError as e:
        print("ERROR: cannot remove folder test_nucleo : %s" % (e))


""" def test_check_if_metrics_file_are_same():
    """
    General tests for checking if the metrics file is the same
    """
    print("\n### Check if files are the same from alignment metrics calculation ###\n")

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

    # Todo: info.txt, md metrics, trimming report
 """

def test_output_json():
    """
    General tests for output json
    """
    assert os.path.exists(OUTPUT_JSON_FILENAME)
    OUTPUT_JSON = json.loads(open(OUTPUT_JSON_FILENAME, 'r').read())
    # Todo: use constant instead of magic number
    assert len(OUTPUT_JSON) == 20


def compare_picard_metrics_files(output, expected):
    """
    Remove lines starting with `#` in picard metrics and use difflib to print differences if any and then assert
    """
    lines_result = open(output, "r").readlines()
    lines_result = list(filter(predicate, lines_result))
    lines_expected = open(expected, "r").readlines()
    lines_expected = list(filter(predicate, lines_expected))
    print("\n".join(difflib.ndiff(lines_result, lines_expected)))
    assert all([a == b for a, b in zip(lines_result, lines_expected)])


def predicate(line):
    """
    Remove lines starting with `#`
    """
    if "#" in line:
        return False
    return True
