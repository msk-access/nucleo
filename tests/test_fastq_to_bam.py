#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `fastq_to_bam` package."""

import os
import subprocess
import shutil
import difflib
import json

RESULT_FILE_NAME = [
    "info.txt",
    "composite-umi-frequencies.txt",
    "merged_fastq_R1_umi-clipped.fastq.gz_trimming_report.txt",
    "merged_fastq_R2_umi-clipped.fastq.gz_trimming_report.txt",
    "test_standard_md.bam",
    "test_standard_md.bai",
    "test_standard_md.metrics",
    "test_standard_md.bed",
    "test_standard_md_abra_fm_bqsr.bam",
    "test_standard_md_abra_fm_bqsr.bai",
    "test_standard_md_abra_fm_bqsr_alignment_metrics.txt",
    "test_standard_md_abra_fm_bqsr-pileup-without-duplicates.txt",
    "test_standard_md_abra_fm_bqsr-intervals.txt",
    "test_standard_md_abra_fm_bqsr-intervals-without-duplicates.txt",
    "first-pass-insertions.txt",
    "first-pass-alt-alleles.txt",
    "second-pass-insertions.txt",
    "collapsed_R1_.fastq",
    "collapsed_R2_.fastq",
    "second-pass-alt-alleles.txt",
    "test_collapsed_R1_fastq.gz",
    "test_collapsed_R2_fastq.gz",
    "test_unfiltered.bed",
    "test_unfiltered_abra_fm.bam",
    "test_unfiltered_abra_fm.bai",
    "test_unfiltered_abra_fm_alignment_metrics.txt",
    "test_unfiltered_abra_fm-simplex.bam",
    "test_unfiltered_abra_fm-simplex.bai",
    "test_unfiltered_abra_fm-duplex.bam",
    "test_unfiltered_abra_fm-duplex.bai",
    "test_unfiltered_abra_fm-duplex_alignment_metrics.txt",
    "test_unfiltered_abra_fm-simplex_alignment_metrics.txt",
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
            "fastq_to_bam.cwl",
            "test_fastq_to_bam/test_input/inputs.json",
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
        shutil.rmtree("test_fastq_to_bam")
    except OSError as e:
        print("ERROR: cannot remove folder test_fastq_to_bam : %s" % (e))


def test_check_if_metrics_file_are_same():
    """
    General tests for checking if the metrics file is the same
    """
    print("\n### Check if files are the same from alignment metrics calculation ###\n")

    compare_picard_metrics_files(
        "test_standard_md_abra_fm_bqsr_alignment_metrics.txt",
        "test_fastq_to_bam/test_output/test_standard_md_abra_fm_bqsr_alignment_metrics.txt",
    )
    compare_picard_metrics_files(
        "test_unfiltered_abra_fm_alignment_metrics.txt",
        "test_fastq_to_bam/test_output/test_unfiltered_abra_fm_alignment_metrics.txt",
    )
    compare_picard_metrics_files(
        "test_unfiltered-abra_fm_simplex_alignment_metrics.txt",
        "test_fastq_to_bam/test_output/test_unfiltered_abra_fm-simplex_alignment_metrics.txt",
    )
    compare_picard_metrics_files(
        "test_unfiltered-abra_fm_duplex_alignment_metrics.txt",
        "test_fastq_to_bam/test_output/test_unfiltered_abra_fm-duplex_alignment_metrics.txt",
    )

    # Todo: info.txt, md metrics, trimming report


def test_output_json():
    """
    General tests for output json
    """
    assert os.path.exists(OUTPUT_JSON_FILENAME)
    OUTPUT_JSON = json.loads(open(OUTPUT_JSON_FILENAME, 'r').read())
    # Todo: use constant instead of magic number
    assert len(OUTPUT_JSON) == 26


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
