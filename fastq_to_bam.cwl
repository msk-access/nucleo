class: Workflow
cwlVersion: v1.0
id: fastq_to_bam
label: fastq_to_bam.cwl
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: fastq1
    type: File
    'sbg:x': 0
    'sbg:y': 855.625
  - id: fastq2
    type: File
    'sbg:x': 0
    'sbg:y': 748.671875
  - id: reference
    type: File
    'sbg:x': 124.96875
    'sbg:y': 727.671875
  - id: known_sites_1
    type: File
    'sbg:x': 124.96875
    'sbg:y': 983.578125
  - id: bed_file
    type: File
    'sbg:x': 0
    'sbg:y': 962.578125
outputs:
  - id: composite_umi_frequencies
    outputSource:
      - marianas_process_loop_umi_cwl/composite_umi_frequencies
    type: File
    'sbg:x': 432.03369140625
    'sbg:y': 855.625
  - id: clipping_info
    outputSource:
      - marianas_process_loop_umi_cwl/clipping_info
    type: File
    'sbg:x': 432.03369140625
    'sbg:y': 962.578125
  - id: md_bam
    outputSource:
      - standard_bam_processing_cwl/md_bam
    type: File
    'sbg:x': 705.62744140625
    'sbg:y': 522.71875
  - id: clstats2
    outputSource:
      - standard_bam_processing_cwl/clstats2
    type: File
    'sbg:x': 705.62744140625
    'sbg:y': 629.671875
  - id: clstats1
    outputSource:
      - standard_bam_processing_cwl/clstats1
    type: File
    'sbg:x': 705.62744140625
    'sbg:y': 736.625
  - id: bqsr_bam
    outputSource:
      - standard_bam_processing_cwl/bqsr_bam
    type: File?
    'sbg:x': 705.62744140625
    'sbg:y': 843.578125
  - id: unfiltered-bam
    outputSource:
      - bam_collapsing/unfiltered-bam
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 0
  - id: simplex-bam
    outputSource:
      - bam_collapsing/simplex-bam
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 106.953125
  - id: second_pass_insertions
    outputSource:
      - bam_collapsing/second_pass_insertions
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 213.90625
  - id: second_pass_alt_alleles
    outputSource:
      - bam_collapsing/second_pass_alt_alleles
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 320.859375
  - id: pileup_without_duplicates
    outputSource:
      - bam_collapsing/pileup_without_duplicates
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 427.8125
  - id: intervals_without_duplicates
    outputSource:
      - bam_collapsing/intervals_without_duplicates
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 534.765625
  - id: intervals
    outputSource:
      - bam_collapsing/intervals
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 641.71875
  - id: gzip_read1
    outputSource:
      - bam_collapsing/gzip_read1
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 855.625
  - id: gzip_read2
    outputSource:
      - bam_collapsing/gzip_read2
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 748.671875
  - id: first_pass_insertions
    outputSource:
      - bam_collapsing/first_pass_insertions
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 962.578125
  - id: duplex-bam
    outputSource:
      - bam_collapsing/duplex-bam
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 1069.53125
  - id: collapsed_fastq_2
    outputSource:
      - bam_collapsing/collapsed_fastq_2
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 1176.484375
  - id: collapsed_fastq_1
    outputSource:
      - bam_collapsing/collapsed_fastq_1
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 1283.4375
  - id: alt_allele_file
    outputSource:
      - bam_collapsing/alt_allele_file
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 1390.390625
  - id: alignment_metrics_unfiltered
    outputSource:
      - bam_collapsing/alignment_metrics_unfiltered
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 1497.34375
  - id: alignment_metrics_simplex
    outputSource:
      - bam_collapsing/alignment_metrics_simplex
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 1604.296875
  - id: alignment_metrics_duplex
    outputSource:
      - bam_collapsing/alignment_metrics_duplex
    type: File
    'sbg:x': 1221.5023193359375
    'sbg:y': 1711.25
steps:
  - id: marianas_process_loop_umi_cwl
    in:
      - id: fastq1
        source: fastq1
      - id: fastq2
        source: fastq2
    out:
      - id: processed_fastq_1
      - id: processed_fastq_2
      - id: clipping_info
      - id: composite_umi_frequencies
    run: >-
      command_line_tools/marianas_process_loop_umi_1.8.1/marianas_process_loop_umi.cwl
    label: marianas_process_loop_umi.cwl
    'sbg:x': 124.96875
    'sbg:y': 855.625
  - id: standard_bam_processing_cwl
    in:
      - id: fastq2
        source: marianas_process_loop_umi_cwl/processed_fastq_1
      - id: reference
        source: reference
      - id: known_sites_1
        source: known_sites_1
      - id: fastq1
        source: marianas_process_loop_umi_cwl/processed_fastq_2
    out:
      - id: clstats2
      - id: clstats1
      - id: bqsr_bam
      - id: md_bam
    run: standard_bam_processing/standard_bam_processing.cwl
    label: standard_bam_processing.cwl
    'sbg:x': 432.03369140625
    'sbg:y': 727.671875
  - id: bam_collapsing
    in:
      - id: reference_fasta
        source: reference
      - id: bed_file
        source: bed_file
      - id: bam
        source: standard_bam_processing_cwl/bqsr_bam
    out:
      - id: second_pass_insertions
      - id: second_pass_alt_alleles
      - id: collapsed_fastq_2
      - id: collapsed_fastq_1
      - id: pileup_without_duplicates
      - id: intervals_without_duplicates
      - id: intervals
      - id: first_pass_insertions
      - id: alt_allele_file
      - id: first_pass_output_dir
      - id: gzip_read1
      - id: gzip_read2
      - id: unfiltered-bam
      - id: simplex-bam
      - id: duplex-bam
      - id: alignment_metrics_unfiltered
      - id: alignment_metrics_simplex
      - id: alignment_metrics_duplex
    run: bam_collapsing/bam_collapsing.cwl
    label: bam_collapsing
    'sbg:x': 705.62744140625
    'sbg:y': 1069.53125
requirements:
  - class: SubworkflowFeatureRequirement
