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
    'sbg:y': 1335.9375
  - id: fastq2
    type: File
    'sbg:x': 0
    'sbg:y': 1229.0625
  - id: reference
    type: File
    'sbg:x': 317.09375
    'sbg:y': 353.0625
  - id: known_sites_1
    type: File
    'sbg:x': 317.09375
    'sbg:y': 1463.8125
  - id: bed_file
    type: File
    'sbg:x': 0
    'sbg:y': 1442.8125
  - id: read_group_sequencing_platform
    type: string
    'sbg:x': 317.09375
    'sbg:y': 566.8125
  - id: read_group_sample_name
    type: string
    'sbg:x': 317.09375
    'sbg:y': 673.6875
  - id: read_group_platform_unit
    type: string
    'sbg:x': 317.09375
    'sbg:y': 780.5625
  - id: read_group_library
    type: int
    'sbg:x': 317.09375
    'sbg:y': 887.4375
  - id: read_group_identifier
    type: string
    'sbg:x': 317.09375
    'sbg:y': 994.3125
  - id: sort_order
    type: string?
    'sbg:x': 317.09375
    'sbg:y': 246.1875
  - id: wobble
    type: int
    'sbg:x': 0
    'sbg:y': 267.1875
  - id: sort_first_pass_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 374.0625
  - id: read_group_sequnecing_center
    type: string
    'sbg:x': 317.09375
    'sbg:y': 459.9375
  - id: output_name_collapsed_gzip_R2
    type: string?
    'sbg:x': 0
    'sbg:y': 587.8125
  - id: output_name_collapsed_gzip_R1
    type: string?
    'sbg:x': 0
    'sbg:y': 694.6875
  - id: mismatches
    type: int
    'sbg:x': 0
    'sbg:y': 801.5625
  - id: min_map_quality
    type: int
    'sbg:x': 0
    'sbg:y': 908.4375
  - id: min_consensus_percent
    type: int
    'sbg:x': 0
    'sbg:y': 1015.3125
  - id: min_base_quality
    type: int
    'sbg:x': 0
    'sbg:y': 1122.1875
  - id: aln_output_file_name
    type: string?
    'sbg:x': 636.203125
    'sbg:y': 1015.3125
  - id: picard_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 480.9375
  - id: output
    type: string?
    label: standard_aln_output_file_name
    'sbg:x': 317.09375
    'sbg:y': 1208.0625
  - id: output_file_name
    type: string?
    label: standard_picard_addrg_output_filename
    'sbg:x': 317.09375
    'sbg:y': 1101.1875
outputs:
  - id: composite_umi_frequencies
    outputSource:
      - marianas_process_loop_umi_cwl/composite_umi_frequencies
    type: File
    'sbg:x': 636.203125
    'sbg:y': 801.5625
  - id: clipping_info
    outputSource:
      - marianas_process_loop_umi_cwl/clipping_info
    type: File
    'sbg:x': 636.203125
    'sbg:y': 908.4375
  - id: md_bam
    outputSource:
      - standard_bam_processing_cwl/md_bam
    type: File
    label: mark_duplicates_bam
    'sbg:x': 1093.02294921875
    'sbg:y': 515.25
  - id: clstats2
    outputSource:
      - standard_bam_processing_cwl/clstats2
    type: File
    label: trimming_stats_read2
    'sbg:x': 1093.02294921875
    'sbg:y': 622.125
  - id: clstats1
    outputSource:
      - standard_bam_processing_cwl/clstats1
    type: File
    label: trimming_stats_read1
    'sbg:x': 1093.02294921875
    'sbg:y': 729
  - id: bqsr_bam
    outputSource:
      - standard_bam_processing_cwl/bqsr_bam
    type: File?
    label: standard_processed_bam
    'sbg:x': 1093.02294921875
    'sbg:y': 835.875
  - id: unfiltered-bam
    outputSource:
      - bam_collapsing/unfiltered-bam
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 0
  - id: simplex-bam
    outputSource:
      - bam_collapsing/simplex-bam
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 106.875
  - id: second_pass_insertions
    outputSource:
      - bam_collapsing/second_pass_insertions
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 213.75
  - id: second_pass_alt_alleles
    outputSource:
      - bam_collapsing/second_pass_alt_alleles
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 320.625
  - id: pileup_without_duplicates
    outputSource:
      - bam_collapsing/pileup_without_duplicates
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 427.5
  - id: intervals_without_duplicates
    outputSource:
      - bam_collapsing/intervals_without_duplicates
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 534.375
  - id: intervals
    outputSource:
      - bam_collapsing/intervals
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 641.25
  - id: gzip_read1
    outputSource:
      - bam_collapsing/gzip_read1
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 855
  - id: gzip_read2
    outputSource:
      - bam_collapsing/gzip_read2
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 748.125
  - id: first_pass_insertions
    outputSource:
      - bam_collapsing/first_pass_insertions
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 961.875
  - id: duplex-bam
    outputSource:
      - bam_collapsing/duplex-bam
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 1068.75
  - id: collapsed_fastq_2
    outputSource:
      - bam_collapsing/collapsed_fastq_2
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 1175.625
  - id: alt_allele_file
    outputSource:
      - bam_collapsing/alt_allele_file
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 1389.375
  - id: alignment_metrics_unfiltered
    outputSource:
      - bam_collapsing/alignment_metrics_unfiltered
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 1496.25
  - id: alignment_metrics_simplex
    outputSource:
      - bam_collapsing/alignment_metrics_simplex
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 1603.125
  - id: alignment_metrics_duplex
    outputSource:
      - bam_collapsing/alignment_metrics_duplex
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 1710
  - id: collapsed_fastq_1
    outputSource:
      - bam_collapsing/collapsed_fastq_1
    type: File
    'sbg:x': 1741.572998046875
    'sbg:y': 1282.5
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
    'sbg:x': 317.09375
    'sbg:y': 1335.9375
  - id: standard_bam_processing_cwl
    in:
      - id: fastq2
        source: marianas_process_loop_umi_cwl/processed_fastq_1
      - id: reference
        source: reference
      - id: known_sites_1
        source: known_sites_1
      - id: paired
        default: true
      - id: gzip
        default: true
      - id: create_bam_index
        default: true
      - id: assume_sorted
        default: true
      - id: bam_index
        default: true
      - id: option_bedgraph
        default: true
      - id: fastq1
        source: marianas_process_loop_umi_cwl/processed_fastq_2
      - id: read_group_sequnecing_center
        source: read_group_sequnecing_center
      - id: read_group_sequencing_platform
        source: read_group_sequencing_platform
      - id: read_group_sample_name
        source: read_group_sample_name
      - id: read_group_platform_unit
        source: read_group_platform_unit
      - id: read_group_library
        source: read_group_library
      - id: read_group_identifier
        source: read_group_identifier
      - id: P
        default: true
      - id: M
        default: true
      - id: create_bam_index_1
        default: true
      - id: sort_order
        source: sort_order
      - id: output
        source: output
      - id: output_file_name
        source: output_file_name
    out:
      - id: clstats2
      - id: clstats1
      - id: bqsr_bam
      - id: md_bam
    run: standard_bam_processing/standard_bam_processing.cwl
    label: standard_bam_processing.cwl
    'sbg:x': 636.203125
    'sbg:y': 610.6875
  - id: bam_collapsing
    in:
      - id: reference_fasta
        source: reference
      - id: bed_file
        source: bed_file
      - id: bam
        source: standard_bam_processing_cwl/bqsr_bam
      - id: min_map_quality
        source: min_map_quality
      - id: min_base_quality
        source: min_base_quality
      - id: mismatches
        source: mismatches
      - id: wobble
        source: wobble
      - id: min_consensus_percent
        source: min_consensus_percent
      - id: sort_first_pass_output_file_name
        source: sort_first_pass_output_file_name
      - id: output_name_collapsed_gzip_R1
        source: output_name_collapsed_gzip_R1
      - id: output_name_collapsed_gzip_R2
        source: output_name_collapsed_gzip_R2
      - id: read_group_sequnecing_center
        source: read_group_sequnecing_center
      - id: read_group_sequencing_platform
        source: read_group_sequencing_platform
      - id: read_group_sample_name
        source: read_group_sample_name
      - id: read_group_platform_unit
        source: read_group_platform_unit
      - id: read_group_library
        source: read_group_library
      - id: read_group_identifier
        source: read_group_identifier
      - id: picard_output_file_name
        source: picard_output_file_name
      - id: aln_output_file_name
        source: aln_output_file_name
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
    'sbg:x': 1093.02294921875
    'sbg:y': 1068.75
requirements:
  - class: SubworkflowFeatureRequirement
