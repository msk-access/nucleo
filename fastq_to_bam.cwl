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
    'sbg:y': 1279.125
  - id: fastq2
    type: File
    'sbg:x': 0
    'sbg:y': 1172.546875
  - id: reference
    type: File
    'sbg:x': 317.109375
    'sbg:y': 458.1796875
  - id: known_sites_1
    type: File
    'sbg:x': 317.109375
    'sbg:y': 1353.5703125
  - id: bed_file
    type: File
    'sbg:x': 0
    'sbg:y': 1385.703125
  - id: read_group_sequencing_platform
    type: string
    'sbg:x': 317.109375
    'sbg:y': 671.3984375
  - id: read_group_sample_name
    type: string
    'sbg:x': 317.109375
    'sbg:y': 778.1328125
  - id: read_group_platform_unit
    type: string
    'sbg:x': 317.109375
    'sbg:y': 884.8671875
  - id: read_group_library
    type: int
    'sbg:x': 317.109375
    'sbg:y': 991.6015625
  - id: read_group_identifier
    type: string
    'sbg:x': 317.109375
    'sbg:y': 1098.3359375
  - id: sort_order
    type: string?
    'sbg:x': 317.109375
    'sbg:y': 351.7734375
  - id: wobble
    type: int
    'sbg:x': 0
    'sbg:y': 319.8125
  - id: sort_first_pass_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 426.21875
  - id: read_group_sequnecing_center
    type: string
    'sbg:x': 317.109375
    'sbg:y': 564.6640625
  - id: output_name_collapsed_gzip_R2
    type: string?
    'sbg:x': 0
    'sbg:y': 532.875
  - id: output_name_collapsed_gzip_R1
    type: string?
    'sbg:x': 0
    'sbg:y': 639.609375
  - id: mismatches
    type: int
    'sbg:x': 0
    'sbg:y': 746.09375
  - id: min_map_quality
    type: int
    'sbg:x': 0
    'sbg:y': 852.578125
  - id: min_consensus_percent
    type: int
    'sbg:x': 0
    'sbg:y': 959.234375
  - id: min_base_quality
    type: int
    'sbg:x': 0
    'sbg:y': 1065.890625
  - id: aln_output_file_name
    type: string?
    'sbg:x': 636.234375
    'sbg:y': 1012.6171875
outputs:
  - id: composite_umi_frequencies
    outputSource:
      - marianas_process_loop_umi_cwl/composite_umi_frequencies
    type: File
    'sbg:x': 636.234375
    'sbg:y': 799.3046875
  - id: clipping_info
    outputSource:
      - marianas_process_loop_umi_cwl/clipping_info
    type: File
    'sbg:x': 636.234375
    'sbg:y': 905.9609375
  - id: md_bam
    outputSource:
      - standard_bam_processing_cwl/md_bam
    type: File
    label: mark_duplicates_bam
    'sbg:x': 1072.93017578125
    'sbg:y': 520.28125
  - id: clstats2
    outputSource:
      - standard_bam_processing_cwl/clstats2
    type: File
    label: trimming_stats_read2
    'sbg:x': 1072.93017578125
    'sbg:y': 626.9375
  - id: clstats1
    outputSource:
      - standard_bam_processing_cwl/clstats1
    type: File
    label: trimming_stats_read1
    'sbg:x': 1072.93017578125
    'sbg:y': 733.671875
  - id: bqsr_bam
    outputSource:
      - standard_bam_processing_cwl/bqsr_bam
    type: File?
    label: standard_processed_bam
    'sbg:x': 1072.93017578125
    'sbg:y': 840.328125
  - id: unfiltered-bam
    outputSource:
      - bam_collapsing/unfiltered-bam
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 0
  - id: simplex-bam
    outputSource:
      - bam_collapsing/simplex-bam
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 106.40625
  - id: second_pass_insertions
    outputSource:
      - bam_collapsing/second_pass_insertions
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 212.984375
  - id: second_pass_alt_alleles
    outputSource:
      - bam_collapsing/second_pass_alt_alleles
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 319.5625
  - id: pileup_without_duplicates
    outputSource:
      - bam_collapsing/pileup_without_duplicates
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 426.140625
  - id: intervals_without_duplicates
    outputSource:
      - bam_collapsing/intervals_without_duplicates
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 532.71875
  - id: intervals
    outputSource:
      - bam_collapsing/intervals
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 639.125
  - id: gzip_read1
    outputSource:
      - bam_collapsing/gzip_read1
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 852.34375
  - id: gzip_read2
    outputSource:
      - bam_collapsing/gzip_read2
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 745.609375
  - id: first_pass_insertions
    outputSource:
      - bam_collapsing/first_pass_insertions
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 959
  - id: duplex-bam
    outputSource:
      - bam_collapsing/duplex-bam
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 1065.578125
  - id: collapsed_fastq_2
    outputSource:
      - bam_collapsing/collapsed_fastq_2
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 1172.15625
  - id: alt_allele_file
    outputSource:
      - bam_collapsing/alt_allele_file
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 1385.3125
  - id: alignment_metrics_unfiltered
    outputSource:
      - bam_collapsing/alignment_metrics_unfiltered
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 1491.96875
  - id: alignment_metrics_simplex
    outputSource:
      - bam_collapsing/alignment_metrics_simplex
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 1598.703125
  - id: alignment_metrics_duplex
    outputSource:
      - bam_collapsing/alignment_metrics_duplex
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 1705.4375
  - id: collapsed_fastq_1
    outputSource:
      - bam_collapsing/collapsed_fastq_1
    type: File
    'sbg:x': 1709.0699462890625
    'sbg:y': 1278.734375
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
    'sbg:x': 317.109375
    'sbg:y': 1225.9921875
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
    out:
      - id: clstats2
      - id: clstats1
      - id: bqsr_bam
      - id: md_bam
    run: standard_bam_processing/standard_bam_processing.cwl
    label: standard_bam_processing.cwl
    'sbg:x': 636.234375
    'sbg:y': 622.6484375
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
    'sbg:x': 1072.93017578125
    'sbg:y': 1065.984375
requirements:
  - class: SubworkflowFeatureRequirement
