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
    'sbg:y': 907.2421875
  - id: fastq2
    type: File
    'sbg:x': 0
    'sbg:y': 800.5078125
  - id: reference
    type: File
    'sbg:x': 337.359375
    'sbg:y': 405.9375
  - id: known_sites_1
    type: File
    'sbg:x': 337.359375
    'sbg:y': 1301.8125
  - id: bed_file
    type: File
    'sbg:x': 0
    'sbg:y': 1227.4453125
  - id: read_group_sample_name
    type: string
    'sbg:x': 337.359375
    'sbg:y': 512.671875
  - id: read_group_platform_unit
    type: string
    'sbg:x': 337.359375
    'sbg:y': 619.40625
  - id: read_group_library
    type: int
    'sbg:x': 337.359375
    'sbg:y': 726.140625
  - id: read_group_identifier
    type: string
    'sbg:x': 337.359375
    'sbg:y': 832.875
  - id: sort_first_pass_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 480.3046875
  - id: output_name_collapsed_gzip_R2
    type: string?
    'sbg:x': 0
    'sbg:y': 587.0390625
  - id: output_name_collapsed_gzip_R1
    type: string?
    'sbg:x': 0
    'sbg:y': 693.7734375
  - id: collapsing_aln_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 1120.7109375
  - id: collapsing_picard_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 1013.9765625
  - id: output
    type: string?
    label: standard_aln_output_file_name
    'sbg:x': 337.359375
    'sbg:y': 1046.34375
  - id: output_file_name
    type: string?
    label: standard_picard_addrg_output_filename
    'sbg:x': 337.359375
    'sbg:y': 939.609375
outputs:
  - id: composite_umi_frequencies
    outputSource:
      - marianas_process_loop_umi_cwl/composite_umi_frequencies
    type: File
    'sbg:x': 644.42431640625
    'sbg:y': 853.875
  - id: clipping_info
    outputSource:
      - marianas_process_loop_umi_cwl/clipping_info
    type: File
    'sbg:x': 644.42431640625
    'sbg:y': 960.609375
  - id: md_bam
    outputSource:
      - standard_bam_processing_cwl/md_bam
    type: File
    label: mark_duplicates_bam
    'sbg:x': 1028.618896484375
    'sbg:y': 521.2734375
  - id: clstats2
    outputSource:
      - standard_bam_processing_cwl/clstats2
    type: File
    label: trimming_stats_read2
    'sbg:x': 1028.618896484375
    'sbg:y': 628.0078125
  - id: clstats1
    outputSource:
      - standard_bam_processing_cwl/clstats1
    type: File
    label: trimming_stats_read1
    'sbg:x': 1028.618896484375
    'sbg:y': 734.7421875
  - id: bqsr_bam
    outputSource:
      - standard_bam_processing_cwl/bqsr_bam
    type: File?
    label: standard_processed_bam
    'sbg:x': 1028.618896484375
    'sbg:y': 841.4765625
  - id: unfiltered-bam
    outputSource:
      - bam_collapsing/unfiltered-bam
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 0
  - id: simplex-bam
    outputSource:
      - bam_collapsing/simplex-bam
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 106.734375
  - id: second_pass_insertions
    outputSource:
      - bam_collapsing/second_pass_insertions
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 213.46875
  - id: second_pass_alt_alleles
    outputSource:
      - bam_collapsing/second_pass_alt_alleles
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 320.203125
  - id: pileup_without_duplicates
    outputSource:
      - bam_collapsing/pileup_without_duplicates
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 426.9375
  - id: intervals_without_duplicates
    outputSource:
      - bam_collapsing/intervals_without_duplicates
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 533.671875
  - id: intervals
    outputSource:
      - bam_collapsing/intervals
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 640.40625
  - id: gzip_read1
    outputSource:
      - bam_collapsing/gzip_read1
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 853.875
  - id: gzip_read2
    outputSource:
      - bam_collapsing/gzip_read2
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 747.140625
  - id: first_pass_insertions
    outputSource:
      - bam_collapsing/first_pass_insertions
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 960.609375
  - id: duplex-bam
    outputSource:
      - bam_collapsing/duplex-bam
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 1067.34375
  - id: collapsed_fastq_2
    outputSource:
      - bam_collapsing/collapsed_fastq_2
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 1174.078125
  - id: alt_allele_file
    outputSource:
      - bam_collapsing/alt_allele_file
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 1387.546875
  - id: alignment_metrics_unfiltered
    outputSource:
      - bam_collapsing/alignment_metrics_unfiltered
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 1494.28125
  - id: alignment_metrics_simplex
    outputSource:
      - bam_collapsing/alignment_metrics_simplex
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 1601.015625
  - id: alignment_metrics_duplex
    outputSource:
      - bam_collapsing/alignment_metrics_duplex
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 1707.75
  - id: collapsed_fastq_1
    outputSource:
      - bam_collapsing/collapsed_fastq_1
    type: File
    'sbg:x': 1658.921875
    'sbg:y': 1280.8125
steps:
  - id: marianas_process_loop_umi_cwl
    in:
      - id: fastq1
        source: fastq1
      - id: fastq2
        source: fastq2
      - id: umi_length
        default: 3
    out:
      - id: processed_fastq_1
      - id: processed_fastq_2
      - id: clipping_info
      - id: composite_umi_frequencies
    run: >-
      command_line_tools/marianas_process_loop_umi_1.8.1/marianas_process_loop_umi.cwl
    label: marianas_process_loop_umi.cwl
    'sbg:x': 337.359375
    'sbg:y': 1174.078125
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
      - id: validation_stringency
        default: LENIENT
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
        default: MSKCC
      - id: read_group_sequencing_platform
        default: ILLUMINA
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
        default: coordinate
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
    'sbg:x': 644.42431640625
    'sbg:y': 684.0078125
  - id: bam_collapsing
    in:
      - id: reference_fasta
        source: reference
      - id: bed_file
        source: bed_file
      - id: bam
        source: standard_bam_processing_cwl/bqsr_bam
      - id: min_map_quality
        default: 1
      - id: min_base_quality
        default: 20
      - id: mismatches
        default: 0
      - id: wobble
        default: 1
      - id: min_consensus_percent
        default: 90
      - id: key
        default:
          - '6,6n'
          - '8,8n'
      - id: sort_first_pass_output_file_name
        source: sort_first_pass_output_file_name
      - id: output_name_collapsed_gzip_R1
        source: output_name_collapsed_gzip_R1
      - id: output_name_collapsed_gzip_R2
        source: output_name_collapsed_gzip_R2
      - id: read_group_sequnecing_center
        default: MSKCC
      - id: read_group_sequencing_platform
        default: ILLUMINA
      - id: read_group_sample_name
        source: read_group_sample_name
      - id: read_group_platform_unit
        source: read_group_platform_unit
      - id: read_group_library
        source: read_group_library
      - id: read_group_identifier
        source: read_group_identifier
      - id: picard_output_file_name
        source: collapsing_picard_output_file_name
      - id: aln_output_file_name
        source: collapsing_aln_output_file_name
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
    'sbg:x': 1028.618896484375
    'sbg:y': 1067.34375
requirements:
  - class: SubworkflowFeatureRequirement
