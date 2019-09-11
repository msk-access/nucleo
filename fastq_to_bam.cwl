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
    'sbg:y': 906.1015625
  - id: fastq2
    type: File
    'sbg:x': 0
    'sbg:y': 799.5234375
  - id: reference
    type: File
    secondaryFiles:
      - ^.dict
      - .fai
      - .amb
      - .ann
      - .bwt
      - .fai
      - .pac
      - .index
      - .rbwt
      - .rpac
      - .rsa
      - .sa
    'sbg:x': 337.34375
    'sbg:y': 351.921875
  - id: known_sites_1
    type: File
    secondaryFiles:
      - .idx
    'sbg:x': 337.34375
    'sbg:y': 1353.8046875
  - id: bed_file
    type: File
    'sbg:x': 0
    'sbg:y': 1226.1484375
  - id: read_group_sample_name
    type: string
    'sbg:x': 337.34375
    'sbg:y': 458.4765625
  - id: read_group_platform_unit
    type: string
    'sbg:x': 337.34375
    'sbg:y': 565.2109375
  - id: read_group_library
    type: int
    'sbg:x': 337.34375
    'sbg:y': 671.9453125
  - id: read_group_identifier
    type: string
    'sbg:x': 337.34375
    'sbg:y': 778.6796875
  - id: sort_first_pass_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 479.4765625
  - id: output_name_collapsed_gzip_R2
    type: string?
    'sbg:x': 0
    'sbg:y': 586.1328125
  - id: output_name_collapsed_gzip_R1
    type: string?
    'sbg:x': 0
    'sbg:y': 692.8671875
  - id: collapsing_aln_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 1119.4921875
  - id: collapsing_picard_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 1012.7578125
  - id: output
    type: string?
    label: standard_aln_output_file_name
    'sbg:x': 337.34375
    'sbg:y': 992.0703125
  - id: output_file_name
    type: string?
    label: standard_picard_addrg_output_filename
    'sbg:x': 337.34375
    'sbg:y': 885.4140625
  - id: known_sites_2
    type: File?
    secondaryFiles:
      - .idx
    'sbg:x': 337.34375
    'sbg:y': 1247.2265625
outputs:
  - id: composite_umi_frequencies
    outputSource:
      - marianas_process_loop_umi_cwl/composite_umi_frequencies
    type: File
    'sbg:x': 644.4086303710938
    'sbg:y': 852.734375
  - id: clipping_info
    outputSource:
      - marianas_process_loop_umi_cwl/clipping_info
    type: File
    'sbg:x': 644.4086303710938
    'sbg:y': 959.390625
  - id: md_bam
    outputSource:
      - standard_bam_processing_cwl/md_bam
    type: File
    label: mark_duplicates_bam
    secondaryFiles:
      - ^.bai
    'sbg:x': 1043.011474609375
    'sbg:y': 520.421875
  - id: clstats2
    outputSource:
      - standard_bam_processing_cwl/clstats2
    type: File
    label: trimming_stats_read2
    'sbg:x': 1043.011474609375
    'sbg:y': 627.078125
  - id: clstats1
    outputSource:
      - standard_bam_processing_cwl/clstats1
    type: File
    label: trimming_stats_read1
    'sbg:x': 1043.011474609375
    'sbg:y': 733.8125
  - id: bqsr_bam
    outputSource:
      - standard_bam_processing_cwl/bqsr_bam
    type: File?
    label: standard_processed_bam
    secondaryFiles:
      - ^.bai
    'sbg:x': 1043.011474609375
    'sbg:y': 840.46875
  - id: unfiltered-bam
    outputSource:
      - bam_collapsing/unfiltered-bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1673.314453125
    'sbg:y': 0
  - id: simplex-bam
    outputSource:
      - bam_collapsing/simplex-bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1673.314453125
    'sbg:y': 106.4765625
  - id: second_pass_insertions
    outputSource:
      - bam_collapsing/second_pass_insertions
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 213.0546875
  - id: second_pass_alt_alleles
    outputSource:
      - bam_collapsing/second_pass_alt_alleles
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 319.6328125
  - id: pileup_without_duplicates
    outputSource:
      - bam_collapsing/pileup_without_duplicates
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 426.2109375
  - id: intervals_without_duplicates
    outputSource:
      - bam_collapsing/intervals_without_duplicates
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 532.7890625
  - id: intervals
    outputSource:
      - bam_collapsing/intervals
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 639.265625
  - id: gzip_read1
    outputSource:
      - bam_collapsing/gzip_read1
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 852.5546875
  - id: gzip_read2
    outputSource:
      - bam_collapsing/gzip_read2
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 745.8203125
  - id: first_pass_insertions
    outputSource:
      - bam_collapsing/first_pass_insertions
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 959.2109375
  - id: duplex-bam
    outputSource:
      - bam_collapsing/duplex-bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1673.314453125
    'sbg:y': 1065.7890625
  - id: collapsed_fastq_2
    outputSource:
      - bam_collapsing/collapsed_fastq_2
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 1172.3671875
  - id: alt_allele_file
    outputSource:
      - bam_collapsing/alt_allele_file
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 1385.5234375
  - id: alignment_metrics_unfiltered
    outputSource:
      - bam_collapsing/alignment_metrics_unfiltered
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 1492.1796875
  - id: alignment_metrics_simplex
    outputSource:
      - bam_collapsing/alignment_metrics_simplex
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 1598.9140625
  - id: alignment_metrics_duplex
    outputSource:
      - bam_collapsing/alignment_metrics_duplex
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 1705.6484375
  - id: collapsed_fastq_1
    outputSource:
      - bam_collapsing/collapsed_fastq_1
    type: File
    'sbg:x': 1673.314453125
    'sbg:y': 1278.9453125
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
    'sbg:x': 337.34375
    'sbg:y': 1119.6484375
  - id: standard_bam_processing_cwl
    in:
      - id: fastq2
        source: marianas_process_loop_umi_cwl/processed_fastq_1
      - id: reference
        source: reference
      - id: known_sites_1
        source: known_sites_1
      - id: known_sites_2
        source: known_sites_2
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
    'sbg:x': 644.4086303710938
    'sbg:y': 675.9453125
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
    'sbg:x': 1043.011474609375
    'sbg:y': 1066.125
requirements:
  - class: SubworkflowFeatureRequirement
