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
    'sbg:y': 908.96875
  - id: fastq2
    type: File
    'sbg:x': 0
    'sbg:y': 802.03125
  - id: reference
    type: File
    secondaryFiles:
      - ^.dict
      - .fai
      - .amb
      - .ann
      - .bwt
      - .pac
      - .index
      - .rbwt
      - .rpac
      - .rsa
      - .sa
    'sbg:x': 337.3125
    'sbg:y': 246.34375
  - id: known_sites_1
    type: File
    secondaryFiles:
      - .idx
    'sbg:x': 337.3125
    'sbg:y': 1250.78125
  - id: bed_file
    type: File
    'sbg:x': 0
    'sbg:y': 1229.78125
  - id: read_group_sample_name
    type: string
    'sbg:x': 337.3125
    'sbg:y': 353.28125
  - id: read_group_platform_unit
    type: string
    'sbg:x': 337.3125
    'sbg:y': 460.21875
  - id: read_group_library
    type: int
    'sbg:x': 337.3125
    'sbg:y': 567.15625
  - id: read_group_identifier
    type: string
    'sbg:x': 337.3125
    'sbg:y': 674.09375
  - id: sort_first_pass_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 481.21875
  - id: output_name_collapsed_gzip_R2
    type: string?
    'sbg:x': 0
    'sbg:y': 588.15625
  - id: output_name_collapsed_gzip_R1
    type: string?
    'sbg:x': 0
    'sbg:y': 695.09375
  - id: collapsing_aln_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 1122.84375
  - id: collapsing_picard_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 1015.90625
  - id: output
    type: string?
    label: standard_aln_output_file_name
    'sbg:x': 337.3125
    'sbg:y': 887.96875
  - id: output_file_name
    type: string?
    label: standard_picard_addrg_output_filename
    'sbg:x': 337.3125
    'sbg:y': 781.03125
  - id: known_sites_2
    type: File?
    secondaryFiles:
      - .idx
    'sbg:x': 337.3125
    'sbg:y': 1143.84375
  - id: adapter2
    type: string?
    'sbg:x': 337.3125
    'sbg:y': 1357.71875
  - id: adapter
    type: string?
    'sbg:x': 337.3125
    'sbg:y': 1464.65625
outputs:
  - id: composite_umi_frequencies
    outputSource:
      - marianas_process_loop_umi_cwl/composite_umi_frequencies
    type: File
    'sbg:x': 644.3773803710938
    'sbg:y': 855.5
  - id: clipping_info
    outputSource:
      - marianas_process_loop_umi_cwl/clipping_info
    type: File
    'sbg:x': 644.3773803710938
    'sbg:y': 962.4375
  - id: md_bam
    outputSource:
      - standard_bam_processing_cwl/md_bam
    type: File
    label: mark_duplicates_bam
    secondaryFiles:
      - ^.bai
    'sbg:x': 1068.72607421875
    'sbg:y': 522.625
  - id: clstats2
    outputSource:
      - standard_bam_processing_cwl/clstats2
    type: File
    label: trimming_stats_read2
    'sbg:x': 1068.72607421875
    'sbg:y': 629.5625
  - id: clstats1
    outputSource:
      - standard_bam_processing_cwl/clstats1
    type: File
    label: trimming_stats_read1
    'sbg:x': 1068.72607421875
    'sbg:y': 736.5
  - id: bqsr_bam
    outputSource:
      - standard_bam_processing_cwl/bqsr_bam
    type: File?
    label: standard_processed_bam
    secondaryFiles:
      - ^.bai
    'sbg:x': 1068.72607421875
    'sbg:y': 843.4375
  - id: unfiltered-bam
    outputSource:
      - bam_collapsing/unfiltered-bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1699.029052734375
    'sbg:y': 0
  - id: simplex-bam
    outputSource:
      - bam_collapsing/simplex-bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1699.029052734375
    'sbg:y': 106.9375
  - id: second_pass_insertions
    outputSource:
      - bam_collapsing/second_pass_insertions
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 213.875
  - id: second_pass_alt_alleles
    outputSource:
      - bam_collapsing/second_pass_alt_alleles
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 320.8125
  - id: pileup_without_duplicates
    outputSource:
      - bam_collapsing/pileup_without_duplicates
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 427.75
  - id: intervals_without_duplicates
    outputSource:
      - bam_collapsing/intervals_without_duplicates
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 534.6875
  - id: intervals
    outputSource:
      - bam_collapsing/intervals
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 641.625
  - id: gzip_read1
    outputSource:
      - bam_collapsing/gzip_read1
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 855.5
  - id: gzip_read2
    outputSource:
      - bam_collapsing/gzip_read2
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 748.5625
  - id: first_pass_insertions
    outputSource:
      - bam_collapsing/first_pass_insertions
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 962.4375
  - id: duplex-bam
    outputSource:
      - bam_collapsing/duplex-bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1699.029052734375
    'sbg:y': 1069.375
  - id: collapsed_fastq_2
    outputSource:
      - bam_collapsing/collapsed_fastq_2
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 1176.3125
  - id: alt_allele_file
    outputSource:
      - bam_collapsing/alt_allele_file
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 1390.1875
  - id: alignment_metrics_unfiltered
    outputSource:
      - bam_collapsing/alignment_metrics_unfiltered
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 1497.125
  - id: alignment_metrics_simplex
    outputSource:
      - bam_collapsing/alignment_metrics_simplex
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 1604.0625
  - id: alignment_metrics_duplex
    outputSource:
      - bam_collapsing/alignment_metrics_duplex
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 1711
  - id: collapsed_fastq_1
    outputSource:
      - bam_collapsing/collapsed_fastq_1
    type: File
    'sbg:x': 1699.029052734375
    'sbg:y': 1283.25
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
    'sbg:x': 337.3125
    'sbg:y': 1015.90625
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
      - id: window_size
        default: '800,700'
      - id: soft_clip_contig
        default: '100,30,80,15'
      - id: scoring_gap_alignments
        default: '8,32,48,1'
      - id: maximum_mixmatch_rate
        default: 0.1
      - id: maximum_average_depth
        default: 1000
      - id: contig_anchor
        default: '10,1'
      - id: consensus_sequence
        default: true
      - id: stringency
        default: 3
      - id: quality
        default: 1
      - id: length
        default: 25
      - id: adapter2
        source: adapter2
      - id: adapter
        source: adapter
      - id: bqsr_read_filter
        default: GoodCigarReadFilter
    out:
      - id: clstats2
      - id: clstats1
      - id: bqsr_bam
      - id: md_bam
      - id: output_file
    run: standard_bam_processing/standard_bam_processing.cwl
    label: standard_bam_processing.cwl
    'sbg:x': 644.3773803710938
    'sbg:y': 664.5625
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
    'sbg:x': 1068.72607421875
    'sbg:y': 1069.375
requirements:
  - class: SubworkflowFeatureRequirement
