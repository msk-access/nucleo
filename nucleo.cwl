class: Workflow
cwlVersion: v1.0
id: nucleo
doc: >-
  This workflow takes a READ1 and READ2 fastq.gz file generated for MSK-ACCESS
  assay and generated four different Binary Alignment Map file along with
  alignment metrics for each.
label: nucleo
$namespaces:
  s: 'https://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: reference_sequence
    type: File
    doc: >-
      Reference sequence file.  Please include ".fai", "^.dict", ".amb" , ".sa",
      ".bwt", ".pac", ".ann" as  secondary files if they are not present in the
      same location as the  ".fasta" file
    secondaryFiles:
      - .amb
      - .fai
      - .sa
      - ^.dict
      - .ann
      - .bwt
      - .pac
    'sbg:x': 0
    'sbg:y': 1389.375
  - id: gatk_base_recalibrator_known_sites
    type:
      type: array
      items: File
      inputBinding:
        prefix: '--known-sites'
    secondaryFiles:
      - .idx
    'sbg:x': 0
    'sbg:y': 2992.5
  - id: sequencing-center
    type: string?
    'sbg:x': 0
    'sbg:y': 1068.75
  - id: run-date
    type: string?
    'sbg:x': 0
    'sbg:y': 1282.5
  - id: sample
    type: string
    'sbg:x': 0
    'sbg:y': 1175.625
  - id: read-structures
    type: 'string[]?'
    'sbg:x': 0
    'sbg:y': 1496.25
  - id: read-group-id
    type: string
    'sbg:x': 0
    'sbg:y': 1603.125
  - id: platform-unit
    type: string
    'sbg:x': 0
    'sbg:y': 1710
  - id: platform-model
    type: string?
    'sbg:x': 0
    'sbg:y': 1816.875
  - id: platform
    type: string?
    'sbg:x': 0
    'sbg:y': 1923.75
  - id: library
    type: string
    'sbg:x': 0
    'sbg:y': 2351.25
  - id: validation_stringency
    type: string?
    'sbg:x': 0
    'sbg:y': 0
  - id: UBG_picard_SamToFastq_R1_output_fastq
    type: string
    'sbg:x': 0
    'sbg:y': 213.75
  - id: UBG_picard_SamToFastq_R2_output_fastq
    type: string
    'sbg:x': 0
    'sbg:y': 106.875
  - id: fastp_read2_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 5343.75
  - id: fastp_read2_adapter_sequence
    type: string?
    'sbg:x': 0
    'sbg:y': 5450.625
  - id: fastp_read1_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 5557.5
  - id: fastp_read1_adapter_sequence
    type: string?
    'sbg:x': 0
    'sbg:y': 5664.375
  - id: fastp_minimum_read_length
    type: int?
    'sbg:x': 0
    'sbg:y': 5771.25
  - id: fastp_html_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 6198.75
  - id: fastp_json_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 6091.875
  - id: sort_order
    type: string?
    'sbg:x': 0
    'sbg:y': 961.875
  - id: bwa_mem_T
    type: int?
    'sbg:x': 0
    'sbg:y': 6733.125
  - id: bwa_mem_Y
    type: boolean?
    'sbg:x': 0
    'sbg:y': 6626.25
  - id: UBG_picard_addRG_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 427.5
  - id: create_bam_index
    type: boolean?
    'sbg:x': 0
    'sbg:y': 6519.375
  - id: UBG_bwa_mem_output
    type: string
    'sbg:x': 0
    'sbg:y': 641.25
  - id: bwa_mem_K
    type: int?
    'sbg:x': 0
    'sbg:y': 6840
  - id: UBG_gatk_merge_bam_alignment_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 534.375
  - id: optical_duplicate_pixel_distance
    type: int?
    'sbg:x': 0
    'sbg:y': 2137.5
  - id: gatk_mark_duplicates_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 2565
  - id: gatk_mark_duplicates_duplication_metrics_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 2671.875
  - id: bedtools_merge_distance_between_features
    type: int?
    'sbg:x': 0
    'sbg:y': 6946.875
  - id: bedtools_genomecov_option_bedgraph
    type: boolean?
    'sbg:x': 0
    'sbg:y': 7053.75
  - id: apply_bqsr_output_file_name
    type: string
    'sbg:x': 611.53125
    'sbg:y': 4542.1875
  - id: abra2_window_size
    type: string?
    'sbg:x': 0
    'sbg:y': 8015.625
  - id: abra2_soft_clip_contig
    type: string?
    'sbg:x': 0
    'sbg:y': 8122.5
  - id: abra2_scoring_gap_alignments
    type: string?
    'sbg:x': 0
    'sbg:y': 8229.375
  - id: UBG_abra2_output_bams
    type:
      - string
      - type: array
        items: string
    'sbg:x': 0
    'sbg:y': 748.125
  - id: abra2_no_edge_complex_indel
    type: boolean?
    'sbg:x': 0
    'sbg:y': 8443.125
  - id: abra2_maximum_mixmatch_rate
    type: float?
    'sbg:x': 0
    'sbg:y': 8550
  - id: abra2_maximum_average_depth
    type: int?
    'sbg:x': 0
    'sbg:y': 8656.875
  - id: abra2_consensus_sequence
    type: boolean?
    'sbg:x': 0
    'sbg:y': 8870.625
  - id: abra2_bam_index
    type: boolean?
    'sbg:x': 0
    'sbg:y': 8977.5
  - id: fgbio_fastq_to_bam_input
    type:
      type: array
      items:
        items: File
        type: array
    'sbg:x': 0
    'sbg:y': 4595.625
  - id: UBG_picard_fixmateinformation_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 320.625
  - id: merge_sam_files_sort_order
    type: string?
    'sbg:x': 0
    'sbg:y': 2244.375
  - id: gatk_merge_sam_files_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 2458.125
  - id: fgbio_collect_duplex_seq_metrics_intervals
    type: File?
    'sbg:x': 0
    'sbg:y': 4809.375
  - id: fgbio_group_reads_by_umi_strategy
    type: string?
    'sbg:x': 0
    'sbg:y': 3313.125
  - id: fgbio_group_reads_by_umi_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 3420
  - id: fgbio_group_reads_by_umi_family_size_histogram
    type: string
    'sbg:x': 0
    'sbg:y': 3526.875
  - id: fgbio_collect_duplex_seq_metrics_output_prefix
    type: string?
    'sbg:x': 0
    'sbg:y': 4702.5
  - id: fgbio_collect_duplex_seq_metrics_duplex_umi_counts
    type: boolean?
    'sbg:x': 0
    'sbg:y': 4916.25
  - id: fgbio_call_duplex_consensus_reads_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 5023.125
  - id: fgbio_call_duplex_consensus_reads_min_reads
    type: 'int[]?'
    'sbg:x': 0
    'sbg:y': 5130
  - id: BC_gatk_sam_to_fastq_output_name_R2
    type: string
    'sbg:x': 0
    'sbg:y': 7374.375
  - id: BC_gatk_sam_to_fastq_output_name_R1
    type: string
    'sbg:x': 0
    'sbg:y': 7481.25
  - id: BC_picard_fixmate_information_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 7160.625
  - id: BC_picard_addRG_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 7267.5
  - id: BC_gatk_merge_bam_alignment_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 7588.125
  - id: fgbio_postprocessing_output_file_name_simplex
    type: string
    'sbg:x': 0
    'sbg:y': 3206.25
  - id: gatk_collect_alignment_summary_metrics_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 2885.625
  - id: abra2_contig_anchor
    type: string?
    'sbg:x': 0
    'sbg:y': 8763.75
  - id: BC_abra2_output_bams
    type:
      - string
      - type: array
        items: string
    'sbg:x': 0
    'sbg:y': 7801.875
  - id: fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex
    type: boolean?
    'sbg:x': 0
    'sbg:y': 3633.75
  - id: fgbio_filter_consensus_read_output_file_name_simplex_duplex
    type: string
    'sbg:x': 0
    'sbg:y': 3740.625
  - id: fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics
    type: string
    'sbg:x': 0
    'sbg:y': 3847.5
  - id: fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics
    type: string?
    'sbg:x': 0
    'sbg:y': 3954.375
  - id: fgbio_filter_consensus_read_output_file_name_duplex
    type: string?
    'sbg:x': 0
    'sbg:y': 4061.25
  - id: fgbio_filter_consensus_read_min_reads_duplex
    type: 'int[]?'
    'sbg:x': 0
    'sbg:y': 4275
  - id: fgbio_filter_consensus_read_min_base_quality_simplex_duplex
    type: int?
    'sbg:x': 0
    'sbg:y': 4381.875
  - id: fgbio_filter_consensus_read_min_base_quality_duplex
    type: int?
    'sbg:x': 0
    'sbg:y': 4488.75
  - id: BC_bwa_mem_output
    type: string
    'sbg:x': 0
    'sbg:y': 7695
  - id: fgbio_filter_consensus_read_min_reads_simplex_duplex
    type: 'int[]?'
    'sbg:x': 0
    'sbg:y': 4168.125
  - id: picard_addRG_sort_order
    type: string?
    'sbg:x': 0
    'sbg:y': 2030.625
  - id: disable_trim_poly_g
    type: boolean?
    'sbg:x': 0
    'sbg:y': 6305.625
  - id: disable_quality_filtering
    type: boolean?
    'sbg:x': 0
    'sbg:y': 6412.5
  - id: gatk_base_recalibrator_add_output_sam_program_record
    type: boolean?
    'sbg:x': 0
    'sbg:y': 3099.375
  - id: gatk_collect_aln_summary_metrics_bqsr_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 2778.75
  - id: abra2_no_sort
    type: boolean?
    'sbg:x': 0
    'sbg:y': 8336.25
  - id: base_recalibrator_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 7908.75
  - id: temporary_directory
    type: string?
    'sbg:x': 0
    'sbg:y': 855
  - id: fgbio_async_io
    type: string?
    'sbg:x': 0
    'sbg:y': 5236.875
  - id: fastp_max_len_read2
    type: int?
    'sbg:x': 0
    'sbg:y': 5878.125
  - id: fastp_max_len_read1
    type: int?
    'sbg:x': 0
    'sbg:y': 5985
outputs:
  - id: fastp_html_output
    outputSource:
      - uncollapsed_bam_generation/fastp_html_output
    type: File
    'sbg:x': 1896.48779296875
    'sbg:y': 4108.1875
  - id: fastp_json_output
    outputSource:
      - uncollapsed_bam_generation/fastp_json_output
    type: File
    'sbg:x': 1896.48779296875
    'sbg:y': 4001.3125
  - id: gatk_collect_alignment_summary_metrics_txt_uncollapsed
    outputSource:
      - >-
        gatk_collect_alignment_summary_metrics_4_1_8_0/gatk_collect_alignment_summary_metrics_txt
    type: File
    'sbg:x': 3899.757080078125
    'sbg:y': 4488.75
  - id: indel_realignment_bam
    outputSource:
      - uncollapsed_bam_generation/indel_realignment_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1896.48779296875
    'sbg:y': 3894.4375
  - id: picard_mark_duplicates_metrics
    outputSource:
      - uncollapsed_bam_generation/picard_mark_duplicates_metrics
    type: File
    'sbg:x': 1896.48779296875
    'sbg:y': 3787.5625
  - id: gatk_collect_alignment_summary_metrics_txt_simplex
    outputSource:
      - bam_collapsing/gatk_collect_alignment_summary_metrics_txt_simplex
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 3773.0625
  - id: gatk_collect_alignment_summary_metrics_txt_duplex
    outputSource:
      - bam_collapsing/gatk_collect_alignment_summary_metrics_txt_duplex
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 3879.9375
  - id: gatk_collect_alignment_summary_metrics_txt_collapsed
    outputSource:
      - bam_collapsing/gatk_collect_alignment_summary_metrics_txt_collapsed
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 3986.8125
  - id: fgbio_postprocessing_simplex_bam
    outputSource:
      - bam_collapsing/fgbio_postprocessing_simplex_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3383.319580078125
    'sbg:y': 4242.5625
  - id: fgbio_group_reads_by_umi_histogram
    outputSource:
      - bam_collapsing/fgbio_group_reads_by_umi_histogram
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 4349.4375
  - id: fgbio_group_reads_by_umi_bam
    outputSource:
      - bam_collapsing/fgbio_group_reads_by_umi_bam
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 4456.3125
  - id: fgbio_filter_consensus_reads_duplex_bam
    outputSource:
      - bam_collapsing/fgbio_filter_consensus_reads_duplex_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3383.319580078125
    'sbg:y': 4563.1875
  - id: fgbio_collect_duplex_seq_metrics_umi_counts
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_umi_counts
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 4670.0625
  - id: fgbio_collect_duplex_seq_metrics_family_size
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_family_size
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 4776.9375
  - id: fgbio_collect_duplex_seq_metrics_duplex_yield_metrics
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_yield_metrics
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 4883.8125
  - id: fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 4990.6875
  - id: fgbio_collect_duplex_seq_metrics_duplex_qc
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_qc
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 5097.5625
  - id: fgbio_collect_duplex_seq_metrics_duplex_family_size
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_family_size
    type: File
    'sbg:x': 3383.319580078125
    'sbg:y': 5204.4375
  - id: fgbio_collapsed_bam
    outputSource:
      - bam_collapsing/fgbio_collapsed_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3383.319580078125
    'sbg:y': 5311.3125
  - id: uncollapsed_bam
    outputSource:
      - base_quality_recalibration/gatk_apply_bqsr_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3383.319580078125
    'sbg:y': 3666.1875
steps:
  - id: uncollapsed_bam_generation
    in:
      - id: sequencing-center
        default: MSKCC
        source: sequencing-center
      - id: sample
        source: sample
      - id: run-date
        source: run-date
      - id: read-structures
        default:
          - 3M2S+T
          - 3M2S+T
        source:
          - read-structures
      - id: read-group-id
        source: read-group-id
      - id: platform-unit
        source: platform-unit
      - id: platform-model
        default: novaseq
        source: platform-model
      - id: platform
        default: ILLUMINA
        source: platform
      - id: library
        source: library
      - id: validation_stringency
        default: LENIENT
        source: validation_stringency
      - id: R1_output_fastq
        source: UBG_picard_SamToFastq_R1_output_fastq
      - id: R2_output_fastq
        source: UBG_picard_SamToFastq_R2_output_fastq
      - id: reference_sequence
        source: reference_sequence
      - id: fastp_read2_output_file_name
        source: fastp_read2_output_file_name
      - id: fastp_read2_adapter_sequence
        default: AGATCGGAAGAGC
        source: fastp_read2_adapter_sequence
      - id: fastp_read1_output_file_name
        source: fastp_read1_output_file_name
      - id: fastp_read1_adapter_sequence
        default: GATCGGAAGAGC
        source: fastp_read1_adapter_sequence
      - id: fastp_minimum_read_length
        default: 25
        source: fastp_minimum_read_length
      - id: fastp_json_output_file_name
        source: fastp_json_output_file_name
      - id: fastp_html_output_file_name
        source: fastp_html_output_file_name
      - id: bwa_mem_Y
        default: true
        source: bwa_mem_Y
      - id: bwa_mem_T
        source: bwa_mem_T
      - id: sort_order
        default: coordinate
        source: sort_order
      - id: picard_addRG_output_file_name
        source: UBG_picard_addRG_output_file_name
      - id: bwa_mem_output
        source: UBG_bwa_mem_output
      - id: bwa_mem_K
        source: bwa_mem_K
      - id: create_bam_index
        default: true
        source: create_bam_index
      - id: gatk_merge_bam_alignment_output_file_name
        source: UBG_gatk_merge_bam_alignment_output_file_name
      - id: optical_duplicate_pixel_distance
        default: 2500
        source: optical_duplicate_pixel_distance
      - id: gatk_mark_duplicates_output_file_name
        source: gatk_mark_duplicates_output_file_name
      - id: gatk_mark_duplicates_duplication_metrics_file_name
        source: gatk_mark_duplicates_duplication_metrics_file_name
      - id: abra2_window_size
        default: '800,700'
        source: abra2_window_size
      - id: abra2_soft_clip_contig
        default: '100,30,80,15'
        source: abra2_soft_clip_contig
      - id: abra2_scoring_gap_alignments
        default: '8,32,48,1'
        source: abra2_scoring_gap_alignments
      - id: abra2_output_bams
        source:
          - UBG_abra2_output_bams
      - id: abra2_maximum_average_depth
        default: 1000
        source: abra2_maximum_average_depth
      - id: abra2_bam_index
        source: abra2_bam_index
      - id: abra2_contig_anchor
        default: '10,1'
        source: abra2_contig_anchor
      - id: abra2_consensus_sequence
        default: true
        source: abra2_consensus_sequence
      - id: bedtools_merge_distance_between_features
        default: 10
        source: bedtools_merge_distance_between_features
      - id: abra2_maximum_mixmatch_rate
        default: 0.1
        source: abra2_maximum_mixmatch_rate
      - id: bedtools_genomecov_option_bedgraph
        default: true
        source: bedtools_genomecov_option_bedgraph
      - id: picard_fixmateinformation_output_file_name
        source: UBG_picard_fixmateinformation_output_file_name
      - id: abra2_no_sort
        default: true
        source: abra2_no_sort
      - id: abra2_no_edge_complex_indel
        default: true
        source: abra2_no_edge_complex_indel
      - id: merge_sam_files_sort_order
        default: queryname
        source: merge_sam_files_sort_order
      - id: gatk_merge_sam_files_output_file_name
        source: gatk_merge_sam_files_output_file_name
      - id: fgbio_fastq_to_bam_input
        source:
          - fgbio_fastq_to_bam_input
      - id: picard_addRG_sort_order
        default: queryname
        source: picard_addRG_sort_order
      - id: disable_trim_poly_g
        default: true
        source: disable_trim_poly_g
      - id: disable_quality_filtering
        default: true
        source: disable_quality_filtering
      - id: temporary_directory
        source: temporary_directory
      - id: fgbio_async_io
        default: 'true'
        source: fgbio_async_io
      - id: fastp_max_len_read1
        default: 101
        source: fastp_max_len_read1
      - id: fastp_max_len_read2
        default: 101
        source: fastp_max_len_read2
    out:
      - id: gatk_sam_to_fastq_unpaired_fastq
      - id: fastp_unpaired2_output
      - id: fastp_unpaired1_output
      - id: fastp_json_output
      - id: fastp_html_output
      - id: picard_mark_duplicates_metrics
      - id: indel_realignment_bam
    run: uncollapsed_bam_generation/uncollapsed_bam_generation.cwl
    label: Uncollapsed BAM Generation
    'sbg:x': 611.53125
    'sbg:y': 4057.3125
  - id: bam_collapsing
    in:
      - id: fgbio_group_reads_by_umi_input
        source: uncollapsed_bam_generation/indel_realignment_bam
      - id: fgbio_group_reads_by_umi_strategy
        default: paired
        source: fgbio_group_reads_by_umi_strategy
      - id: fgbio_group_reads_by_umi_output_file_name
        source: fgbio_group_reads_by_umi_output_file_name
      - id: fgbio_group_reads_by_umi_family_size_histogram
        source: fgbio_group_reads_by_umi_family_size_histogram
      - id: fgbio_collect_duplex_seq_metrics_intervals
        source: fgbio_collect_duplex_seq_metrics_intervals
      - id: fgbio_collect_duplex_seq_metrics_output_prefix
        source: fgbio_collect_duplex_seq_metrics_output_prefix
      - id: fgbio_collect_duplex_seq_metrics_duplex_umi_counts
        default: true
        source: fgbio_collect_duplex_seq_metrics_duplex_umi_counts
      - id: fgbio_call_duplex_consensus_reads_read_group_id
        source: read-group-id
      - id: fgbio_call_duplex_consensus_reads_output_file_name
        source: fgbio_call_duplex_consensus_reads_output_file_name
      - id: fgbio_call_duplex_consensus_reads_min_reads
        default:
          - 1
          - 1
          - 0
        source:
          - fgbio_call_duplex_consensus_reads_min_reads
      - id: reference_sequence
        source: reference_sequence
      - id: validation_stringency
        default: LENIENT
        source: validation_stringency
      - id: gatk_sam_to_fastq_output_name_R2
        source: BC_gatk_sam_to_fastq_output_name_R2
      - id: gatk_sam_to_fastq_output_name_R1
        source: BC_gatk_sam_to_fastq_output_name_R1
      - id: bwa_mem_Y
        default: true
        source: bwa_mem_Y
      - id: bwa_mem_T
        source: bwa_mem_T
      - id: sort_order
        default: coordinate
        source: sort_order
      - id: picard_addRG_read_group_sequencing_platform
        default: ILLUMINA
        source: platform
      - id: picard_addRG_read_group_sequencing_center
        default: MSKCC
        source: sequencing-center
      - id: picard_addRG_read_group_run_date
        source: run-date
      - id: picard_addRG_read_group_platform_unit
        source: platform-unit
      - id: picard_addRG_read_group_library
        source: library
      - id: picard_addRG_read_group_identifier
        source: read-group-id
      - id: picard_addRG_output_file_name
        source: BC_picard_addRG_output_file_name
      - id: bwa_mem_output
        source: BC_bwa_mem_output
      - id: create_bam_index
        default: true
        source: create_bam_index
      - id: bwa_mem_K
        source: bwa_mem_K
      - id: abra2_window_size
        default: '800,700'
        source: abra2_window_size
      - id: abra2_soft_clip_contig
        default: '100,30,80,15'
        source: abra2_soft_clip_contig
      - id: abra2_scoring_gap_alignments
        default: '8,32,48,1'
        source: abra2_scoring_gap_alignments
      - id: picard_fixmate_information_output_file_name
        source: BC_picard_fixmate_information_output_file_name
      - id: abra2_output_bams
        source:
          - BC_abra2_output_bams
      - id: bedtools_genomecov_option_bedgraph
        default: true
        source: bedtools_genomecov_option_bedgraph
      - id: abra2_no_sort
        default: true
        source: abra2_no_sort
      - id: abra2_no_edge_complex_indel
        default: true
        source: abra2_no_edge_complex_indel
      - id: abra2_maximum_mixmatch_rate
        default: 0.1
        source: abra2_maximum_mixmatch_rate
      - id: abra2_maximum_average_depth
        default: 1000
        source: abra2_maximum_average_depth
      - id: bedtools_merge_distance_between_features
        default: 10
        source: bedtools_merge_distance_between_features
      - id: abra2_contig_anchor
        default: '10,1'
        source: abra2_contig_anchor
      - id: abra2_consensus_sequence
        default: true
        source: abra2_consensus_sequence
      - id: gatk_merge_bam_alignment_output_file_name
        source: BC_gatk_merge_bam_alignment_output_file_name
      - id: fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex
        default: true
        source: fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex
      - id: fgbio_filter_consensus_read_reverse_per_base_tags_duplex
        default: true
      - id: fgbio_filter_consensus_read_min_base_quality_duplex
        default: 30
        source: fgbio_filter_consensus_read_min_base_quality_duplex
      - id: fgbio_filter_consensus_read_min_base_quality_simplex_duplex
        default: 30
        source: fgbio_filter_consensus_read_min_base_quality_simplex_duplex
      - id: fgbio_filter_consensus_read_min_reads_duplex
        default:
          - 2
          - 1
          - 1
        source:
          - fgbio_filter_consensus_read_min_reads_duplex
      - id: fgbio_filter_consensus_read_min_reads_simplex_duplex
        default:
          - 3
          - 3
          - 0
        source:
          - fgbio_filter_consensus_read_min_reads_simplex_duplex
      - id: fgbio_filter_consensus_read_output_file_name_simplex_duplex
        source: fgbio_filter_consensus_read_output_file_name_simplex_duplex
      - id: fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics
        source: fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics
      - id: fgbio_postprocessing_output_file_name_simplex
        source: fgbio_postprocessing_output_file_name_simplex
      - id: fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics
        source: fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics
      - id: fgbio_filter_consensus_read_output_file_name_duplex
        source: fgbio_filter_consensus_read_output_file_name_duplex
      - id: picard_addRG_read_group_sample_name
        source: sample
      - id: gatk_collect_alignment_summary_metrics_output_file_name
        source: gatk_collect_alignment_summary_metrics_output_file_name
      - id: picard_addRG_sort_order
        default: queryname
        source: picard_addRG_sort_order
      - id: temporary_directory
        source: temporary_directory
      - id: async_io
        default: 'true'
        source: fgbio_async_io
    out:
      - id: fgbio_group_reads_by_umi_histogram
      - id: fgbio_collect_duplex_seq_metrics_umi_counts
      - id: fgbio_collect_duplex_seq_metrics_family_size
      - id: fgbio_collect_duplex_seq_metrics_duplex_yield_metrics
      - id: fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt
      - id: fgbio_collect_duplex_seq_metrics_duplex_qc
      - id: fgbio_collect_duplex_seq_metrics_duplex_family_size
      - id: gatk_sam_to_fastq_unpaired_fastq
      - id: gatk_sam_to_fastq_second_end_fastq
      - id: gatk_sam_to_fastq_fastq
      - id: gatk_collect_alignment_summary_metrics_txt_simplex
      - id: gatk_collect_alignment_summary_metrics_txt_duplex
      - id: fgbio_postprocessing_simplex_bam
      - id: fgbio_filter_consensus_reads_duplex_bam
      - id: fgbio_collapsed_bam
      - id: gatk_collect_alignment_summary_metrics_txt_collapsed
      - id: fgbio_group_reads_by_umi_bam
      - id: fgbio_filter_consensus_reads_simplex_duplex_bam
    run: bam_collapsing/bam_collapsing.cwl
    label: bam_collapsing
    'sbg:x': 1896.48779296875
    'sbg:y': 4804.9375
  - id: base_quality_recalibration
    in:
      - id: input
        source: uncollapsed_bam_generation/indel_realignment_bam
      - id: reference
        source: reference_sequence
      - id: known_sites
        source:
          - gatk_base_recalibrator_known_sites
      - id: base_recalibrator_output_file_name
        source: base_recalibrator_output_file_name
      - id: add_output_sam_program_record
        default: true
        source: gatk_base_recalibrator_add_output_sam_program_record
      - id: apply_bqsr_create_output_bam_index
        source: create_bam_index
      - id: apply_bqsr_output_file_name
        source: apply_bqsr_output_file_name
      - id: temporary_directory
        source: temporary_directory
    out:
      - id: gatk_apply_bqsr_bam
    run: subworkflows/base_quality_recalibration/base_quality_recalibration.cwl
    label: base_quality_recalibration
    'sbg:x': 1896.48779296875
    'sbg:y': 4264.0625
  - id: gatk_collect_alignment_summary_metrics_4_1_8_0
    in:
      - id: input
        source: base_quality_recalibration/gatk_apply_bqsr_bam
      - id: output_file_name
        source: gatk_collect_aln_summary_metrics_bqsr_output_file_name
      - id: reference
        source: reference_sequence
      - id: temporary_directory
        source: temporary_directory
    out:
      - id: gatk_collect_alignment_summary_metrics_txt
    run: >-
      command_line_tools/gatk_collect_alignment_summary_metrics_4.1.8.0/gatk_collect_alignment_summary_metrics_4.1.8.0.cwl
    label: GATK-CollectAlignmentSummaryMetrics
    'sbg:x': 3383.319580078125
    'sbg:y': 4114.6875
requirements:
  - class: SubworkflowFeatureRequirement
$schemas:
  - 'http://schema.org/version/latest/schemaorg-current-http.rdf'
's:author':
  - class: 's:Person'
    's:email': 'mailto:shahr2@mskcc.org'
    's:identifier': ''
    's:name': Ronak Shah
's:citation': ''
's:codeRepository': 'https://github.com/msk-access/nucleo'
's:contributor':
  - class: 's:Person'
    's:email': 'mailto:shahr2@mskcc.org'
    's:identifier': 'https://orcid.org/0000-0001-9042-6213'
    's:name': Ronak Shah
's:dateCreated': '2020-11-23'
's:license': 'https://spdx.org/licenses/Apache-2.0'
