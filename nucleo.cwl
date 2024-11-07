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
    'sbg:y': 1495.8125
  - id: gatk_base_recalibrator_known_sites
    type:
      type: array
      items: File
      inputBinding:
        prefix: '--known-sites'
    secondaryFiles:
      - .idx
    'sbg:x': 0
    'sbg:y': 3098.46875
  - id: sequencing-center
    type: string?
    'sbg:x': 0
    'sbg:y': 1175.28125
  - id: run-date
    type: string?
    'sbg:x': 0
    'sbg:y': 1388.96875
  - id: sample
    type: string
    'sbg:x': 0
    'sbg:y': 1282.125
  - id: read-structures
    type: 'string[]?'
    'sbg:x': 0
    'sbg:y': 1602.65625
  - id: read-group-id
    type: string
    'sbg:x': 0
    'sbg:y': 1709.5
  - id: platform-unit
    type: string
    'sbg:x': 0
    'sbg:y': 1816.34375
  - id: platform-model
    type: string?
    'sbg:x': 0
    'sbg:y': 1923.1875
  - id: platform
    type: string?
    'sbg:x': 0
    'sbg:y': 2030.03125
  - id: library
    type: string
    'sbg:x': 0
    'sbg:y': 2457.40625
  - id: validation_stringency
    type: string?
    'sbg:x': 0
    'sbg:y': 0
  - id: UBG_picard_SamToFastq_R1_output_fastq
    type: string
    'sbg:x': 0
    'sbg:y': 213.6875
  - id: UBG_picard_SamToFastq_R2_output_fastq
    type: string
    'sbg:x': 0
    'sbg:y': 106.84375
  - id: fastp_read2_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 5449.03125
  - id: fastp_read2_adapter_sequence
    type: string?
    'sbg:x': 0
    'sbg:y': 5555.875
  - id: fastp_read1_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 5662.71875
  - id: fastp_read1_adapter_sequence
    type: string?
    'sbg:x': 0
    'sbg:y': 5769.5625
  - id: fastp_minimum_read_length
    type: int?
    'sbg:x': 0
    'sbg:y': 5876.40625
  - id: fastp_html_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 6303.78125
  - id: fastp_json_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 6196.9375
  - id: sort_order
    type: string?
    'sbg:x': 0
    'sbg:y': 1068.4375
  - id: bwa_mem_T
    type: int?
    'sbg:x': 0
    'sbg:y': 6838
  - id: bwa_mem_Y
    type: boolean?
    'sbg:x': 0
    'sbg:y': 6731.15625
  - id: UBG_picard_addRG_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 427.375
  - id: create_bam_index
    type: boolean?
    'sbg:x': 0
    'sbg:y': 6624.3125
  - id: UBG_bwa_mem_output
    type: string
    'sbg:x': 0
    'sbg:y': 641.0625
  - id: bwa_mem_K
    type: int?
    'sbg:x': 0
    'sbg:y': 6944.84375
  - id: UBG_gatk_merge_bam_alignment_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 534.21875
  - id: optical_duplicate_pixel_distance
    type: int?
    'sbg:x': 0
    'sbg:y': 2243.71875
  - id: gatk_mark_duplicates_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 2671.09375
  - id: gatk_mark_duplicates_duplication_metrics_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 2777.9375
  - id: bedtools_merge_distance_between_features
    type: int?
    'sbg:x': 0
    'sbg:y': 7051.6875
  - id: bedtools_genomecov_option_bedgraph
    type: boolean?
    'sbg:x': 0
    'sbg:y': 7158.53125
  - id: apply_bqsr_output_file_name
    type: string
    'sbg:x': 611.515625
    'sbg:y': 4594.28125
  - id: abra2_window_size
    type: string?
    'sbg:x': 0
    'sbg:y': 8120.125
  - id: abra2_soft_clip_contig
    type: string?
    'sbg:x': 0
    'sbg:y': 8226.96875
  - id: abra2_scoring_gap_alignments
    type: string?
    'sbg:x': 0
    'sbg:y': 8333.8125
  - id: UBG_abra2_output_bams
    type:
      - string
      - type: array
        items: string
    'sbg:x': 0
    'sbg:y': 854.75
  - id: abra2_no_edge_complex_indel
    type: boolean?
    'sbg:x': 0
    'sbg:y': 8547.5
  - id: abra2_maximum_mixmatch_rate
    type: float?
    'sbg:x': 0
    'sbg:y': 8654.34375
  - id: abra2_maximum_average_depth
    type: int?
    'sbg:x': 0
    'sbg:y': 8761.1875
  - id: abra2_consensus_sequence
    type: boolean?
    'sbg:x': 0
    'sbg:y': 8974.875
  - id: abra2_bam_index
    type: boolean?
    'sbg:x': 0
    'sbg:y': 9081.71875
  - id: fgbio_fastq_to_bam_input
    type:
      type: array
      items:
        items: File
        type: array
    'sbg:x': 0
    'sbg:y': 4701.125
  - id: UBG_picard_fixmateinformation_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 320.53125
  - id: merge_sam_files_sort_order
    type: string?
    'sbg:x': 0
    'sbg:y': 2350.5625
  - id: gatk_merge_sam_files_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 2564.25
  - id: fgbio_collect_duplex_seq_metrics_intervals
    type: File?
    'sbg:x': 0
    'sbg:y': 4914.8125
  - id: fgbio_group_reads_by_umi_strategy
    type: string?
    'sbg:x': 0
    'sbg:y': 3419
  - id: fgbio_group_reads_by_umi_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 3525.84375
  - id: fgbio_group_reads_by_umi_family_size_histogram
    type: string
    'sbg:x': 0
    'sbg:y': 3632.6875
  - id: fgbio_collect_duplex_seq_metrics_output_prefix
    type: string?
    'sbg:x': 0
    'sbg:y': 4807.96875
  - id: fgbio_collect_duplex_seq_metrics_duplex_umi_counts
    type: boolean?
    'sbg:x': 0
    'sbg:y': 5021.65625
  - id: fgbio_call_duplex_consensus_reads_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 5128.5
  - id: fgbio_call_duplex_consensus_reads_min_reads
    type: 'int[]?'
    'sbg:x': 0
    'sbg:y': 5235.34375
  - id: BC_gatk_sam_to_fastq_output_name_R2
    type: string
    'sbg:x': 0
    'sbg:y': 7479.0625
  - id: BC_gatk_sam_to_fastq_output_name_R1
    type: string
    'sbg:x': 0
    'sbg:y': 7585.90625
  - id: BC_picard_fixmate_information_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 7265.375
  - id: BC_picard_addRG_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 7372.21875
  - id: BC_gatk_merge_bam_alignment_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 7692.75
  - id: fgbio_postprocessing_output_file_name_simplex
    type: string
    'sbg:x': 0
    'sbg:y': 3312.15625
  - id: gatk_collect_alignment_summary_metrics_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 2991.625
  - id: abra2_contig_anchor
    type: string?
    'sbg:x': 0
    'sbg:y': 8868.03125
  - id: BC_abra2_output_bams
    type:
      - string
      - type: array
        items: string
    'sbg:x': 0
    'sbg:y': 7906.4375
  - id: fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex
    type: boolean?
    'sbg:x': 0
    'sbg:y': 3739.53125
  - id: fgbio_filter_consensus_read_output_file_name_simplex_duplex
    type: string
    'sbg:x': 0
    'sbg:y': 3846.375
  - id: fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics
    type: string
    'sbg:x': 0
    'sbg:y': 3953.21875
  - id: fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics
    type: string?
    'sbg:x': 0
    'sbg:y': 4060.0625
  - id: fgbio_filter_consensus_read_output_file_name_duplex
    type: string?
    'sbg:x': 0
    'sbg:y': 4166.90625
  - id: fgbio_filter_consensus_read_min_reads_duplex
    type: 'int[]?'
    'sbg:x': 0
    'sbg:y': 4380.59375
  - id: fgbio_filter_consensus_read_min_base_quality_simplex_duplex
    type: int?
    'sbg:x': 0
    'sbg:y': 4487.4375
  - id: fgbio_filter_consensus_read_min_base_quality_duplex
    type: int?
    'sbg:x': 0
    'sbg:y': 4594.28125
  - id: BC_bwa_mem_output
    type: string
    'sbg:x': 0
    'sbg:y': 7799.59375
  - id: fgbio_filter_consensus_read_min_reads_simplex_duplex
    type: 'int[]?'
    'sbg:x': 0
    'sbg:y': 4273.75
  - id: picard_addRG_sort_order
    type: string?
    'sbg:x': 0
    'sbg:y': 2136.875
  - id: disable_trim_poly_g
    type: boolean?
    'sbg:x': 0
    'sbg:y': 6410.625
  - id: disable_quality_filtering
    type: boolean?
    'sbg:x': 0
    'sbg:y': 6517.46875
  - id: gatk_base_recalibrator_add_output_sam_program_record
    type: boolean?
    'sbg:x': 0
    'sbg:y': 3205.3125
  - id: gatk_collect_aln_summary_metrics_bqsr_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 2884.78125
  - id: abra2_no_sort
    type: boolean?
    'sbg:x': 0
    'sbg:y': 8440.65625
  - id: base_recalibrator_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 8013.28125
  - id: temporary_directory
    type: string?
    'sbg:x': 0
    'sbg:y': 961.59375
  - id: fgbio_async_io
    type: string?
    'sbg:x': 0
    'sbg:y': 5342.1875
  - id: fastp_max_len_read2
    type: int?
    'sbg:x': 0
    'sbg:y': 5983.25
  - id: fastp_max_len_read1
    type: int?
    'sbg:x': 0
    'sbg:y': 6090.09375
  - id: UBG_abra2_targets
    type: File
    'sbg:x': 0
    'sbg:y': 747.90625
outputs:
  - id: fastp_html_output
    outputSource:
      - uncollapsed_bam_generation/fastp_html_output
    type: File
    'sbg:x': 1884.178466796875
    'sbg:y': 4160.28125
  - id: fastp_json_output
    outputSource:
      - uncollapsed_bam_generation/fastp_json_output
    type: File
    'sbg:x': 1884.178466796875
    'sbg:y': 4053.4375
  - id: gatk_collect_alignment_summary_metrics_txt_uncollapsed
    outputSource:
      - >-
        gatk_collect_alignment_summary_metrics_4_1_8_0/gatk_collect_alignment_summary_metrics_txt
    type: File
    'sbg:x': 3887.41650390625
    'sbg:y': 4540.859375
  - id: indel_realignment_bam
    outputSource:
      - uncollapsed_bam_generation/indel_realignment_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1884.178466796875
    'sbg:y': 3946.59375
  - id: picard_mark_duplicates_metrics
    outputSource:
      - uncollapsed_bam_generation/picard_mark_duplicates_metrics
    type: File
    'sbg:x': 1884.178466796875
    'sbg:y': 3839.75
  - id: gatk_collect_alignment_summary_metrics_txt_simplex
    outputSource:
      - bam_collapsing/gatk_collect_alignment_summary_metrics_txt_simplex
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 3825.375
  - id: gatk_collect_alignment_summary_metrics_txt_duplex
    outputSource:
      - bam_collapsing/gatk_collect_alignment_summary_metrics_txt_duplex
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 3932.21875
  - id: gatk_collect_alignment_summary_metrics_txt_collapsed
    outputSource:
      - bam_collapsing/gatk_collect_alignment_summary_metrics_txt_collapsed
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 4039.0625
  - id: fgbio_postprocessing_simplex_bam
    outputSource:
      - bam_collapsing/fgbio_postprocessing_simplex_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3371.01025390625
    'sbg:y': 4294.75
  - id: fgbio_group_reads_by_umi_histogram
    outputSource:
      - bam_collapsing/fgbio_group_reads_by_umi_histogram
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 4401.59375
  - id: fgbio_group_reads_by_umi_bam
    outputSource:
      - bam_collapsing/fgbio_group_reads_by_umi_bam
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 4508.4375
  - id: fgbio_filter_consensus_reads_duplex_bam
    outputSource:
      - bam_collapsing/fgbio_filter_consensus_reads_duplex_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3371.01025390625
    'sbg:y': 4615.28125
  - id: fgbio_collect_duplex_seq_metrics_umi_counts
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_umi_counts
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 4722.125
  - id: fgbio_collect_duplex_seq_metrics_family_size
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_family_size
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 4828.96875
  - id: fgbio_collect_duplex_seq_metrics_duplex_yield_metrics
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_yield_metrics
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 4935.8125
  - id: fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 5042.65625
  - id: fgbio_collect_duplex_seq_metrics_duplex_qc
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_qc
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 5149.5
  - id: fgbio_collect_duplex_seq_metrics_duplex_family_size
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_family_size
    type: File
    'sbg:x': 3371.01025390625
    'sbg:y': 5256.34375
  - id: fgbio_collapsed_bam
    outputSource:
      - bam_collapsing/fgbio_collapsed_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3371.01025390625
    'sbg:y': 5363.1875
  - id: uncollapsed_bam
    outputSource:
      - base_quality_recalibration/gatk_apply_bqsr_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3371.01025390625
    'sbg:y': 3718.53125
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
      - id: abra2_maximum_mixmatch_rate
        default: 0.1
        source: abra2_maximum_mixmatch_rate
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
        default: 95
        source: fastp_max_len_read1
      - id: fastp_max_len_read2
        default: 95
        source: fastp_max_len_read2
      - id: abra2_targets
        source: UBG_abra2_targets
        default: '/work/access/production/resources/msk-access/v2.0/regions_of_interest/versions/v1.0/MSK-ACCESS-v2_targetsAllwFP.bed'
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
    'sbg:x': 611.515625
    'sbg:y': 4116.4375
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
    'sbg:x': 1884.178466796875
    'sbg:y': 4856.96875
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
    'sbg:x': 1884.178466796875
    'sbg:y': 4316.125
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
    'sbg:x': 3371.01025390625
    'sbg:y': 4166.90625
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
