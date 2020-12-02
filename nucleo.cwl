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
    'sbg:y': 1280.8125
  - id: gatk_base_recalibrator_known_sites
    type:
      type: array
      items: File
      inputBinding:
        prefix: '--known-sites'
    secondaryFiles:
      - .idx
    'sbg:x': 0
    'sbg:y': 2775.09375
  - id: sequencing-center
    type: string?
    'sbg:x': 0
    'sbg:y': 960.609375
  - id: run-date
    type: string?
    'sbg:x': 0
    'sbg:y': 1174.078125
  - id: sample
    type: string
    'sbg:x': 0
    'sbg:y': 1067.34375
  - id: read-structures
    type: 'string[]?'
    'sbg:x': 0
    'sbg:y': 1387.546875
  - id: read-group-id
    type: string
    'sbg:x': 0
    'sbg:y': 1494.28125
  - id: platform-unit
    type: string
    'sbg:x': 0
    'sbg:y': 1601.015625
  - id: platform-model
    type: string?
    'sbg:x': 0
    'sbg:y': 1707.75
  - id: platform
    type: string?
    'sbg:x': 0
    'sbg:y': 1814.484375
  - id: library
    type: string
    'sbg:x': 0
    'sbg:y': 2241.421875
  - id: validation_stringency
    type: string?
    'sbg:x': 0
    'sbg:y': 0
  - id: UBG_picard_SamToFastq_R1_output_fastq
    type: string
    'sbg:x': 0
    'sbg:y': 213.46875
  - id: UBG_picard_SamToFastq_R2_output_fastq
    type: string
    'sbg:x': 0
    'sbg:y': 106.734375
  - id: fastp_read2_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 4909.78125
  - id: fastp_read2_adapter_sequence
    type: string?
    'sbg:x': 0
    'sbg:y': 5016.515625
  - id: fastp_read1_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 5123.25
  - id: fastp_read1_adapter_sequence
    type: string?
    'sbg:x': 0
    'sbg:y': 5229.984375
  - id: fastp_minimum_read_length
    type: int?
    'sbg:x': 0
    'sbg:y': 5336.71875
  - id: fastp_html_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 5550.1875
  - id: fastp_json_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 5443.453125
  - id: sort_order
    type: string?
    'sbg:x': 0
    'sbg:y': 853.875
  - id: bwa_mem_T
    type: int?
    'sbg:x': 0
    'sbg:y': 5870.390625
  - id: bwa_mem_Y
    type: boolean?
    'sbg:x': 0
    'sbg:y': 5763.65625
  - id: UBG_picard_addRG_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 426.9375
  - id: create_bam_index
    type: boolean?
    'sbg:x': 0
    'sbg:y': 5656.921875
  - id: UBG_bwa_mem_output
    type: string
    'sbg:x': 0
    'sbg:y': 640.40625
  - id: bwa_mem_K
    type: int?
    'sbg:x': 0
    'sbg:y': 5977.125
  - id: UBG_gatk_merge_bam_alignment_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 533.671875
  - id: optical_duplicate_pixel_distance
    type: int?
    'sbg:x': 0
    'sbg:y': 2027.953125
  - id: gatk_mark_duplicates_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 2454.890625
  - id: gatk_mark_duplicates_duplication_metrics_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 2561.625
  - id: bedtools_merge_distance_between_features
    type: int?
    'sbg:x': 0
    'sbg:y': 6083.859375
  - id: bedtools_genomecov_option_bedgraph
    type: boolean?
    'sbg:x': 0
    'sbg:y': 6190.59375
  - id: base_recalibrator_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 7044.46875
  - id: apply_bqsr_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 7151.203125
  - id: abra2_window_size
    type: string?
    'sbg:x': 0
    'sbg:y': 7257.9375
  - id: abra2_soft_clip_contig
    type: string?
    'sbg:x': 0
    'sbg:y': 7364.671875
  - id: abra2_scoring_gap_alignments
    type: string?
    'sbg:x': 0
    'sbg:y': 7471.40625
  - id: UBG_abra2_output_bams
    type:
      - string
      - type: array
        items: string
    'sbg:x': 0
    'sbg:y': 747.140625
  - id: abra2_no_edge_complex_indel
    type: boolean?
    'sbg:x': 0
    'sbg:y': 7578.140625
  - id: abra2_maximum_mixmatch_rate
    type: float?
    'sbg:x': 0
    'sbg:y': 7684.875
  - id: abra2_maximum_average_depth
    type: int?
    'sbg:x': 0
    'sbg:y': 7791.609375
  - id: abra2_consensus_sequence
    type: boolean?
    'sbg:x': 0
    'sbg:y': 8005.078125
  - id: abra2_bam_index
    type: boolean?
    'sbg:x': 0
    'sbg:y': 8111.8125
  - id: fgbio_fastq_to_bam_input
    type:
      type: array
      items:
        items: File
        type: array
    'sbg:x': 0
    'sbg:y': 4269.375
  - id: UBG_picard_fixmateinformation_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 320.203125
  - id: merge_sam_files_sort_order
    type: string?
    'sbg:x': 0
    'sbg:y': 2134.6875
  - id: gatk_merge_sam_files_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 2348.15625
  - id: fgbio_collect_duplex_seq_metrics_intervals
    type: File?
    'sbg:x': 0
    'sbg:y': 4482.84375
  - id: fgbio_group_reads_by_umi_strategy
    type: string?
    'sbg:x': 0
    'sbg:y': 2988.5625
  - id: fgbio_group_reads_by_umi_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 3095.296875
  - id: fgbio_group_reads_by_umi_family_size_histogram
    type: string
    'sbg:x': 0
    'sbg:y': 3202.03125
  - id: fgbio_collect_duplex_seq_metrics_output_prefix
    type: string?
    'sbg:x': 0
    'sbg:y': 4376.109375
  - id: fgbio_collect_duplex_seq_metrics_duplex_umi_counts
    type: boolean?
    'sbg:x': 0
    'sbg:y': 4589.578125
  - id: fgbio_call_duplex_consensus_reads_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 4696.3125
  - id: fgbio_call_duplex_consensus_reads_min_reads
    type: 'int[]?'
    'sbg:x': 0
    'sbg:y': 4803.046875
  - id: BC_gatk_sam_to_fastq_output_name_R2
    type: string
    'sbg:x': 0
    'sbg:y': 6510.796875
  - id: BC_gatk_sam_to_fastq_output_name_R1
    type: string
    'sbg:x': 0
    'sbg:y': 6617.53125
  - id: BC_picard_fixmate_information_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 6297.328125
  - id: BC_picard_addRG_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 6404.0625
  - id: BC_gatk_merge_bam_alignment_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 6724.265625
  - id: fgbio_postprocessing_output_file_name_simplex
    type: string
    'sbg:x': 0
    'sbg:y': 2881.828125
  - id: gatk_collect_alignment_summary_metrics_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 2668.359375
  - id: abra2_contig_anchor
    type: string?
    'sbg:x': 0
    'sbg:y': 7898.34375
  - id: BC_abra2_output_bams
    type:
      - string
      - type: array
        items: string
    'sbg:x': 0
    'sbg:y': 6937.734375
  - id: fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex
    type: boolean?
    'sbg:x': 0
    'sbg:y': 3308.765625
  - id: fgbio_filter_consensus_read_output_file_name_simplex_duplex
    type: string
    'sbg:x': 0
    'sbg:y': 3415.5
  - id: fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics
    type: string
    'sbg:x': 0
    'sbg:y': 3522.234375
  - id: fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics
    type: string?
    'sbg:x': 0
    'sbg:y': 3628.96875
  - id: fgbio_filter_consensus_read_output_file_name_duplex
    type: string?
    'sbg:x': 0
    'sbg:y': 3735.703125
  - id: fgbio_filter_consensus_read_min_reads_duplex
    type: 'int[]?'
    'sbg:x': 0
    'sbg:y': 3949.171875
  - id: fgbio_filter_consensus_read_min_base_quality_simplex_duplex
    type: int?
    'sbg:x': 0
    'sbg:y': 4055.90625
  - id: fgbio_filter_consensus_read_min_base_quality_duplex
    type: int?
    'sbg:x': 0
    'sbg:y': 4162.640625
  - id: BC_bwa_mem_output
    type: string
    'sbg:x': 0
    'sbg:y': 6831
  - id: fgbio_filter_consensus_read_min_reads_simplex_duplex
    type: 'int[]?'
    'sbg:x': 0
    'sbg:y': 3842.4375
  - id: picard_addRG_sort_order
    type: string?
    'sbg:x': 0
    'sbg:y': 1921.21875
outputs:
  - id: fastp_html_output
    outputSource:
      - uncollapsed_bam_generation/fastp_html_output
    type: File
    'sbg:x': 1916.0830078125
    'sbg:y': 3912.375
  - id: fastp_json_output
    outputSource:
      - uncollapsed_bam_generation/fastp_json_output
    type: File
    'sbg:x': 1916.0830078125
    'sbg:y': 3805.640625
  - id: gatk_collect_alignment_summary_metrics_txt
    outputSource:
      - uncollapsed_bam_generation/gatk_collect_alignment_summary_metrics_txt
    type: File
    'sbg:x': 1916.0830078125
    'sbg:y': 3698.90625
  - id: indel_realignment_bam
    outputSource:
      - uncollapsed_bam_generation/indel_realignment_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1916.0830078125
    'sbg:y': 3592.171875
  - id: picard_mark_duplicates_metrics
    outputSource:
      - uncollapsed_bam_generation/picard_mark_duplicates_metrics
    type: File
    'sbg:x': 1916.0830078125
    'sbg:y': 3485.4375
  - id: uncollapsed_bam
    outputSource:
      - uncollapsed_bam_generation/uncollapsed_bam
    type: File?
    secondaryFiles:
      - ^.bai
    'sbg:x': 1916.0830078125
    'sbg:y': 3378.703125
  - id: gatk_collect_alignment_summary_metrics_txt_simplex
    outputSource:
      - bam_collapsing/gatk_collect_alignment_summary_metrics_txt_simplex
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 3362.1328125
  - id: gatk_collect_alignment_summary_metrics_txt_duplex
    outputSource:
      - bam_collapsing/gatk_collect_alignment_summary_metrics_txt_duplex
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 3468.8671875
  - id: gatk_collect_alignment_summary_metrics_txt_collapsed
    outputSource:
      - bam_collapsing/gatk_collect_alignment_summary_metrics_txt_collapsed
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 3575.6015625
  - id: fgbio_postprocessing_simplex_bam
    outputSource:
      - bam_collapsing/fgbio_postprocessing_simplex_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3349.681884765625
    'sbg:y': 3682.3359375
  - id: fgbio_group_reads_by_umi_histogram
    outputSource:
      - bam_collapsing/fgbio_group_reads_by_umi_histogram
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 3789.0703125
  - id: fgbio_group_reads_by_umi_bam
    outputSource:
      - bam_collapsing/fgbio_group_reads_by_umi_bam
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 3895.8046875
  - id: fgbio_filter_consensus_reads_duplex_bam
    outputSource:
      - bam_collapsing/fgbio_filter_consensus_reads_duplex_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3349.681884765625
    'sbg:y': 4002.5390625
  - id: fgbio_collect_duplex_seq_metrics_umi_counts
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_umi_counts
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 4109.2734375
  - id: fgbio_collect_duplex_seq_metrics_family_size
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_family_size
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 4216.0078125
  - id: fgbio_collect_duplex_seq_metrics_duplex_yield_metrics
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_yield_metrics
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 4322.7421875
  - id: fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 4429.4765625
  - id: fgbio_collect_duplex_seq_metrics_duplex_qc
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_qc
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 4536.2109375
  - id: fgbio_collect_duplex_seq_metrics_duplex_family_size
    outputSource:
      - bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_family_size
    type: File
    'sbg:x': 3349.681884765625
    'sbg:y': 4642.9453125
  - id: fgbio_collapsed_bam
    outputSource:
      - bam_collapsing/fgbio_collapsed_bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 3349.681884765625
    'sbg:y': 4749.6796875
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
        default: 30
        source: bwa_mem_T
      - id: sort_order
        default: coordinate
        source: sort_order
      - id: picard_addRG_output_file_name
        source: UBG_picard_addRG_output_file_name
      - id: bwa_mem_output
        source: UBG_bwa_mem_output
      - id: bwa_mem_K
        default: 1000000
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
        default: true
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
      - id: abra2_no_edge_complex_indel
        default: true
        source: abra2_no_edge_complex_indel
      - id: gatk_base_recalibrator_known_sites
        source:
          - gatk_base_recalibrator_known_sites
      - id: base_recalibrator_output_file_name
        source: base_recalibrator_output_file_name
      - id: apply_bqsr_output_file_name
        source: apply_bqsr_output_file_name
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
    out:
      - id: gatk_sam_to_fastq_unpaired_fastq
      - id: fastp_unpaired2_output
      - id: fastp_unpaired1_output
      - id: fastp_json_output
      - id: fastp_html_output
      - id: picard_mark_duplicates_metrics
      - id: indel_realignment_bam
      - id: uncollapsed_bam
      - id: gatk_collect_alignment_summary_metrics_txt
    run: uncollapsed_bam_generation/uncollapsed_bam_generation.cwl
    label: Uncollapsed BAM Generation
    'sbg:x': 611.53125
    'sbg:y': 3705.90625
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
        default: 30
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
        default: 1000000
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
    'sbg:x': 1916.0830078125
    'sbg:y': 4376.109375
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
