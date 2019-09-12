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
    'sbg:y': 962.15625
  - id: fastq2
    type: File
    'sbg:x': 0
    'sbg:y': 855.25
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
    'sbg:x': 373.8125
    'sbg:y': 406.625
  - id: known_sites_1
    type: File
    secondaryFiles:
      - .idx
    'sbg:x': 373.8125
    'sbg:y': 1196.96875
  - id: bed_file
    type: File
    'sbg:x': 0
    'sbg:y': 1282.875
  - id: read_group_sample_name
    type: string
    'sbg:x': 373.8125
    'sbg:y': 513.53125
  - id: read_group_platform_unit
    type: string
    'sbg:x': 373.8125
    'sbg:y': 620.4375
  - id: read_group_library
    type: int
    'sbg:x': 373.8125
    'sbg:y': 727.34375
  - id: read_group_identifier
    type: string
    'sbg:x': 373.8125
    'sbg:y': 834.25
  - id: sort_first_pass_output_file_name
    type: string
    'sbg:x': 0
    'sbg:y': 534.53125
  - id: output_name_collapsed_gzip_R2
    type: string?
    'sbg:x': 0
    'sbg:y': 641.4375
  - id: output_name_collapsed_gzip_R1
    type: string?
    'sbg:x': 0
    'sbg:y': 748.34375
  - id: collapsing_aln_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 1175.96875
  - id: collapsing_picard_output_file_name
    type: string?
    'sbg:x': 0
    'sbg:y': 1069.0625
  - id: standard_aln_output_file_name
    type: string?
    label: standard_aln_output_file_name
    'sbg:x': 373.8125
    'sbg:y': 299.71875
  - id: standard_picard_addrg_output_filename
    type: string?
    label: standard_picard_addrg_output_filename
    'sbg:x': 0
    'sbg:y': 427.625
  - id: known_sites_2
    type: File?
    secondaryFiles:
      - .idx
    'sbg:x': 373.8125
    'sbg:y': 1090.0625
outputs:
  - id: composite_umi_frequencies
    outputSource:
      - marianas_process_loop_umi_cwl/composite_umi_frequencies
    type: File
    'sbg:x': 680.87744140625
    'sbg:y': 855.25
  - id: clipping_info
    outputSource:
      - marianas_process_loop_umi_cwl/clipping_info
    type: File
    'sbg:x': 680.87744140625
    'sbg:y': 962.15625
  - id: md_bam
    outputSource:
      - standard_bam_processing_cwl/md_bam
    type: File
    label: mark_duplicates_bam
    secondaryFiles:
      - ^.bai
    'sbg:x': 1139.8087158203125
    'sbg:y': 515.4375
  - id: clstats2
    outputSource:
      - standard_bam_processing_cwl/clstats2
    type: File
    label: trimming_stats_read2
    'sbg:x': 1139.8087158203125
    'sbg:y': 622.34375
  - id: clstats1
    outputSource:
      - standard_bam_processing_cwl/clstats1
    type: File
    label: trimming_stats_read1
    'sbg:x': 1139.8087158203125
    'sbg:y': 729.25
  - id: bqsr_bam
    outputSource:
      - standard_bam_processing_cwl/bqsr_bam
    type: File?
    label: standard_processed_bam
    secondaryFiles:
      - ^.bai
    'sbg:x': 1139.8087158203125
    'sbg:y': 836.15625
  - id: unfiltered-bam
    outputSource:
      - bam_collapsing/unfiltered-bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1788.5467529296875
    'sbg:y': 0
  - id: simplex-bam
    outputSource:
      - bam_collapsing/simplex-bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1788.5467529296875
    'sbg:y': 106.90625
  - id: second_pass_insertions
    outputSource:
      - bam_collapsing/second_pass_insertions
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 213.8125
  - id: second_pass_alt_alleles
    outputSource:
      - bam_collapsing/second_pass_alt_alleles
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 320.71875
  - id: pileup_without_duplicates
    outputSource:
      - bam_collapsing/pileup_without_duplicates
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 427.625
  - id: intervals_without_duplicates
    outputSource:
      - bam_collapsing/intervals_without_duplicates
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 534.53125
  - id: intervals
    outputSource:
      - bam_collapsing/intervals
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 641.4375
  - id: gzip_read1
    outputSource:
      - bam_collapsing/gzip_read1
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 855.25
  - id: gzip_read2
    outputSource:
      - bam_collapsing/gzip_read2
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 748.34375
  - id: first_pass_insertions
    outputSource:
      - bam_collapsing/first_pass_insertions
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 962.15625
  - id: duplex-bam
    outputSource:
      - bam_collapsing/duplex-bam
    type: File
    secondaryFiles:
      - ^.bai
    'sbg:x': 1788.5467529296875
    'sbg:y': 1069.0625
  - id: collapsed_fastq_2
    outputSource:
      - bam_collapsing/collapsed_fastq_2
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 1175.96875
  - id: alt_allele_file
    outputSource:
      - bam_collapsing/alt_allele_file
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 1389.78125
  - id: alignment_metrics_unfiltered
    outputSource:
      - bam_collapsing/alignment_metrics_unfiltered
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 1496.6875
  - id: alignment_metrics_simplex
    outputSource:
      - bam_collapsing/alignment_metrics_simplex
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 1603.59375
  - id: alignment_metrics_duplex
    outputSource:
      - bam_collapsing/alignment_metrics_duplex
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 1710.5
  - id: collapsed_fastq_1
    outputSource:
      - bam_collapsing/collapsed_fastq_1
    type: File
    'sbg:x': 1788.5467529296875
    'sbg:y': 1282.875
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
    'sbg:x': 373.8125
    'sbg:y': 962.15625
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
        source: standard_aln_output_file_name
      - id: output_file_name
        source: standard_picard_addrg_output_filename
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
      - id: ignore_bad_assembly
        default: true
      - id: contig_anchor
        default: '10,1'
      - id: consensus_sequence
        default: true
      - id: suppress_warn
        default: true
      - id: stringency
        default: 3
      - id: quality
        default: 1
      - id: length
        default: 25
      - id: adapter2
        default: AGATCGGAAGAGC
      - id: adapter
        default: GATCGGAAGAGC
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
    'sbg:x': 680.87744140625
    'sbg:y': 664.34375
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
      - id: output_file
    run: bam_collapsing/bam_collapsing.cwl
    label: bam_collapsing
    'sbg:x': 1139.8087158203125
    'sbg:y': 1069.0625
requirements:
  - class: SubworkflowFeatureRequirement
