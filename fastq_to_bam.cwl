class: Workflow
cwlVersion: v1.0
id: fastq_to_bam
doc: >-
  This workflow takes a READ1 and READ2 fastq.gz file generated for MSK-ACCESS
  assay and generated four different Binary Alignment Map file along with
  alignment metrics for each.
label: fastq_to_bam.cwl
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: fastq1
    type: File
    doc: Gzipped Fastq File for READ1
    'sbg:x': 0
    'sbg:y': 2208.2734375
  - id: fastq2
    type: File
    doc: Gzipped Fastq File for READ2
    'sbg:x': 0
    'sbg:y': 2101.6171875
  - id: reference
    type: File
    doc: >-
      The reference sequence in a single reference sequence in FASTA format,
      with all contigs in the same file, validated according to the FASTA
      standard. It has multiple secondary file associated with it ending in
      ".dict, .fai, .amb, .ann, .bwt, .pac, .index, .rbwt, .rpac, .rsa, .sa"
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
    'sbg:x': 373.78125
    'sbg:y': 426.78125
  - id: known_sites_1
    type: File
    doc: >-
      A database of known polymorphic sites on VCF format, Ex: DBSNP or
      Mills_and_1000G. Note: ".vcf.idx" secaondary file should be present where
      the ".vcf" file is located
    secondaryFiles:
      - .idx
    'sbg:x': 373.78125
    'sbg:y': 2282.71875
  - id: bed_file
    type: File
    doc: >-
      Targets in BED file format used by Waltz to generate the PileUp for
      collapsing of the BAM file.The genotype from positions in this bed file
      will be used as the consensus base if min_consensus_percent threshold is
      not reached. Otherwise, the reference base from the supplied
      reference_fasta will be used
    'sbg:x': 0
    'sbg:y': 2528.3984375
  - id: read_group_sample_name
    type: string
    doc: >-
      SM = Sample

      The name of the sample sequenced in this read group. GATK tools treat all
      read groups with the same SM value as containing sequencing data for the
      same sample, and this is also the name that will be used for the sample
      column in the VCF file. Therefore it's critical that the SM field be
      specified correctly. When sequencing pools of samples, use a pool name
      instead of an individual sample name. Note, when we say pools, we mean
      samples that are not individually barcoded. In the case of multiplexing
      (often confused with pooling) where you know which reads come from each
      sample and you have simply run the samples together in one lane, you can
      keep the SM tag as the sample name and not the "pooled name".
    'sbg:x': 373.78125
    'sbg:y': 746.9453125
  - id: read_group_platform_unit
    type: string
    doc: >-
      PU = Platform Unit

      The PU holds three types of information, the
      {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}. The {FLOWCELL_BARCODE} refers
      to the unique identifier for a particular flow cell. The {LANE} indicates
      the lane of the flow cell and the {SAMPLE_BARCODE} is a
      sample/library-specific identifier. Although the PU is not required by
      GATK but takes precedence over ID for base recalibration if it is present.
      In the example shown earlier, two read group fields, ID and PU,
      appropriately differentiate flow cell lane, marked by .2, a factor that
      contributes to batch effects.
    'sbg:x': 373.78125
    'sbg:y': 853.6796875
  - id: read_group_library
    type: int
    doc: >-
      LB = DNA preparation library identifier

      MarkDuplicates uses the LB field to determine which read groups might
      contain molecular duplicates, in case the same DNA library was sequenced
      on multiple lanes.
    'sbg:x': 373.78125
    'sbg:y': 960.4140625
  - id: read_group_identifier
    type: string
    doc: >-
      ID = Read group identifier

      This tag identifies which read group each read belongs to, so each read
      group's ID must be unique. It is referenced both in the read group
      definition line in the file header (starting with @RG) and in the RG:Z tag
      for each read record.
    'sbg:x': 373.78125
    'sbg:y': 1067.1484375
  - id: sort_first_pass_output_file_name
    type: string
    doc: Name for the Marianas Duplex Collapsing First Pass output TXT file.
    'sbg:x': 0
    'sbg:y': 1247.9765625
  - id: output_name_collapsed_gzip_R2
    type: string?
    doc: Name of the output collapsed READ1 gzip fastq file.
    'sbg:x': 0
    'sbg:y': 1354.671875
  - id: output_name_collapsed_gzip_R1
    type: string?
    doc: Name of the output collapsed READ1 gzip fastq file.
    'sbg:x': 0
    'sbg:y': 1461.40625
  - id: collapsing_aln_output_file_name
    type: string?
    doc: Name of the SAM format output file created by bwa mem for collapsing step.
    'sbg:x': 0
    'sbg:y': 2421.703125
  - id: collapsing_picard_output_file_name
    type: string?
    doc: >-
      Name of the BAM format output file created by Picard
      AddOrReplaceReadGroups for collapsing step.
    'sbg:x': 0
    'sbg:y': 2314.96875
  - id: standard_aln_output_file_name
    type: string?
    label: standard_aln_output_file_name
    doc: >-
      Name of the SAM format output file created by bwa mem for standard bam
      processing step
    'sbg:x': 373.78125
    'sbg:y': 0
  - id: standard_picard_addrg_output_filename
    type: string?
    label: standard_picard_addrg_output_filename
    doc: >-
      Name of the BAM format output file created by Picard
      AddOrReplaceReadGroups for standard bam processing step.
    'sbg:x': 0
    'sbg:y': 1141.28125
  - id: known_sites_2
    type: File?
    doc: >-
      A database of known polymorphic sites on VCF format, Ex: DBSNP or
      Mills_and_1000G. Note: ".vcf.idx" secaondary file should be present where
      the ".vcf" file is located
    secondaryFiles:
      - .idx
    'sbg:x': 373.78125
    'sbg:y': 2176.0625
  - id: wobble
    type: int?
    'sbg:x': 0
    'sbg:y': 607.7265625
  - id: read_group_sequnecing_center
    type: string?
    doc: RGCN tag for BAM file indicating where the data is sequenced.
    'sbg:x': 373.78125
    'sbg:y': 533.4765625
  - id: read_group_sequencing_platform
    type: string?
    doc: BAM Tag describing the Platform used to generate the sequencing data.
    'sbg:x': 373.78125
    'sbg:y': 640.2109375
  - id: mismatches
    type: int?
    'sbg:x': 0
    'sbg:y': 1568.1015625
  - id: min_map_quality
    type: int?
    'sbg:x': 0
    'sbg:y': 1674.796875
  - id: min_consensus_percent
    type: int?
    'sbg:x': 0
    'sbg:y': 1781.4921875
  - id: min_base_quality
    type: int?
    'sbg:x': 0
    'sbg:y': 1888.1875
  - id: key
    type:
      - 'null'
      - type: array
        items: string
        inputBinding:
          prefix: '-k'
    'sbg:x': 0
    'sbg:y': 1994.921875
  - id: abra_collapsing_number_of_threads
    type: int?
    'sbg:x': 692.921875
    'sbg:y': 1728.0859375
  - id: adapter
    type: string?
    'sbg:x': 373.78125
    'sbg:y': 3136.125
  - id: adapter2
    type: string?
    'sbg:x': 373.78125
    'sbg:y': 3029.46875
  - id: assume_sorted
    type: boolean?
    'sbg:x': 373.78125
    'sbg:y': 2922.8125
  - id: bqsr_read_filter
    type: string?
    'sbg:x': 373.78125
    'sbg:y': 2816.15625
  - id: consensus_sequence
    type: boolean?
    'sbg:x': 373.78125
    'sbg:y': 2709.5
  - id: contig_anchor
    type: string?
    'sbg:x': 373.78125
    'sbg:y': 2602.8046875
  - id: option_bedgraph
    type: boolean?
    'sbg:x': 373.78125
    'sbg:y': 1387.2734375
  - id: number_of_threads
    type: int?
    'sbg:x': 373.78125
    'sbg:y': 1493.96875
  - id: maximum_mixmatch_rate
    type: float?
    'sbg:x': 373.78125
    'sbg:y': 1600.625
  - id: maximum_average_depth
    type: int?
    'sbg:x': 373.78125
    'sbg:y': 1707.3203125
  - id: M
    type: boolean?
    'sbg:x': 373.78125
    'sbg:y': 1962.671875
  - id: length
    type: int?
    'sbg:x': 373.78125
    'sbg:y': 2069.3671875
  - id: P
    type: boolean?
    'sbg:x': 373.78125
    'sbg:y': 1280.578125
  - id: ignore_bad_assembly
    type: boolean?
    'sbg:x': 373.78125
    'sbg:y': 2389.4140625
  - id: scoring_gap_alignments
    type: string?
    'sbg:x': 373.78125
    'sbg:y': 320.0859375
  - id: soft_clip_contig
    type: string?
    'sbg:x': 373.78125
    'sbg:y': 213.3515625
  - id: window_size
    type: string?
    'sbg:x': 0
    'sbg:y': 714.3828125
  - id: validation_stringency
    type: string?
    'sbg:x': 0
    'sbg:y': 821.078125
  - id: trim_galore_number_of_threads
    type: int?
    'sbg:x': 0
    'sbg:y': 927.8125
  - id: stringency
    type: int?
    'sbg:x': 0
    'sbg:y': 1034.546875
  - id: sort_order
    type: string?
    'sbg:x': 373.78125
    'sbg:y': 106.65625
  - id: quality
    type: int?
    'sbg:x': 373.78125
    'sbg:y': 1173.8828125
outputs:
  - id: composite_umi_frequencies
    outputSource:
      - marianas_process_loop_umi_cwl/composite_umi_frequencies
    type: File
    doc: >-
      This is text file consisting of frequencines of unique molecular
      identifier as seen by Marianas ProcessLoopUMIFastq
    'sbg:x': 692.921875
    'sbg:y': 1514.65625
  - id: clipping_info
    outputSource:
      - marianas_process_loop_umi_cwl/clipping_info
    type: File
    doc: >-
      File having information about all the clipped unique molecular identifiers
      from the fastq.gz files by Marianas ProcessLoopUMIFastq
    'sbg:x': 692.921875
    'sbg:y': 1621.3515625
  - id: md_bam
    outputSource:
      - standard_bam_processing_cwl/md_bam
    type: File
    label: mark_duplicates_bam
    doc: >-
      Binary Alignment Map (BAM) File generated after marking duplicate reads
      using Picard MarkDuplicate tool.
    secondaryFiles:
      - ^.bai
    'sbg:x': 1508.0369873046875
    'sbg:y': 1321.328125
  - id: clstats2
    outputSource:
      - standard_bam_processing_cwl/clstats2
    type: File
    label: trimming_stats_read2
    doc: Trimming statistics generated by TrimGalore/Cutadapt for READ2
    'sbg:x': 1508.0369873046875
    'sbg:y': 1428.0234375
  - id: clstats1
    outputSource:
      - standard_bam_processing_cwl/clstats1
    type: File
    label: trimming_stats_read1
    doc: Trimming statistics generated by TrimGalore/Cutadapt for READ1
    'sbg:x': 1508.0369873046875
    'sbg:y': 1534.7578125
  - id: bqsr_bam
    outputSource:
      - standard_bam_processing_cwl/bqsr_bam
    type: File?
    label: standard_processed_bam
    doc: >-
      Base Recalibrated Binary Alignment Map format file generated using GATK
      BaseRecalibrator and ApplyBQSR tool.
    secondaryFiles:
      - ^.bai
    'sbg:x': 1508.0369873046875
    'sbg:y': 1641.453125
  - id: unfiltered-bam
    outputSource:
      - bam_collapsing/unfiltered-bam
    type: File
    doc: >-
      This file is generated after collapsing of reads from the standard bam
      file. This is all duplex,simplex and sigletons as part of the alignment
    secondaryFiles:
      - ^.bai
    'sbg:x': 2190.28076171875
    'sbg:y': 661.25
  - id: simplex-bam
    outputSource:
      - bam_collapsing/simplex-bam
    type: File
    doc: >-
      This SIMPLEX BAM file is generated from Marianas SeparateBams module which
      seprate bam file based on duplex and simple clusters.
    secondaryFiles:
      - ^.bai
    'sbg:x': 2190.28076171875
    'sbg:y': 874.640625
  - id: second_pass_insertions
    outputSource:
      - bam_collapsing/second_pass_insertions
    type: File
    doc: >-
      This file containing inserstion is generated by Marianas
      DuplexUMIBamToCollapsedFastqSecondPass
    'sbg:x': 2190.28076171875
    'sbg:y': 981.296875
  - id: second_pass_alt_alleles
    outputSource:
      - bam_collapsing/second_pass_alt_alleles
    type: File
    doc: >-
      This file containing ALT ALLELES is generated by Marianas
      DuplexUMIBamToCollapsedFastqSecondPass
    'sbg:x': 2190.28076171875
    'sbg:y': 1087.953125
  - id: pileup_without_duplicates
    outputSource:
      - bam_collapsing/pileup_without_duplicates
    type: File
    'sbg:x': 2190.28076171875
    'sbg:y': 1194.609375
  - id: intervals_without_duplicates
    outputSource:
      - bam_collapsing/intervals_without_duplicates
    type: File
    'sbg:x': 2190.28076171875
    'sbg:y': 1301.265625
  - id: intervals
    outputSource:
      - bam_collapsing/intervals
    type: File
    'sbg:x': 2190.28076171875
    'sbg:y': 1407.921875
  - id: gzip_read1
    outputSource:
      - bam_collapsing/gzip_read1
    type: File
    doc: >-
      This is the collapsed READ1 gzip fastq file generated after MARIANAS
      collapsing
    'sbg:x': 2190.28076171875
    'sbg:y': 1621.3515625
  - id: gzip_read2
    outputSource:
      - bam_collapsing/gzip_read2
    type: File
    doc: >-
      This is the collapsed READ2 gzip fastq file generated after MARIANAS
      collapsing
    'sbg:x': 2190.28076171875
    'sbg:y': 1514.6171875
  - id: first_pass_insertions
    outputSource:
      - bam_collapsing/first_pass_insertions
    type: File
    doc: >-
      This file containing inserstion is generated by Marianas
      DuplexUMIBamToCollapsedFastqFirstPass
    'sbg:x': 2190.28076171875
    'sbg:y': 1728.046875
  - id: duplex-bam
    outputSource:
      - bam_collapsing/duplex-bam
    type: File
    doc: >-
      This DUPLEX BAM file is generated from Marianas SeparateBams module which
      seprate bam file based on duplex and simple clusters.
    secondaryFiles:
      - ^.bai
    'sbg:x': 2190.28076171875
    'sbg:y': 1834.703125
  - id: collapsed_fastq_2
    outputSource:
      - bam_collapsing/collapsed_fastq_2
    type: File
    doc: This is the collapsed READ2 fastq file generated after MARIANAS collapsing
    'sbg:x': 2190.28076171875
    'sbg:y': 1941.359375
  - id: alt_allele_file
    outputSource:
      - bam_collapsing/alt_allele_file
    type: File
    doc: >-
      This file containing ALT ALLELES is generated by Marianas
      DuplexUMIBamToCollapsedFastqFIRSTPASS
    'sbg:x': 2190.28076171875
    'sbg:y': 2154.671875
  - id: alignment_metrics_unfiltered
    outputSource:
      - bam_collapsing/alignment_metrics_unfiltered
    type: File
    doc: >-
      Alignment metrics TXT file generated by Picard CollectALignmentMetrics for
      Unfilered BAM File.
    'sbg:x': 2190.28076171875
    'sbg:y': 2261.3671875
  - id: alignment_metrics_simplex
    outputSource:
      - bam_collapsing/alignment_metrics_simplex
    type: File
    doc: >-
      Alignment metrics TXT file generated by Picard CollectALignmentMetrics for
      SIMPLEX BAM File.
    'sbg:x': 2190.28076171875
    'sbg:y': 2368.1015625
  - id: alignment_metrics_duplex
    outputSource:
      - bam_collapsing/alignment_metrics_duplex
    type: File
    doc: >-
      Alignment metrics TXT file generated by Picard CollectALignmentMetrics for
      DUPLEX BAM File.
    'sbg:x': 2190.28076171875
    'sbg:y': 2474.8359375
  - id: collapsed_fastq_1
    outputSource:
      - bam_collapsing/collapsed_fastq_1
    type: File
    doc: This is the collapsed READ1 fastq file generated after MARIANAS collapsing
    'sbg:x': 2190.28076171875
    'sbg:y': 2048.015625
  - id: standard_bam_indel_realign_targets
    outputSource:
      - standard_bam_processing_cwl/output_file
    type: File?
    'sbg:x': 1508.0369873046875
    'sbg:y': 1107.8984375
  - id: unfiltered_bam_indel_realigned_targets
    outputSource:
      - bam_collapsing/output_file
    type: File?
    'sbg:x': 2190.28076171875
    'sbg:y': 767.9453125
  - id: standard_bam_alignment_metrics
    outputSource:
      - standard_bam_processing_cwl/standard_bam_alignment_metrics
    type: File
    'sbg:x': 1508.0369873046875
    'sbg:y': 1214.6328125
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
    'sbg:x': 373.78125
    'sbg:y': 1835.015625
  - id: standard_bam_processing_cwl
    in:
      - id: fastq2
        source: marianas_process_loop_umi_cwl/processed_fastq_2
      - id: reference
        source: reference
      - id: known_sites_1
        source: known_sites_1
      - id: known_sites_2
        source: known_sites_2
      - id: validation_stringency
        default: LENIENT
        source: validation_stringency
      - id: create_bam_index
        default: true
      - id: assume_sorted
        default: true
        source: assume_sorted
      - id: bam_index
        default: true
      - id: option_bedgraph
        default: true
        source: option_bedgraph
      - id: fastq1
        source: marianas_process_loop_umi_cwl/processed_fastq_1
      - id: read_group_sequnecing_center
        default: MSKCC
        source: read_group_sequnecing_center
      - id: read_group_sequencing_platform
        default: ILLUMINA
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
        source: P
      - id: M
        default: true
        source: M
      - id: create_bam_index_1
        default: true
      - id: sort_order
        default: coordinate
        source: sort_order
      - id: output
        source: standard_aln_output_file_name
      - id: output_file_name
        source: standard_picard_addrg_output_filename
      - id: window_size
        default: '800,700'
        source: window_size
      - id: soft_clip_contig
        default: '100,30,80,15'
        source: soft_clip_contig
      - id: scoring_gap_alignments
        default: '8,32,48,1'
        source: scoring_gap_alignments
      - id: maximum_mixmatch_rate
        default: 0.1
        source: maximum_mixmatch_rate
      - id: maximum_average_depth
        default: 1000
        source: maximum_average_depth
      - id: ignore_bad_assembly
        default: true
        source: ignore_bad_assembly
      - id: contig_anchor
        default: '10,1'
        source: contig_anchor
      - id: consensus_sequence
        default: true
        source: consensus_sequence
      - id: stringency
        default: 3
        source: stringency
      - id: quality
        default: 1
        source: quality
      - id: length
        default: 25
        source: length
      - id: adapter2
        default: AGATCGGAAGAGC
        source: adapter2
      - id: adapter
        default: GATCGGAAGAGC
        source: adapter
      - id: bqsr_read_filter
        default: GoodCigarReadFilter
        source: bqsr_read_filter
      - id: number_of_threads
        default: 16
        source: number_of_threads
      - id: trim_galore_number_of_threads
        default: 4
        source: trim_galore_number_of_threads
    out:
      - id: clstats2
      - id: clstats1
      - id: bqsr_bam
      - id: md_bam
      - id: output_file
      - id: standard_bam_alignment_metrics
    run: standard_bam_processing/standard_bam_processing.cwl
    label: standard_bam_processing.cwl
    'sbg:x': 692.921875
    'sbg:y': 1148.9609375
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
        source: min_map_quality
      - id: min_base_quality
        default: 20
        source: min_base_quality
      - id: mismatches
        default: 0
        source: mismatches
      - id: wobble
        default: 1
        source: wobble
      - id: min_consensus_percent
        default: 90
        source: min_consensus_percent
      - id: key
        default:
          - '6,6n'
          - '8,8n'
        source:
          - key
      - id: sort_first_pass_output_file_name
        source: sort_first_pass_output_file_name
      - id: output_name_collapsed_gzip_R1
        source: output_name_collapsed_gzip_R1
      - id: output_name_collapsed_gzip_R2
        source: output_name_collapsed_gzip_R2
      - id: read_group_sequnecing_center
        default: MSKCC
        source: read_group_sequnecing_center
      - id: read_group_sequencing_platform
        default: ILLUMINA
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
        source: collapsing_picard_output_file_name
      - id: aln_output_file_name
        source: collapsing_aln_output_file_name
      - id: abra_collapsing_number_of_threads
        default: 16
        source: abra_collapsing_number_of_threads
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
    'sbg:x': 1508.0369873046875
    'sbg:y': 1888.1484375
requirements:
  - class: SubworkflowFeatureRequirement
'sbg:license': Apache Software License 2.0
'sbg:toolAuthor': 'Ronak Shah, Ian Johnson, Shalabh Suman'
