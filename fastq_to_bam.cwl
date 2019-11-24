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
  - id: reference
    type: File
    doc: >-
      The reference sequence in a single reference sequence in FASTA format,
      with all contigs in the same file, validated according to the FASTA
      standard. It has multiple secondary file associated with it ending in
      ".dict, .fai, .amb, .ann, .bwt, .pac, .sa"
    secondaryFiles:
      - ^.dict
      - .fai
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
    'sbg:x': 577.974853515625
    'sbg:y': 427.25
  - id: known_sites_1
    type: File
    doc: >-
      A database of known polymorphic sites on VCF format, Ex: DBSNP or
      Mills_and_1000G. Note: ".vcf.idx" secaondary file should be present where
      the ".vcf" file is located
    secondaryFiles:
      - .idx
    'sbg:x': 577.974853515625
    'sbg:y': 2285.0625
  - id: bed_file
    type: File
    doc: >-
      Targets in BED file format used by Waltz to generate the PileUp for
      collapsing of the BAM file.The genotype from positions in this bed file
      will be used as the consensus base if min_consensus_percent threshold is
      not reached. Otherwise, the reference base from the supplied
      reference_fasta will be used
    'sbg:x': 0
    'sbg:y': 2584.5
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
    'sbg:x': 577.974853515625
    'sbg:y': 747.6875
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
    'sbg:x': 577.974853515625
    'sbg:y': 854.5
  - id: read_group_library
    type: string
    doc: >-
      LB = DNA preparation library identifier

      MarkDuplicates uses the LB field to determine which read groups might
      contain molecular duplicates, in case the same DNA library was sequenced
      on multiple lanes.
    'sbg:x': 577.974853515625
    'sbg:y': 961.3125
  - id: read_group_identifier
    type: string
    doc: >-
      ID = Read group identifier

      This tag identifies which read group each read belongs to, so each read
      group's ID must be unique. It is referenced both in the read group
      definition line in the file header (starting with @RG) and in the RG:Z tag
      for each read record.
    'sbg:x': 577.974853515625
    'sbg:y': 1068.125
  - id: sort_first_pass_output_file_name
    type: string
    doc: Name for the Marianas Duplex Collapsing First Pass output TXT file.
    'sbg:x': 0
    'sbg:y': 1302.75
  - id: output_name_collapsed_gzip_R2
    type: string?
    doc: Name of the output collapsed READ1 gzip fastq file.
    'sbg:x': 0
    'sbg:y': 1409.5625
  - id: output_name_collapsed_gzip_R1
    type: string?
    doc: Name of the output collapsed READ1 gzip fastq file.
    'sbg:x': 0
    'sbg:y': 1516.375
  - id: collapsing_aln_output_file_name
    type: string?
    doc: Name of the SAM format output file created by bwa mem for collapsing step.
    'sbg:x': 0
    'sbg:y': 2477.6875
  - id: collapsing_picard_output_file_name
    type: string?
    doc: >-
      Name of the BAM format output file created by Picard
      AddOrReplaceReadGroups for collapsing step.
    'sbg:x': 0
    'sbg:y': 2370.875
  - id: standard_aln_output_file_name
    type: string?
    label: standard_aln_output_file_name
    doc: >-
      Name of the SAM format output file created by bwa mem for standard bam
      processing step
    'sbg:x': 577.974853515625
    'sbg:y': 0
  - id: standard_picard_addrg_output_filename
    type: string?
    label: standard_picard_addrg_output_filename
    doc: >-
      Name of the BAM format output file created by Picard
      AddOrReplaceReadGroups for standard bam processing step.
    'sbg:x': 0
    'sbg:y': 1195.9375
  - id: known_sites_2
    type: File?
    doc: >-
      A database of known polymorphic sites on VCF format, Ex: DBSNP or
      Mills_and_1000G. Note: ".vcf.idx" secaondary file should be present where
      the ".vcf" file is located
    secondaryFiles:
      - .idx
    'sbg:x': 577.974853515625
    'sbg:y': 2178.25
  - id: wobble
    type: int?
    doc: Allowable left and right shift amount for grouping UMI families
    'sbg:x': 0
    'sbg:y': 555.0625
  - id: read_group_sequencing_center
    type: string?
    doc: RGCN tag for BAM file indicating where the data is sequenced.
    'sbg:x': 577.974853515625
    'sbg:y': 640.875
  - id: read_group_sequencing_platform
    type: string?
    doc: BAM Tag describing the Platform used to generate the sequencing data.
    'sbg:x': 577.974853515625
    'sbg:y': 534.0625
  - id: mismatches
    type: int?
    doc: Allowable mismatch count in UMI bases for grouping UMI families
    'sbg:x': 0
    'sbg:y': 1623.1875
  - id: min_map_quality
    type: int?
    doc: Make sure to use high quality reads.
    'sbg:x': 0
    'sbg:y': 1730
  - id: min_consensus_percent
    type: int?
    doc: >-
      Percentage of bases that must be in agreement at each position in the
      consensus read before masking that base as "N
    'sbg:x': 0
    'sbg:y': 1836.8125
  - id: min_base_quality
    type: int?
    doc: Minimum Base Quality score to be used during collapsing.
    'sbg:x': 0
    'sbg:y': 1943.625
  - id: key
    type:
      - 'null'
      - type: array
        items: string
        inputBinding:
          prefix: '-k'
    label: sort key
    doc: >-
      sort via a key; KEYDEF gives location and type. KEYDEF  is 
      F[.C][OPTS][,F[.C][OPTS]]  for  start  and  stop  position, where F is a
      field number and C a character position in the field; both are origin 1,
      and the stop position defaults to the line's end.  If neither -t nor -b 
      is in effect, characters in a field are counted from the beginning of the
      preceding whitespace.  OPTS is one or more  single-letter ordering
      options  [bdfgiMhnRrV], which override global ordering options for that
      key.   If no key is given, use the entire line as the key.
    'sbg:x': 0
    'sbg:y': 2050.4375
  - id: adapter
    type: string?
    label: Adapter for READ1 for trim_galore
    doc: Adapter sequence to trim READ1.
    'sbg:x': 577.974853515625
    'sbg:y': 3139.5625
  - id: adapter2
    type: string?
    label: Adapter for READ2 for trim_galore
    doc: Adapter sequence to trim READ2.
    'sbg:x': 577.974853515625
    'sbg:y': 3032.75
  - id: assume_sorted
    type: boolean?
    doc: Assume that the given bam file is coordinate sorted for picard tools
    'sbg:x': 577.974853515625
    'sbg:y': 2925.9375
  - id: bqsr_read_filter
    type: 'string[]?'
    doc: GATK READ_FILTER option to apply defferent set of ReadFilter
    'sbg:x': 577.974853515625
    'sbg:y': 2819.125
  - id: consensus_sequence
    type: boolean?
    doc: Use positional consensus sequence when aligning high quality soft clipping
    'sbg:x': 577.974853515625
    'sbg:y': 2712.3125
  - id: contig_anchor
    type: string?
    doc: >-
      Contig anchor [M_bases_at_contig_edge,max_mismatches_near_edge]
      (default:10,2)
    'sbg:x': 577.974853515625
    'sbg:y': 2605.5
  - id: option_bedgraph
    type: boolean?
    doc: >-
      Report depth in BedGraph format. For details, see:
      http://genome.ucsc.edu/goldenPath/help/bedgraph.html
    'sbg:x': 577.974853515625
    'sbg:y': 1388.5625
  - id: number_of_threads
    type: int?
    label: abra_number_of_threads
    doc: Number of threads for parallel exectution of ABRA
    'sbg:x': 577.974853515625
    'sbg:y': 1495.375
  - id: maximum_mixmatch_rate
    type: float?
    doc: >-
      Max allowed mismatch rate when mapping reads back to contigs
      (default:0.05)
    'sbg:x': 577.974853515625
    'sbg:y': 1602.1875
  - id: maximum_average_depth
    type: int?
    doc: >-
      Regions with average depth exceeding this value will be downsampled
      (default: 1000)
    'sbg:x': 577.974853515625
    'sbg:y': 1709
  - id: M
    type: boolean?
    label: BWA mark shorter split hits as secondary
    doc: mark shorter split hits as secondary (for Picard/GATK compatibility)
    'sbg:x': 577.974853515625
    'sbg:y': 1964.625
  - id: length
    type: int?
    label: trim_galore minimum length for read
    doc: Trim_galore minimum length for read
    'sbg:x': 577.974853515625
    'sbg:y': 2071.4375
  - id: P
    type: boolean?
    label: BWA skip pairing
    doc: skip pairing; mate rescue performed unless -S also in use
    'sbg:x': 577.974853515625
    'sbg:y': 1281.75
  - id: ignore_bad_assembly
    type: boolean?
    doc: Use this option to avoid parsing errors for corrupted assemblies
    'sbg:x': 577.974853515625
    'sbg:y': 2391.875
  - id: scoring_gap_alignments
    type: string?
    doc: >-
      Scoring used for contig
      alignments(match,mismatch_penalty,gap_open_penalty, gap_extend_penalty
      (default:8,32,48,1)
    'sbg:x': 577.974853515625
    'sbg:y': 320.4375
  - id: soft_clip_contig
    type: string?
    doc: >-
      Soft clip contig args [max_contigs,min_base_qual,frac_high_qual_bases,
      min_soft_clip_len (default:16,13,80,15)
    'sbg:x': 577.974853515625
    'sbg:y': 213.625
  - id: window_size
    type: string?
    doc: 'Processing window size and overlap (size,overlap) (default: 400,200)'
    'sbg:x': 0
    'sbg:y': 661.875
  - id: validation_stringency
    type: string?
    doc: Picard Validation Stringency while running Picard Tools
    'sbg:x': 0
    'sbg:y': 768.6875
  - id: trim_galore_number_of_threads
    type: int?
    doc: Number of threads to run Trim Galore with Cutadapt
    'sbg:x': 0
    'sbg:y': 875.5
  - id: stringency
    type: int?
    label: trim_galore overlap stringency
    doc: >-
      Overlap with adapter sequence required to trim a sequence. Defaults to a
      very stringent setting of '1', i.e. even a single bp of overlapping
      sequence will be trimmed of the 3' end of any read.
    'sbg:x': 0
    'sbg:y': 1089.125
  - id: sort_order
    type: string?
    doc: How the BAM file should be sorted (default to coordinate)
    'sbg:x': 577.974853515625
    'sbg:y': 106.8125
  - id: quality
    type: int?
    label: trim_galore base quality
    doc: trim_galore quality value for trimming
    'sbg:x': 577.974853515625
    'sbg:y': 1174.9375
  - id: create_bam_index
    type: boolean?
    'sbg:x': 577.974853515625
    'sbg:y': 2498.6875
  - id: fastq1
    type:
      type: array
      items: File
      inputBinding:
        prefix: '--fastq1'
    doc: >-
      Full path to gziped READ1 fastq files, can be specified multiple times,
      please make sure that order between the FASTQ1 and FASTQ2 is always
      maintained
    'sbg:x': 0
    'sbg:y': 2264.0625
  - id: fastq2
    type:
      type: array
      items: File
      inputBinding:
        prefix: '--fastq2'
    doc: >-
      Full path to gziped READ2 fastq files, can be specified multiple times,
      please make sure that order between the FASTQ1 and FASTQ2 is always
      maintained
    'sbg:x': 0
    'sbg:y': 2157.25
  - id: temporary_directory
    type: string?
    label: picard_tools_tmpdir
    'sbg:x': 0
    'sbg:y': 982.3125
outputs:
  - id: composite_umi_frequencies
    outputSource:
      - marianas_process_loop_umi_cwl/composite_umi_frequencies
    type: File
    doc: >-
      This is text file consisting of frequencines of unique molecular
      identifier as seen by Marianas ProcessLoopUMIFastq
    'sbg:x': 897.131103515625
    'sbg:y': 1569.78125
  - id: clipping_info
    outputSource:
      - marianas_process_loop_umi_cwl/clipping_info
    type: File
    doc: >-
      File having information about all the clipped unique molecular identifiers
      from the fastq.gz files by Marianas ProcessLoopUMIFastq
    'sbg:x': 897.131103515625
    'sbg:y': 1676.59375
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
    'sbg:x': 1696.256103515625
    'sbg:y': 1245.875
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
    'sbg:x': 1696.256103515625
    'sbg:y': 1566.3125
  - id: unfiltered-bam
    outputSource:
      - bam_collapsing/unfiltered-bam
    type: File
    doc: >-
      This file is generated after collapsing of reads from the standard bam
      file. This is all duplex,simplex and sigletons as part of the alignment
    secondaryFiles:
      - ^.bai
    'sbg:x': 2527.21142578125
    'sbg:y': 661.875
  - id: simplex-bam
    outputSource:
      - bam_collapsing/simplex-bam
    type: File
    doc: >-
      This SIMPLEX BAM file is generated from Marianas SeparateBams module which
      seprate bam file based on duplex and simple clusters.
    secondaryFiles:
      - ^.bai
    'sbg:x': 2527.21142578125
    'sbg:y': 875.5
  - id: second_pass_insertions
    outputSource:
      - bam_collapsing/second_pass_insertions
    type: File
    doc: >-
      This file containing inserstion is generated by Marianas
      DuplexUMIBamToCollapsedFastqSecondPass
    'sbg:x': 2527.21142578125
    'sbg:y': 982.3125
  - id: second_pass_alt_alleles
    outputSource:
      - bam_collapsing/second_pass_alt_alleles
    type: File
    doc: >-
      This file containing ALT ALLELES is generated by Marianas
      DuplexUMIBamToCollapsedFastqSecondPass
    'sbg:x': 2527.21142578125
    'sbg:y': 1089.125
  - id: pileup_without_duplicates
    outputSource:
      - bam_collapsing/pileup_without_duplicates
    type: File
    'sbg:x': 2527.21142578125
    'sbg:y': 1195.9375
  - id: intervals_without_duplicates
    outputSource:
      - bam_collapsing/intervals_without_duplicates
    type: File
    'sbg:x': 2527.21142578125
    'sbg:y': 1302.75
  - id: intervals
    outputSource:
      - bam_collapsing/intervals
    type: File
    'sbg:x': 2527.21142578125
    'sbg:y': 1409.5625
  - id: gzip_read1
    outputSource:
      - bam_collapsing/gzip_read1
    type: File
    doc: >-
      This is the collapsed READ1 gzip fastq file generated after MARIANAS
      collapsing
    'sbg:x': 2527.21142578125
    'sbg:y': 1623.1875
  - id: gzip_read2
    outputSource:
      - bam_collapsing/gzip_read2
    type: File
    doc: >-
      This is the collapsed READ2 gzip fastq file generated after MARIANAS
      collapsing
    'sbg:x': 2527.21142578125
    'sbg:y': 1516.375
  - id: first_pass_insertions
    outputSource:
      - bam_collapsing/first_pass_insertions
    type: File
    doc: >-
      This file containing inserstion is generated by Marianas
      DuplexUMIBamToCollapsedFastqFirstPass
    'sbg:x': 2527.21142578125
    'sbg:y': 1730
  - id: duplex-bam
    outputSource:
      - bam_collapsing/duplex-bam
    type: File
    doc: >-
      This DUPLEX BAM file is generated from Marianas SeparateBams module which
      seprate bam file based on duplex and simple clusters.
    secondaryFiles:
      - ^.bai
    'sbg:x': 2527.21142578125
    'sbg:y': 1836.8125
  - id: collapsed_fastq_2
    outputSource:
      - bam_collapsing/collapsed_fastq_2
    type: File
    doc: This is the collapsed READ2 fastq file generated after MARIANAS collapsing
    'sbg:x': 2527.21142578125
    'sbg:y': 1943.625
  - id: alt_allele_file
    outputSource:
      - bam_collapsing/alt_allele_file
    type: File
    doc: >-
      This file containing ALT ALLELES is generated by Marianas
      DuplexUMIBamToCollapsedFastqFIRSTPASS
    'sbg:x': 2527.21142578125
    'sbg:y': 2157.25
  - id: alignment_metrics_unfiltered
    outputSource:
      - bam_collapsing/alignment_metrics_unfiltered
    type: File
    doc: >-
      Alignment metrics TXT file generated by Picard CollectALignmentMetrics for
      Unfilered BAM File.
    'sbg:x': 2527.21142578125
    'sbg:y': 2264.0625
  - id: alignment_metrics_simplex
    outputSource:
      - bam_collapsing/alignment_metrics_simplex
    type: File
    doc: >-
      Alignment metrics TXT file generated by Picard CollectALignmentMetrics for
      SIMPLEX BAM File.
    'sbg:x': 2527.21142578125
    'sbg:y': 2370.875
  - id: alignment_metrics_duplex
    outputSource:
      - bam_collapsing/alignment_metrics_duplex
    type: File
    doc: >-
      Alignment metrics TXT file generated by Picard CollectALignmentMetrics for
      DUPLEX BAM File.
    'sbg:x': 2527.21142578125
    'sbg:y': 2477.6875
  - id: collapsed_fastq_1
    outputSource:
      - bam_collapsing/collapsed_fastq_1
    type: File
    doc: This is the collapsed READ1 fastq file generated after MARIANAS collapsing
    'sbg:x': 2527.21142578125
    'sbg:y': 2050.4375
  - id: standard_bam_indel_realign_targets
    outputSource:
      - standard_bam_processing_cwl/output_file
    type: File?
    'sbg:x': 1696.256103515625
    'sbg:y': 1032.25
  - id: unfiltered_bam_indel_realigned_targets
    outputSource:
      - bam_collapsing/output_file
    type: File?
    'sbg:x': 2527.21142578125
    'sbg:y': 768.6875
  - id: standard_bam_alignment_metrics
    outputSource:
      - standard_bam_processing_cwl/standard_bam_alignment_metrics
    type: File
    'sbg:x': 1696.256103515625
    'sbg:y': 1139.0625
  - id: clstats1
    outputSource:
      - standard_bam_processing_cwl/clstats1
    type: File
    'sbg:x': 1696.256103515625
    'sbg:y': 1459.5
  - id: clstats2
    outputSource:
      - standard_bam_processing_cwl/clstats2
    type: File
    'sbg:x': 1696.256103515625
    'sbg:y': 1352.6875
steps:
  - id: marianas_process_loop_umi_cwl
    in:
      - id: fastq1
        source: merge_fastq_0_1_7/mergedfastq1
      - id: fastq2
        source: merge_fastq_0_1_7/mergedfastq2
      - id: umi_length
        default: 3
    out:
      - id: processed_fastq_1
      - id: processed_fastq_2
      - id: clipping_info
      - id: composite_umi_frequencies
    run: >-
      command_line_tools/marianas_process_loop_umi_1.8.1/marianas_process_loop_umi.cwl
    label: Loop UMI Processing
    doc: Remove Loop UMI from the reads and add them to Read Names
    'sbg:x': 577.974853515625
    'sbg:y': 1836.8125
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
      - id: option_bedgraph
        default: true
        source: option_bedgraph
      - id: fastq1
        source: marianas_process_loop_umi_cwl/processed_fastq_1
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
      - id: number_of_threads
        default: 16
        source: number_of_threads
      - id: validation_stringency
        default: LENIENT
        source: validation_stringency
      - id: create_bam_index
        default: true
        source: create_bam_index
      - id: assume_sorted
        default: true
        source: assume_sorted
      - id: M
        default: true
        source: M
      - id: sort_order
        default: coordinate
        source: sort_order
      - id: trim_galore_number_of_threads
        default: 4
        source: trim_galore_number_of_threads
      - id: read_filter
        default:
          - GoodCigarReadFilter
        source:
          - bqsr_read_filter
      - id: read_group_sequencing_center
        default: MSKCC
        source: read_group_sequencing_center
      - id: temporary_directory
        source: temporary_directory
    out:
      - id: bqsr_bam
      - id: md_bam
      - id: output_file
      - id: standard_bam_alignment_metrics
      - id: clstats1
      - id: clstats2
    run: standard_bam_processing/standard_bam_processing.cwl
    label: Best Practices for BAM Generation
    doc: >-
      Using Trimming, Alignment, MarkDuplicate, Realignment and Recalibration to
      generate standard bam file.
    'sbg:x': 897.131103515625
    'sbg:y': 1210.875
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
      - id: sort_order
        default: coordinate
      - id: create_bam_index
        default: true
        source: create_bam_index
      - id: P
        default: true
        source: P
      - id: window_size
        source: window_size
      - id: scoring_gap_alignments
        source: scoring_gap_alignments
      - id: option_bedgraph
        default: true
        source: option_bedgraph
      - id: maximum_average_depth
        source: maximum_average_depth
      - id: maximum_mixmatch_rate
        source: maximum_mixmatch_rate
      - id: ignore_bad_assembly
        default: true
        source: ignore_bad_assembly
      - id: contig_anchor
        source: contig_anchor
      - id: consensus_sequence
        default: true
        source: consensus_sequence
      - id: abra_collapsing_number_of_threads
        default: 12
        source: number_of_threads
      - id: wobble
        default: 1
        source: wobble
      - id: temporary_directory
        source: temporary_directory
      - id: read_group_sequencing_center
        default: MSKCC
        source: read_group_sequencing_center
    out:
      - id: second_pass_insertions
      - id: second_pass_alt_alleles
      - id: collapsed_fastq_2
      - id: collapsed_fastq_1
      - id: pileup_without_duplicates
      - id: intervals_without_duplicates
      - id: intervals
      - id: unfiltered-bam
      - id: output_file
      - id: first_pass_insertions
      - id: alt_allele_file
      - id: first_pass_output_dir
      - id: gzip_read1
      - id: gzip_read2
      - id: simplex-bam
      - id: duplex-bam
      - id: alignment_metrics_unfiltered
      - id: alignment_metrics_simplex
      - id: alignment_metrics_duplex
    run: bam_collapsing/bam_collapsing.cwl
    label: Collapsing reads for error supression
    doc: >-
      Using Marianas to cluster and collapse reads generating unfiltered,
      simplex and duplex BAM files
    'sbg:x': 1696.256103515625
    'sbg:y': 1890.21875
  - id: merge_fastq_0_1_7
    in:
      - id: fastq1
        source:
          - fastq1
      - id: fastq2
        source:
          - fastq2
    out:
      - id: mergedfastq1
      - id: mergedfastq2
    run: command_line_tools/merge_fastq_0.1.7/merge_fastq_0.1.7.cwl
    label: Merge FASTQ.gz
    doc: >-
      Given multiple pair-end fastq data merge them into single pair-end fastq
      w.r.t each READ1 and READ2
    'sbg:x': 373.78125
    'sbg:y': 1562.78125
requirements:
  - class: SubworkflowFeatureRequirement
'sbg:license': Apache Software License 2.0
'sbg:toolAuthor': 'Ronak Shah, Ian Johnson, Shalabh Suman'
