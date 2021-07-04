{
    "$graph": [
        {
            "class": "Workflow",
            "id": "#bam_collapsing.cwl",
            "doc": "This is the workflow is written using Common Workflow Language (CWL) version 1.0 (https://www.commonwl.org/v1.0/) and is used at Memorial Sloan Kettering Cancer Center for producing collapsed bam files from data generated for the NY state-approved MSK-ACCESS assay.",
            "label": "bam_collapsing",
            "inputs": [
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_input",
                    "type": "File",
                    "doc": "Fgbio GroupReadsByUmi: The input BAM file.",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2774.9375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_strategy",
                    "type": "string",
                    "doc": "Fgbio GroupReadsByUmi: The UMI assignment strategy.(identity, edit, adjacency, paired)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2348
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_raw_tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio GroupReadsByUmi: The tag containing the raw UMI. (Default:RX)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2454.734375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio GroupReadsByUmi: The output BAM file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2561.46875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_min_umi_length",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio GroupReadsByUmi: The minimum UMI length. If not specified then all UMIs must have the same length, otherwise discard reads with UMIs shorter than this length and allow for differing UMI lengths.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2668.203125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_include_non_pf_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Fgbio GroupReadsByUmi: Include non-PF reads.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2881.671875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_family_size_histogram",
                    "type": "string",
                    "doc": "Fgbio GroupReadsByUmi: Optional output of tag family size counts.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2988.40625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_edits",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio GroupReadsByUmi:  The allowable number of edits between UMIs.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3095.140625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_assign_tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio GroupReadsByUmi: The output tag for UMI grouping. (Default:MI)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3201.875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_intervals",
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "Fgbio CollectDuplexSeqMetrics:  Optional set of intervals over which to restrict analysis.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3735.546875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_output_prefix",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio CollectDuplexSeqMetrics: Prefix of output files to write.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3308.609375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_min_ba_reads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio CollectDuplexSeqMetrics: Minimum BA reads to call a tag family a \u2018duplex\u2019.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3415.34375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_min_ab_reads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio CollectDuplexSeqMetrics: Minimum AB reads to call a tag family a \u2018duplex\u2019.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3522.078125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_mi_tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio CollectDuplexSeqMetrics: The output tag for UMI grouping. (Default: MI)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3628.8125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_duplex_umi_counts",
                    "type": "boolean",
                    "doc": "Fgbio CollectDuplexSeqMetrics: If true, produce the .duplex_umi_counts.txt file with counts of duplex UMI observations.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3842.28125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_description",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio CollectDuplexSeqMetrics:Description of data set used to label plots. Defaults to sample/library.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3949.015625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_trim",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Fgbio CallDuplexConsensusReads : If true, quality trim input reads in addition to masking low Q bases.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4055.75
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio CallDuplexConsensusReads: The sort order of the output, if :none: then the same as the input.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4162.484375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_read_name_prefix",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio CallDuplexConsensusReads: The prefix all consensus read names",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4269.21875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_read_group_id",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio CallDuplexConsensusReads: The new read group ID for all the consensus reads.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4375.953125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio CallDuplexConsensusReads: Output SAM or BAM file to write consensus reads.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4482.6875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_min_reads",
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "doc": "Fgbio CallDuplexConsensusReads: The minimum number of input reads to a consensus read.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4589.421875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_min_input_base_quality",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio CallDuplexConsensusReads: Ignore bases in raw reads that have Q below this value.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4696.15625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_max_reads_per_strand",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio CallDuplexConsensusReads: The maximum number of reads to use when building a single-strand consensus. If more than this many reads are present in a tag family, the family is randomly downsampled to exactly max-reads reads.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4802.890625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_error_rate_pre_umi",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio CallDuplexConsensusReads: The Phred-scaled error rate for an error prior to the UMIs being integrated.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4909.625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_error_rate_post_umi",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio CallDuplexConsensusReads: The Phred-scaled error rate for an error post the UMIs have been integrated.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5016.359375
                },
                {
                    "id": "#bam_collapsing.cwl/reference_sequence",
                    "type": "File",
                    "doc": "Reference sequence file. Please include \".fai\", \"^.dict\", \".amb\" , \".sa\", \".bwt\", \".pac\" as secondary files if they are not present in the same location as the \".fasta\" file",
                    "secondaryFiles": [
                        ".fai",
                        "^.dict",
                        ".amb",
                        ".sa",
                        ".ann",
                        ".bwt",
                        ".pac"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 320.15625
                },
                {
                    "id": "#bam_collapsing.cwl/validation_stringency",
                    "type": "string",
                    "doc": "GATK Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 0
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_output_name_unpaired",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "GATK SamToFastq: unpaired fastq output file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1600.859375
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_output_name_R2",
                    "type": "string",
                    "doc": "GATK SamToFastq: R2 fastq.gz output file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1707.59375
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_include_non_primary_alignments",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "GATK SamToFastq: If true, include non-primary alignments in the output. Support of non-primary alignments in SamToFastq is not comprehensive, so there may be exceptions if this is set to true and there are paired reads with non-primary alignments.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1921.0625
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_include_non_pf_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "GATK SamToFastq: Include non-PF reads from the SAM file into the output FASTQ files. PF means 'passes filtering'. Reads whose 'not passing quality controls' flag is set are non-PF reads. See GATK Dictionary for more info.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2027.796875
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_output_name_R1",
                    "type": "string",
                    "doc": "GATK SamToFastq: R1 fastq.gz output file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1814.328125
                },
                {
                    "id": "#bam_collapsing.cwl/bwa_mem_Y",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "bwa mem: to force soft-clipping rather than default hard-clipping of supplementary alignments",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5229.78125
                },
                {
                    "id": "#bam_collapsing.cwl/bwa_mem_T",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "bwa mem: Don\u2019t output alignment with score lower than INT. This option only affects output.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5336.484375
                },
                {
                    "id": "#bam_collapsing.cwl/sort_order",
                    "type": "string",
                    "doc": "The order in which the merged reads should be output.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 213.453125
                },
                {
                    "id": "#bam_collapsing.cwl/picard_addRG_read_group_sequencing_platform",
                    "type": "string",
                    "doc": "PicardAddOrReplaceReadGroups: Read-Group platform (e.g. ILLUMINA, SOLID)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 640.28125
                },
                {
                    "id": "#bam_collapsing.cwl/picard_addRG_read_group_sequencing_center",
                    "type": "string",
                    "doc": "PicardAddOrReplaceReadGroups: Read-Group sequencing center name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 747.015625
                },
                {
                    "id": "#bam_collapsing.cwl/picard_addRG_read_group_run_date",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "PicardAddOrReplaceReadGroups: Read-Group Run date",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 960.484375
                },
                {
                    "id": "#bam_collapsing.cwl/picard_addRG_read_group_platform_unit",
                    "type": "string",
                    "doc": "PicardAddOrReplaceReadGroups: Read-Group Platform Unit (eg. run barcode)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1067.21875
                },
                {
                    "id": "#bam_collapsing.cwl/picard_addRG_read_group_library",
                    "type": "string",
                    "doc": "PicardAddOrReplaceReadGroups: Read-Group library",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1173.953125
                },
                {
                    "id": "#bam_collapsing.cwl/picard_addRG_read_group_identifier",
                    "type": "string",
                    "doc": "PicardAddOrReplaceReadGroups: Read-Group ID",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1280.6875
                },
                {
                    "id": "#bam_collapsing.cwl/picard_addRG_read_group_description",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "PicardAddOrReplaceReadGroups: Read-Group Description",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1387.421875
                },
                {
                    "id": "#bam_collapsing.cwl/picard_addRG_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "PicardAddOrReplaceReadGroups output file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1494.140625
                },
                {
                    "id": "#bam_collapsing.cwl/bwa_mem_P",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "bwa mem : In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5443.1875
                },
                {
                    "id": "#bam_collapsing.cwl/bwa_mem_output",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "bwa mem: Output SAM file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5549.890625
                },
                {
                    "id": "#bam_collapsing.cwl/create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 2248.361328125,
                    "https://www.sevenbridges.com/y": 2901.7421875
                },
                {
                    "id": "#bam_collapsing.cwl/bwa_mem_M",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Mark shorter split hits as secondary",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5656.59375
                },
                {
                    "id": "#bam_collapsing.cwl/bwa_mem_K",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "bwa mem : to achieve deterministic alignment results (Note: this is a hidden option)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5763.296875
                },
                {
                    "id": "#bam_collapsing.cwl/abra2_window_size",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Processing window size and overlap (size,overlap) (default: 400,200)",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 2614.8671875
                },
                {
                    "id": "#bam_collapsing.cwl/abra2_soft_clip_contig",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Soft clip contig args [max_contigs,min_base_qual,frac_high_qual_bases,min_soft_clip_len] (default:16,13,80,15)",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 2721.5859375
                },
                {
                    "id": "#bam_collapsing.cwl/abra2_scoring_gap_alignments",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Scoring used for contig alignments(match, mismatch_penalty,gap_open_penalty,gap_extend_penalty) (default:8,32,48,1)",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 2828.3203125
                },
                {
                    "id": "#bam_collapsing.cwl/picard_fixmate_information_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "PicardFixMateInformaiton: The output file to write to",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 426.859375
                },
                {
                    "id": "#bam_collapsing.cwl/abra2_output_bams",
                    "type": [
                        "string",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "doc": "Required list of output sam or bam file (s) separated by comma",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 2935.0390625
                },
                {
                    "id": "#bam_collapsing.cwl/bedtools_genomecov_option_bedgraph",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "option flag parameter to choose output file format. -bg refers to bedgraph format",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 2508.1484375
                },
                {
                    "id": "#bam_collapsing.cwl/abra2_no_sort",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Do not attempt to sort final output",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 3041.7421875
                },
                {
                    "id": "#bam_collapsing.cwl/abra2_no_edge_complex_indel",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Prevent output of complex indels at read start or read end",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 3148.4609375
                },
                {
                    "id": "#bam_collapsing.cwl/abra2_maximum_mixmatch_rate",
                    "type": [
                        "null",
                        "float"
                    ],
                    "doc": "Max allowed mismatch rate when mapping reads back to contigs (default: 0.05)",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 3255.1796875
                },
                {
                    "id": "#bam_collapsing.cwl/abra2_maximum_average_depth",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Regions with average depth exceeding this value will be downsampled (default: 1000)",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 3361.8984375
                },
                {
                    "id": "#bam_collapsing.cwl/bedtools_merge_distance_between_features",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Maximum distance between features allowed for features to be merged.",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 2401.4140625
                },
                {
                    "id": "#bam_collapsing.cwl/abra2_contig_anchor",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Contig anchor [M_bases_at_contig_edge,max_mismatches_near_edge] (default:10,2)",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 3468.6328125
                },
                {
                    "id": "#bam_collapsing.cwl/abra2_consensus_sequence",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Use positional consensus sequence when aligning high quality soft clipping",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 3575.3515625
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "GATK MergeBamAlignment: Output File Name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2134.53125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Fgbio FilterConsensusReads: Reverse [complement] per base tags on reverse  strand reads.- Simplex+Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 1974.3984375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_reverse_per_base_tags_duplex",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Fgbio FilterConsensusReads: Reverse [complement] per base tags on reverse  strand reads. - Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 2081.1328125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_require_single_strand_agreement_simplex_duplex",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Fgbio FilterConsensusReads: Mask (make N) consensus bases where the AB and BA consensus reads disagree (for duplex-sequencing only).",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 2187.8671875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_require_single_strand_agreement_duplex",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Fgbio FilterConsensusReads: Mask (make N) consensus bases where the AB and BA consensus reads disagree (for duplex-sequencing only).",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 2294.6015625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_base_error_rate_duplex",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "float"
                        }
                    ],
                    "doc": "Fgbio FilterConsensusReads: The maximum error rate for a single consensus  base. (Max 3 values) - Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 4109.0859375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_base_error_rate_simplex_duplex",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "float"
                        }
                    ],
                    "doc": "Fgbio FilterConsensusReads: The maximum error rate for a single consensus  base. (Max 3 values) - Simplex + Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 4002.3515625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_no_call_fraction_duplex",
                    "type": [
                        "null",
                        "float"
                    ],
                    "doc": "Fgbio FilterConsensusReads: Maximum fraction of no-calls in the read  after filtering - Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 3895.6171875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_read_error_rate_duplex",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "float"
                        }
                    ],
                    "doc": "Fgbio FilterConsensusReads: The maximum raw-read error rate across the  entire consensus read. (Max 3 values) - Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 3682.1484375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_no_call_fraction_simplex_duplex",
                    "type": [
                        "null",
                        "float"
                    ],
                    "doc": "Fgbio FilterConsensusReads: Maximum fraction of no-calls in the read  after filtering - Simplex + Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 3788.8828125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_read_error_rate_simplex_duplex",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "float"
                        }
                    ],
                    "doc": "The maximum raw-read error rate across the entire consensus read. (Max 3 values) - Simplex + Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 3575.4140625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_base_quality_duplex",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio FilterConsensusReads: Mask (make N) consensus bases with quality  less than this threshold. - Dupelx",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 3468.6796875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_base_quality_simplex_duplex",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio FilterConsensusReads: Mask (make N) consensus bases with quality  less than this threshold. - Simplex + Dupelx",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 3361.9453125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_mean_base_quality_duplex",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio FilterConsensusReads: The minimum mean base quality across the  consensus read - Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 3255.2109375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_mean_base_quality_simplex_duplex",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio FilterConsensusReads: The minimum mean base quality across the  consensus read - Simplex + Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 3148.4765625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_reads_duplex",
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "doc": "Fgbio FilterConsensusReads: The minimum number of reads supporting a  consensus base/read. (Max 3 values) - Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 3041.7421875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_reads_simplex_duplex",
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "doc": "Fgbio FilterConsensusReads: The minimum number of reads supporting a  consensus base/read. (Max 3 values) - Simplex+Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 2935.0078125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_output_file_name_simplex_duplex",
                    "type": "string",
                    "doc": "Fgbio FilterConsensusReads: Output file name Simplex + Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 2401.3359375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio FilterConsensusReads: Output file name Simplex alignment metrics",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 2508.0703125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_postprocessing_output_file_name_simplex",
                    "type": "string",
                    "doc": "Fgbio Postprocessing: Output file name Simplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 1867.6640625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio FilterConsensusReads: Output file name Duplex alignment metrics",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 2614.8046875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_output_file_name_duplex",
                    "type": "string",
                    "doc": "Fgbio FilterConsensusReads: Output file name Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 2721.5390625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_simplex_reads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio FilterConsensusReads: The minimum number of reads supporting a  consensus base/read. (Max 3 values) - Simplex+Duplex",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 2828.2734375
                },
                {
                    "id": "#bam_collapsing.cwl/picard_addRG_read_group_sample_name",
                    "type": "string",
                    "doc": "PicardAddOrReplaceReadGroups: Read-Group sample name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 853.75
                },
                {
                    "id": "#bam_collapsing.cwl/bwa_number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "bwa mem: Number of threads",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5123.078125
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "GATK CollectAlignmentSummaryMetrics: Output file name for metrics on collapsed BAM (Duplex+Simplex+Singletons)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2241.265625
                },
                {
                    "id": "#bam_collapsing.cwl/picard_addRG_sort_order",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 533.5625
                },
                {
                    "id": "#bam_collapsing.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 106.734375
                },
                {
                    "id": "#bam_collapsing.cwl/async_io",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5870.015625
                }
            ],
            "outputs": [
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_histogram",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/fgbio_group_reads_by_umi_histogram"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1079.60107421875,
                    "https://www.sevenbridges.com/y": 2627.90625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_umi_counts",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_umi_counts"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1680.5753173828125,
                    "https://www.sevenbridges.com/y": 2721.5390625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_family_size",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_family_size"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1680.5753173828125,
                    "https://www.sevenbridges.com/y": 2828.2734375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_duplex_yield_metrics",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_duplex_yield_metrics"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1680.5753173828125,
                    "https://www.sevenbridges.com/y": 2935.0078125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_duplex_umi_counts"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1680.5753173828125,
                    "https://www.sevenbridges.com/y": 3041.7421875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_duplex_qc",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_duplex_qc"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1680.5753173828125,
                    "https://www.sevenbridges.com/y": 3148.4765625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_duplex_family_size",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_duplex_family_size"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1680.5753173828125,
                    "https://www.sevenbridges.com/y": 3255.2109375
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_unpaired_fastq",
                    "outputSource": [
                        "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_unpaired_fastq"
                    ],
                    "type": [
                        "null",
                        "File"
                    ],
                    "https://www.sevenbridges.com/x": 2248.361328125,
                    "https://www.sevenbridges.com/y": 2581.5546875
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_second_end_fastq",
                    "outputSource": [
                        "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_second_end_fastq"
                    ],
                    "type": [
                        "null",
                        "File"
                    ],
                    "https://www.sevenbridges.com/x": 2248.361328125,
                    "https://www.sevenbridges.com/y": 2688.2890625
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_fastq",
                    "outputSource": [
                        "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_fastq"
                    ],
                    "type": [
                        "null",
                        "File"
                    ],
                    "https://www.sevenbridges.com/x": 2248.361328125,
                    "https://www.sevenbridges.com/y": 2795.0234375
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_txt_simplex",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_separate_bams/gatk_collect_alignment_summary_metrics_txt_simplex"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 5103.671875,
                    "https://www.sevenbridges.com/y": 2668.171875
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_txt_duplex",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_separate_bams/gatk_collect_alignment_summary_metrics_txt_duplex"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 5103.671875,
                    "https://www.sevenbridges.com/y": 2774.90625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_postprocessing_simplex_bam",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_separate_bams/fgbio_postprocessing_simplex_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 5103.671875,
                    "https://www.sevenbridges.com/y": 2988.375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_reads_duplex_bam",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_separate_bams/fgbio_filter_consensus_reads_duplex_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 5103.671875,
                    "https://www.sevenbridges.com/y": 3201.84375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collapsed_bam",
                    "outputSource": [
                        "#bam_collapsing.cwl/indel_realignment/indel_realignment_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 4160.3173828125,
                    "https://www.sevenbridges.com/y": 3230.7421875
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_txt_collapsed",
                    "outputSource": [
                        "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_4_1_8_0/gatk_collect_alignment_summary_metrics_txt"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 5103.671875,
                    "https://www.sevenbridges.com/y": 2881.640625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_bam",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/fgbio_group_reads_by_umi_bam"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1079.60107421875,
                    "https://www.sevenbridges.com/y": 2734.640625
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_filter_consensus_reads_simplex_duplex_bam",
                    "outputSource": [
                        "#bam_collapsing.cwl/fgbio_separate_bams/fgbio_filter_consensus_reads_simplex_duplex_bam"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 5103.671875,
                    "https://www.sevenbridges.com/y": 3095.109375
                }
            ],
            "steps": [
                {
                    "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0",
                    "in": [
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/input",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_input"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/output_file_name",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/family_size_histogram",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_family_size_histogram"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/raw_tag",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_raw_tag"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/assign_tag",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_assign_tag"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/include_non_pf_reads",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_include_non_pf_reads"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/strategy",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_strategy"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/edits",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_edits"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/min_umi_length",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_min_umi_length"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/temporary_directory",
                            "source": "#bam_collapsing.cwl/temporary_directory"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/async_io",
                            "source": "#bam_collapsing.cwl/async_io"
                        }
                    ],
                    "out": [
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/fgbio_group_reads_by_umi_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/fgbio_group_reads_by_umi_histogram"
                        }
                    ],
                    "run": "#fgbio_group_reads_by_umi_1.2.0.cwl",
                    "label": "fgbio_group_reads_by_umi_1.2.0",
                    "https://www.sevenbridges.com/x": 541.75,
                    "https://www.sevenbridges.com/y": 2865.0078125
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0",
                    "in": [
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/input",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/fgbio_group_reads_by_umi_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/output_prefix",
                            "source": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_output_prefix"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/intervals",
                            "source": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_intervals"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/description",
                            "source": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_description"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/duplex_umi_counts",
                            "source": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_duplex_umi_counts"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/min_ab_reads",
                            "source": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_min_ab_reads"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/min_ba_reads",
                            "source": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_min_ba_reads"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/mi_tag",
                            "source": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_mi_tag"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/temporary_directory",
                            "source": "#bam_collapsing.cwl/temporary_directory"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/async_io",
                            "source": "#bam_collapsing.cwl/async_io"
                        }
                    ],
                    "out": [
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_family_size"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_duplex_family_size"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_duplex_yield_metrics"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_umi_counts"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_duplex_qc"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_collect_duplex_seq_metrics_1_2_0/fgbio_collect_duplex_seq_metrics_duplex_umi_counts"
                        }
                    ],
                    "run": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl",
                    "label": "fgbio_collect_duplex_seq_metrics_1.2.0",
                    "https://www.sevenbridges.com/x": 1079.60107421875,
                    "https://www.sevenbridges.com/y": 2904.375
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0",
                    "in": [
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/input",
                            "source": "#bam_collapsing.cwl/fgbio_group_reads_by_umi_1_2_0/fgbio_group_reads_by_umi_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/output_file_name",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/read_name_prefix",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_read_name_prefix"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/read_group_id",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_read_group_id"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/error_rate_pre_umi",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_error_rate_pre_umi"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/error_rate_post_umi",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_error_rate_post_umi"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/min_input_base_quality",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_min_input_base_quality"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/trim",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_trim"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/sort_order",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_sort_order"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/min_reads",
                            "source": [
                                "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_min_reads"
                            ]
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/max_reads_per_strand",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_max_reads_per_strand"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/temporary_directory",
                            "source": "#bam_collapsing.cwl/temporary_directory"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/async_io",
                            "source": "#bam_collapsing.cwl/async_io"
                        }
                    ],
                    "out": [
                        {
                            "id": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/fgbio_call_duplex_consensus_reads_bam"
                        }
                    ],
                    "run": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl",
                    "label": "fgbio_call_duplex_consensus_reads_1.2.0",
                    "https://www.sevenbridges.com/x": 1079.60107421875,
                    "https://www.sevenbridges.com/y": 3158.109375
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0",
                    "in": [
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/fastq",
                            "source": "#bam_collapsing.cwl/gatk_sam_to_fastq_output_name_R1"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/input",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/fgbio_call_duplex_consensus_reads_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/include_non_pf_reads",
                            "source": "#bam_collapsing.cwl/gatk_sam_to_fastq_include_non_pf_reads"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/include_non_primary_alignments",
                            "source": "#bam_collapsing.cwl/gatk_sam_to_fastq_include_non_primary_alignments"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/reference_sequence",
                            "source": "#bam_collapsing.cwl/reference_sequence"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/second_end_fastq",
                            "source": "#bam_collapsing.cwl/gatk_sam_to_fastq_output_name_R2"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/unpaired_fastq",
                            "source": "#bam_collapsing.cwl/gatk_sam_to_fastq_output_name_unpaired"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/validation_stringency",
                            "source": "#bam_collapsing.cwl/validation_stringency"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/temporary_directory",
                            "source": "#bam_collapsing.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_fastq"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_unpaired_fastq"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_second_end_fastq"
                        }
                    ],
                    "run": "#gatk_sam_to_fastq_4.1.8.0.cwl",
                    "label": "GATK-SamToFastq",
                    "https://www.sevenbridges.com/x": 1680.5753173828125,
                    "https://www.sevenbridges.com/y": 2558.8203125
                },
                {
                    "id": "#bam_collapsing.cwl/alignment",
                    "in": [
                        {
                            "id": "#bam_collapsing.cwl/alignment/create_bam_index",
                            "source": "#bam_collapsing.cwl/create_bam_index"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/output_file_name",
                            "source": "#bam_collapsing.cwl/picard_addRG_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/read_group_description",
                            "source": "#bam_collapsing.cwl/picard_addRG_read_group_description"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/read_group_identifier",
                            "source": "#bam_collapsing.cwl/picard_addRG_read_group_identifier"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/read_group_library",
                            "source": "#bam_collapsing.cwl/picard_addRG_read_group_library"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/read_group_platform_unit",
                            "source": "#bam_collapsing.cwl/picard_addRG_read_group_platform_unit"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/read_group_run_date",
                            "source": "#bam_collapsing.cwl/picard_addRG_read_group_run_date"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/read_group_sample_name",
                            "source": "#bam_collapsing.cwl/picard_addRG_read_group_sample_name"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/read_group_sequencing_center",
                            "source": "#bam_collapsing.cwl/picard_addRG_read_group_sequencing_center"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/read_group_sequencing_platform",
                            "source": "#bam_collapsing.cwl/picard_addRG_read_group_sequencing_platform"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/sort_order",
                            "source": "#bam_collapsing.cwl/picard_addRG_sort_order"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/validation_stringency",
                            "source": "#bam_collapsing.cwl/validation_stringency"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/reference",
                            "source": "#bam_collapsing.cwl/reference_sequence"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/reads",
                            "source": [
                                "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_fastq",
                                "#bam_collapsing.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_second_end_fastq"
                            ]
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/output",
                            "source": "#bam_collapsing.cwl/bwa_mem_output"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/P",
                            "source": "#bam_collapsing.cwl/bwa_mem_P"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/M",
                            "source": "#bam_collapsing.cwl/bwa_mem_M"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/T",
                            "source": "#bam_collapsing.cwl/bwa_mem_T"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/Y",
                            "source": "#bam_collapsing.cwl/bwa_mem_Y"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/K",
                            "source": "#bam_collapsing.cwl/bwa_mem_K"
                        },
                        {
                            "id": "#bam_collapsing.cwl/alignment/bwa_number_of_threads",
                            "source": "#bam_collapsing.cwl/bwa_number_of_threads"
                        }
                    ],
                    "out": [
                        {
                            "id": "#bam_collapsing.cwl/alignment/picard_add_or_replace_read_groups_bam"
                        }
                    ],
                    "run": "#alignment.cwl",
                    "label": "alignment",
                    "https://www.sevenbridges.com/x": 2248.361328125,
                    "https://www.sevenbridges.com/y": 3148.4609375
                },
                {
                    "id": "#bam_collapsing.cwl/indel_realignment",
                    "in": [
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/window_size",
                            "source": "#bam_collapsing.cwl/abra2_window_size"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/soft_clip_contig",
                            "source": "#bam_collapsing.cwl/abra2_soft_clip_contig"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/scoring_gap_alignments",
                            "source": "#bam_collapsing.cwl/abra2_scoring_gap_alignments"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/reference_fasta",
                            "source": "#bam_collapsing.cwl/reference_sequence"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/no_sort",
                            "source": "#bam_collapsing.cwl/abra2_no_sort"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/maximum_mixmatch_rate",
                            "source": "#bam_collapsing.cwl/abra2_maximum_mixmatch_rate"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/maximum_average_depth",
                            "source": "#bam_collapsing.cwl/abra2_maximum_average_depth"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/input_bam",
                            "source": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0/gatk_merge_bam_alignment_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/contig_anchor",
                            "source": "#bam_collapsing.cwl/abra2_contig_anchor"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/consensus_sequence",
                            "source": "#bam_collapsing.cwl/abra2_consensus_sequence"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/bam_index",
                            "source": "#bam_collapsing.cwl/create_bam_index"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/option_bedgraph",
                            "source": "#bam_collapsing.cwl/bedtools_genomecov_option_bedgraph"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/no_edge_complex_indel",
                            "source": "#bam_collapsing.cwl/abra2_no_edge_complex_indel"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/distance_between_features",
                            "source": "#bam_collapsing.cwl/bedtools_merge_distance_between_features"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/output_bams",
                            "source": [
                                "#bam_collapsing.cwl/abra2_output_bams"
                            ]
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/validation_stringency",
                            "source": "#bam_collapsing.cwl/validation_stringency"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/sort_order",
                            "source": "#bam_collapsing.cwl/sort_order"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/output_file_name",
                            "source": "#bam_collapsing.cwl/picard_fixmate_information_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/create_bam_index",
                            "source": "#bam_collapsing.cwl/create_bam_index"
                        },
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/temporary_directory",
                            "source": "#bam_collapsing.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#bam_collapsing.cwl/indel_realignment/indel_realignment_bam"
                        }
                    ],
                    "run": "#indel_realignment.cwl",
                    "label": "indel_realignment",
                    "https://www.sevenbridges.com/x": 3467.848876953125,
                    "https://www.sevenbridges.com/y": 1627.9296875
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0",
                    "in": [
                        {
                            "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0/unmapped_bam",
                            "source": "#bam_collapsing.cwl/fgbio_call_duplex_consensus_reads_1_2_0/fgbio_call_duplex_consensus_reads_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0/reference",
                            "source": "#bam_collapsing.cwl/reference_sequence"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0/output_file_name",
                            "source": "#bam_collapsing.cwl/gatk_merge_bam_alignment_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0/aligned_bam",
                            "source": [
                                "#bam_collapsing.cwl/alignment/picard_add_or_replace_read_groups_bam"
                            ],
                            "valueFrom": "${ return [ self ]; }"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0/sort_order",
                            "source": "#bam_collapsing.cwl/sort_order"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0/validation_stringency",
                            "source": "#bam_collapsing.cwl/validation_stringency"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0/create_index",
                            "source": "#bam_collapsing.cwl/create_bam_index"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0/temporary_directory",
                            "source": "#bam_collapsing.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#bam_collapsing.cwl/gatk_merge_bam_alignment_4_1_8_0/gatk_merge_bam_alignment_bam"
                        }
                    ],
                    "run": "#gatk_merge_bam_alignment_4.1.8.0.cwl",
                    "label": "GATK-MergeBamAlignment",
                    "https://www.sevenbridges.com/x": 3001.30419921875,
                    "https://www.sevenbridges.com/y": 2245.6796875
                },
                {
                    "id": "#bam_collapsing.cwl/fgbio_separate_bams",
                    "in": [
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/reference_fasta",
                            "source": "#bam_collapsing.cwl/reference_sequence"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/input",
                            "source": "#bam_collapsing.cwl/indel_realignment/indel_realignment_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/reverse_per_base_tags_simplex_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/require_single_strand_agreement_simplex_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_require_single_strand_agreement_simplex_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/output_file_name_simplex_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_output_file_name_simplex_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/min_reads_simplex_duplex",
                            "source": [
                                "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_reads_simplex_duplex"
                            ]
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/min_mean_base_quality_simplex_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_mean_base_quality_simplex_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/max_base_error_rate_simplex_duplex",
                            "source": [
                                "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_base_error_rate_simplex_duplex"
                            ]
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/max_no_call_fraction_simplex_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_no_call_fraction_simplex_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/min_base_quality_simplex_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_base_quality_simplex_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/max_read_error_rate_simplex_duplex",
                            "source": [
                                "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_read_error_rate_simplex_duplex"
                            ]
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/reverse_per_base_tags_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_reverse_per_base_tags_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/require_single_strand_agreement_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_require_single_strand_agreement_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/output_file_name_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_output_file_name_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/min_reads_duplex",
                            "source": [
                                "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_reads_duplex"
                            ]
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/min_mean_base_quality_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_mean_base_quality_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/min_base_quality_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_base_quality_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/max_read_error_rate_duplex",
                            "source": [
                                "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_read_error_rate_duplex"
                            ]
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/max_no_call_fraction_duplex",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_no_call_fraction_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/max_base_error_rate_duplex",
                            "source": [
                                "#bam_collapsing.cwl/fgbio_filter_consensus_read_max_base_error_rate_duplex"
                            ]
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/validation_stringency",
                            "source": "#bam_collapsing.cwl/validation_stringency"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/output_file_name_duplex_aln_metrics",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/create_index",
                            "source": "#bam_collapsing.cwl/create_bam_index"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/output_file_name_simplex_aln_metrics",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/output_file_name_simpex",
                            "source": "#bam_collapsing.cwl/fgbio_postprocessing_output_file_name_simplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/min_simplex_reads",
                            "source": "#bam_collapsing.cwl/fgbio_filter_consensus_read_min_simplex_reads"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/temporary_directory",
                            "source": "#bam_collapsing.cwl/temporary_directory"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/async_io",
                            "source": "#bam_collapsing.cwl/async_io"
                        }
                    ],
                    "out": [
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/fgbio_filter_consensus_reads_duplex_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/fgbio_postprocessing_simplex_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/gatk_collect_alignment_summary_metrics_txt_duplex"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/fgbio_filter_consensus_reads_simplex_duplex_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/fgbio_separate_bams/gatk_collect_alignment_summary_metrics_txt_simplex"
                        }
                    ],
                    "run": "#fgbio_separate_bams.cwl",
                    "label": "fgbio_separate_bams",
                    "https://www.sevenbridges.com/x": 4160.3173828125,
                    "https://www.sevenbridges.com/y": 2935.0078125
                },
                {
                    "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_4_1_8_0",
                    "in": [
                        {
                            "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_4_1_8_0/input",
                            "source": "#bam_collapsing.cwl/indel_realignment/indel_realignment_bam"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_4_1_8_0/output_file_name",
                            "source": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_4_1_8_0/reference",
                            "source": "#bam_collapsing.cwl/reference_sequence"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_4_1_8_0/validation_stringency",
                            "source": "#bam_collapsing.cwl/validation_stringency"
                        },
                        {
                            "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_4_1_8_0/temporary_directory",
                            "source": "#bam_collapsing.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#bam_collapsing.cwl/gatk_collect_alignment_summary_metrics_4_1_8_0/gatk_collect_alignment_summary_metrics_txt"
                        }
                    ],
                    "run": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl",
                    "label": "GATK-CollectAlignmentSummaryMetrics",
                    "https://www.sevenbridges.com/x": 4160.3173828125,
                    "https://www.sevenbridges.com/y": 2611.2734375
                }
            ],
            "requirements": [
                {
                    "class": "SubworkflowFeatureRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Ronak Shah"
                }
            ],
            "https://schema.org/citation": "",
            "https://schema.org/codeRepository": "https://github.com/msk-access/standard_bam_processing",
            "https://schema.org/contributor": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9042-6213",
                    "https://schema.org/name": "Ronak Shah"
                }
            ],
            "https://schema.org/dateCreated": "2020-09-25",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "$namespaces": {
                "s": "https://schema.org/",
                "sbg": "https://www.sevenbridges.com/"
            }
        },
        {
            "class": "CommandLineTool",
            "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl",
            "baseCommand": [
                "fgbio"
            ],
            "inputs": [
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--input",
                        "shellQuote": false
                    },
                    "doc": "The input SAM or BAM file."
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output SAM or BAM file to write consensus reads."
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/read_name_prefix",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--read-name-prefix"
                    },
                    "doc": "The prefix all consensus read names"
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/read_group_id",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--read-group-id"
                    },
                    "doc": "The new read group ID for all the consensus reads."
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/error_rate_pre_umi",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--error-rate-pre-umi"
                    },
                    "doc": "The Phred-scaled error rate for an error prior to the UMIs being integrated."
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/error_rate_post_umi",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--error-rate-post-umi"
                    },
                    "doc": "The Phred-scaled error rate for an error post the UMIs have been integrated."
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/min_input_base_quality",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--min-input-base-quality"
                    },
                    "doc": "Ignore bases in raw reads that have Q below this value."
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/trim",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--trim"
                    },
                    "doc": "If true, quality trim input reads in addition to masking low Q bases"
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--sort-order"
                    },
                    "doc": "The sort order of the output, if :none: then the same as the input."
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/min_reads",
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--min-reads",
                        "itemSeparator": " ",
                        "shellQuote": false
                    },
                    "doc": "The minimum number of input reads to a consensus read."
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/max_reads_per_strand",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--max-reads-per-strand"
                    },
                    "doc": "The maximum number of reads to use when building a single-strand consensus. If more than this many reads are present in a tag family, the family is randomly downsampled to exactly max-reads reads."
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null."
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/async_io",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "separate": false,
                        "prefix": "--async-io="
                    },
                    "doc": "'Use asynchronous I/O where possible, e.g. for SAM and BAM files [=true|false].'"
                }
            ],
            "outputs": [
                {
                    "id": "#fgbio_call_duplex_consensus_reads_1.2.0.cwl/fgbio_call_duplex_consensus_reads_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n    return  inputs.input.basename.replace(/.bam/,'_cons.bam');\n}"
                    }
                }
            ],
            "doc": "Calls duplex consensus sequences from reads generated from the same double-stranded source molecule. Prior to running this tool, read must have been grouped with GroupReadsByUmi using the paired strategy. Doing so will apply (by default) MI tags to all reads of the form */A and */B where the /A and /B suffixes with the same identifier denote reads that are derived from opposite strands of the same source duplex molecule.\n\nReads from the same unique molecule are first partitioned by source strand and assembled into single strand consensus molecules as described by CallMolecularConsensusReads. Subsequently, for molecules that have at least one observation of each strand, duplex consensus reads are assembled by combining the evidence from the two single strand consensus reads.\n\nBecause of the nature of duplex sequencing, this tool does not support fragment reads - if found in the input they are ignored. Similarly, read pairs for which consensus reads cannot be generated for one or other read (R1 or R2) are omitted from the output.\n\nConsensus reads have a number of additional optional tags set in the resulting BAM file. The tag names follow a pattern where the first letter (a, b or c) denotes that the tag applies to the first single strand consensus (a), second single-strand consensus (b) or the final duplex consensus (c). The second letter is intended to capture the meaning of the tag (e.g. d=depth, m=min depth, e=errors/error-rate) and is upper case for values that are one per read and lower case for values that are one per base.",
            "label": "fgbio_call_duplex_consensus_reads_1.2.0",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx10G\"\n  }\n  else {\n      return \"-Xmx10G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "valueFrom": "-XX:-UseGCOverheadLimit"
                },
                {
                    "position": 1,
                    "valueFrom": "CallDuplexConsensusReads"
                },
                {
                    "position": 0,
                    "prefix": "--tmp-dir=",
                    "separate": false,
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 2,
                    "prefix": "--output",
                    "shellQuote": false,
                    "valueFrom": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n      return  inputs.input.basename.replace(/.bam/,'_cons.bam');\n}"
                },
                {
                    "position": 2,
                    "prefix": "--threads",
                    "valueFrom": "${\n    if(inputs.number_of_threads)\n        return inputs.number_of_threads\n    return runtime.cores\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 20000,
                    "coresMin": 16
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/fgbio:1.2.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "fgbio CallDuplexConsensusReads",
                    "http://usefulinc.com/ns/doap#revision": "1.2.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl",
            "baseCommand": [
                "fgbio"
            ],
            "inputs": [
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--input"
                    },
                    "doc": "Input BAM file generated by GroupReadByUmi."
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/output_prefix",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Prefix of output files to write."
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/intervals",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--intervals"
                    },
                    "doc": "Optional set of intervals over which to restrict analysis. [Optional]."
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/description",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--description"
                    },
                    "doc": "Description of data set used to label plots. Defaults to sample/library. [Optional]."
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/duplex_umi_counts",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--duplex-umi-counts"
                    },
                    "doc": "If true, produce the .duplex_umi_counts.txt file with counts of duplex UMI observations. [Optional]."
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/min_ab_reads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--min-ab-reads"
                    },
                    "doc": "Minimum AB reads to call a tag family a 'duplex'. [Optional]."
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/min_ba_reads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--min-ba-reads"
                    },
                    "doc": "Minimum BA reads to call a tag family a 'duplex'. [Optional]."
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/umi_tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--umi-tag"
                    },
                    "doc": "The tag containing the raw UMI. [Optional]."
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/mi_tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--mi-tag"
                    },
                    "doc": "The output tag for UMI grouping. [Optional]."
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null."
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/async_io",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "separate": false,
                        "prefix": "--async-io="
                    },
                    "doc": "'Use asynchronous I/O where possible, e.g. for SAM and BAM files [=true|false].'"
                }
            ],
            "outputs": [
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/fgbio_collect_duplex_seq_metrics_family_size",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n     if(inputs.output_prefix){\n         return  inputs.output_prefix + '.family_sizes.txt'\n     }\n     else{\n         return inputs.input.basename.replace('.bam','.family_sizes.txt')\n     }\n}"
                    }
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/fgbio_collect_duplex_seq_metrics_duplex_family_size",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_prefix){\n        return  inputs.output_prefix + '.duplex_family_sizes.txt'\n    }\n    else{\n        return inputs.input.basename.replace('.bam','.duplex_family_sizes.txt')\n    }\n}"
                    }
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/fgbio_collect_duplex_seq_metrics_duplex_yield_metrics",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_prefix){\n        return  inputs.output_prefix + '.duplex_yield_metrics.txt'\n    }\n    else{\n        return inputs.input.basename.replace('.bam','.duplex_yield_metrics.txt')\n    }\n}"
                    }
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/fgbio_collect_duplex_seq_metrics_umi_counts",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_prefix){\n        return  inputs.output_prefix + '.umi_counts.txt'\n    }\n    else{\n        return inputs.input.basename.replace('.bam','.umi_counts.txt')\n    }\n}"
                    }
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/fgbio_collect_duplex_seq_metrics_duplex_qc",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_prefix){\n        return  inputs.output_prefix + '.duplex_qc.pdf'\n    }\n    else{\n        return inputs.input.basename.replace('.bam','.duplex_qc.pdf')\n    }\n}"
                    }
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_1.2.0.cwl/fgbio_collect_duplex_seq_metrics_duplex_umi_counts",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "${\n    if(inputs.duplex_umi_counts){\n        if(inputs.output_prefix){\n            return  inputs.output_prefix + '.duplex_umi_counts.txt'\n        }\n        else{\n            return inputs.input.basename.replace('.bam','.duplex_umi_counts.txt')\n        }\n    }\n}"
                    }
                }
            ],
            "doc": "Collects a suite of metrics to QC duplex sequencing data.\nInputs ------\nThe input to this tool must be a BAM file that is either:\n1. The exact BAM output by the 'GroupReadsByUmi' tool (in the sort-order it was produced in) 2. A BAM file that has MI tags present on all reads (usually set by 'GroupReadsByUmi' and has been sorted with\n   'SortBam' into 'TemplateCoordinate' order.\n\nCalculation of metrics may be restricted to a set of regions using the '--intervals' parameter. This can significantly affect results as off-target reads in duplex sequencing experiments often have very different properties than on-target reads due to the lack of enrichment.\nSeveral metrics are calculated related to the fraction of tag families that have duplex coverage. The definition of \"duplex\" is controlled by the '--min-ab-reads' and '--min-ba-reads' parameters. The default is to treat any tag family with at least one observation of each strand as a duplex, but this could be made more stringent, e.g. by setting '--min-ab-reads=3 --min-ba-reads=3'. If different thresholds are used then '--min-ab-reads' must be the higher value.\nOutputs -------\nThe following output files are produced:\n1. <output>.family_sizes.txt: metrics on the frequency of different types of families of different sizes 2. <output>.duplex_family_sizes.txt: metrics on the frequency of duplex tag families by the number of observations\n   from each strand\n3. <output>.duplex_yield_metrics.txt: summary QC metrics produced using 5%, 10%, 15%...100% of the data 4. <output>.umi_counts.txt: metrics on the frequency of observations of UMIs within reads and tag families 5. <output>.duplex_qc.pdf: a series of plots generated from the preceding metrics files for visualization 6. <output>.duplex_umi_counts.txt: (optional) metrics on the frequency of observations of duplex UMIs within reads\n   and tag families. This file is only produced if the '--duplex-umi-counts' option is used as it requires significantly\n   more memory to track all pairs of UMIs seen when a large number of UMI sequences are present.\n\nWithin the metrics files the prefixes 'CS', 'SS' and 'DS' are used to mean:\n* CS: tag families where membership is defined solely on matching genome coordinates and strand * SS: single-stranded tag families where membership is defined by genome coordinates, strand and UMI; ie. 50/A and\n  50/B are considered different tag families.\n* DS: double-stranded tag families where membership is collapsed across single-stranded tag families from the same\n  double-stranded source molecule; i.e. 50/A and 50/B become one family\n\nRequirements ------------\nFor plots to be generated R must be installed and the ggplot2 package installed with suggested dependencies. Successfully executing the following in R will ensure a working installation:\ninstall.packages(\"ggplot2\", repos=\"http://cran.us.r-project.org\", dependencies=TRUE)",
            "label": "fgbio_collect_duplex_seq_metrics_1.2.0",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx12G\"\n  }\n  else {\n      return \"-Xmx12G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "valueFrom": "-XX:-UseGCOverheadLimit"
                },
                {
                    "position": 1,
                    "valueFrom": "CollectDuplexSeqMetrics"
                },
                {
                    "position": 0,
                    "prefix": "--tmp-dir=",
                    "separate": false,
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 2,
                    "prefix": "--output",
                    "valueFrom": "${\n    if(inputs.output_prefix){\n        return inputs.output_prefix\n    }\n    else{\n        return inputs.input.basename.replace(/.bam/,'')\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 16000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/fgbio:1.2.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "fgbio CollectDuplexSeqMetrics",
                    "http://usefulinc.com/ns/doap#revision": "1.2.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#fgbio_group_reads_by_umi_1.2.0.cwl",
            "baseCommand": [
                "fgbio"
            ],
            "inputs": [
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--input",
                        "shellQuote": false
                    },
                    "doc": "The input BAM file."
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "The output SAM or BAM file to be written."
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/family_size_histogram",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--family-size-histogram"
                    },
                    "doc": "Optional output of tag family size counts."
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/raw_tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--raw-tag"
                    },
                    "doc": "The tag containing the raw UMI."
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/assign_tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--assign-tag"
                    },
                    "doc": "The output tag for UMI grouping."
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/min_map_q",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--min-map-q"
                    },
                    "doc": "Minimum mapping quality."
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/include_non_pf_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--include-non-pf-reads"
                    }
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/strategy",
                    "type": "string",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--strategy"
                    },
                    "doc": "The UMI assignment strategy. (identity,edit,adjacency,paired)"
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/edits",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--edits"
                    },
                    "doc": "The allowable number of edits between UMIs."
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/min_umi_length",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--min-umi-length"
                    },
                    "doc": "The minimum UMI length. If not specified then all UMIs must have the same length, otherwise discard reads with UMIs shorter than this length and allow for differing UMI lengths."
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null."
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/async_io",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "separate": false,
                        "prefix": "--async-io="
                    },
                    "doc": "'Use asynchronous I/O where possible, e.g. for SAM and BAM files [=true|false].'"
                }
            ],
            "outputs": [
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/fgbio_group_reads_by_umi_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n    return  inputs.input.basename.replace(/.bam/,'_group.bam');\n}"
                    }
                },
                {
                    "id": "#fgbio_group_reads_by_umi_1.2.0.cwl/fgbio_group_reads_by_umi_histogram",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "${\n    if(inputs.family_size_histogram)\n        return inputs.family_size_histogram\n}"
                    }
                }
            ],
            "doc": "Groups reads together that appear to have come from the same original molecule. Reads are grouped by template, and then templates are sorted by the 5\u2019 mapping positions of the reads from the template, used from earliest mapping position to latest. Reads that have the same end positions are then sub-grouped by UMI sequence.\n\nAccepts reads in any order (including unsorted) and outputs reads sorted by:\n\nThe lower genome coordinate of the two outer ends of the templates\nThe sequencing library\nThe assigned UMI tag\nRead Name\nReads are aggressively filtered out so that only high quality reads/mappings are taken forward. Single-end reads must have mapping quality >= min-map-q. Paired-end reads must have both reads mapped to the same chromosome with both reads having mapping quality >= min-mapq. (Note: the MQ tag is required on reads with mapped mates).\n\nThis is done with the expectation that the next step is building consensus reads, where it is undesirable to either:\n\nAssign reads together that are really from different source molecules\nBuild two groups from reads that are really from the same molecule\nErrors in mapping reads could lead to both and therefore are minimized.\n\nGrouping of UMIs is performed by one of three strategies:\n\n1. identity: only reads with identical UMI sequences are grouped together. This strategy may be useful for evaluating data, but should generally be avoided as it will generate multiple UMI groups per original molecule in the presence of errors.\n2. edit: reads are clustered into groups such that each read within a group has at least one other read in the group with <= edits differences and there are inter-group pairings with <= edits differences. Effective when there are small numbers of reads per UMI, but breaks down at very high coverage of UMIs.\n3. adjacency: a version of the directed adjacency method described in umi_tools that allows for errors between UMIs but only when there is a count gradient.\n4. paired: similar to adjacency but for methods that produce template with a pair of UMIs such that a read with A-B is related to but not identical to a read with B-A. Expects the pair of UMIs to be stored in a single tag, separated by a hyphen (e.g. ACGT-CCGG). The molecular IDs produced have more structure than for single UMI strategies, and are of the form {base}/{AB|BA}. E.g. two UMI pairs would be mapped as follows AAAA-GGGG -> 1/AB, GGGG-AAAA -> 1/BA.\nedit, adjacency and paired make use of the --edits parameter to control the matching of non-identical UMIs.\n\nBy default, all UMIs must be the same length. If --min-umi-length=len is specified then reads that have a UMI shorter than len will be discarded, and when comparing UMIs of different lengths, the first len bases will be compared, where len is the length of the shortest UMI. The UMI length is the number of [ACGT] bases in the UMI (i.e. does not count dashes and other non-ACGT characters). This option is not implemented for reads with UMI pairs (i.e. using the paired assigner).",
            "label": "fgbio_group_reads_by_umi_1.2.0",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx12G\"\n  }\n  else {\n      return \"-Xmx12G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "valueFrom": "-XX:-UseGCOverheadLimit"
                },
                {
                    "position": 1,
                    "valueFrom": "GroupReadsByUmi"
                },
                {
                    "position": 0,
                    "prefix": "--tmp-dir=",
                    "separate": false,
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 2,
                    "prefix": "--output",
                    "shellQuote": false,
                    "valueFrom": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n      return  inputs.input.basename.replace(/.bam/,'_group.bam');\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 16000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/fgbio:1.2.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "fgbio GroupReadsByUmi",
                    "http://usefulinc.com/ns/doap#revision": "1.2.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl",
            "baseCommand": [
                "gatk",
                "CollectAlignmentSummaryMetrics"
            ],
            "inputs": [
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-I"
                    },
                    "doc": "Input file (bam or sam).  Required."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "File to write the output to.  Required."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/reference",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-R"
                    },
                    "doc": "Reference sequence file. Note that while this argument is not required, without it only a small subset of the metrics will be calculated. Note also that if a reference sequence is provided, it must be accompanied by a sequence dictionary.  Default value: null.",
                    "secondaryFiles": [
                        "^.fasta.fai",
                        "^.dict"
                    ]
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/adaptor_sequence",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ADAPTER_SEQUENCE"
                    },
                    "doc": "List of adapter sequences to use when processing the alignment metrics.  This argument may be specified 0 or more times. Default value: [AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG]."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/metrics_acciumulation_level",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--METRIC_ACCUMULATION_LEVEL"
                    },
                    "doc": "The level(s) at which to accumulate metrics. Default value: [ALL_READS]. This option can be set to 'null' to clear the default value. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} This option may be specified 0 or more times. This option can be set to 'null' to clear the default list."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/expected_pair_orientations",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--EXPECTED_PAIR_ORIENTATIONS"
                    },
                    "doc": "Paired-end reads that do not have this expected orientation will be considered chimeric. This argument may be specified 0 or more times. Default value: [FR]. Possible values: {FR, RF, TANDEM}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/is_bisulfite_sequenced",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--IS_BISULFITE_SEQUENCED"
                    },
                    "doc": "Whether the SAM or BAM file consists of bisulfite sequenced reads.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/max_insert_size",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MAX_INSERT_SIZE"
                    },
                    "doc": "Paired-end reads above this insert size will be considered chimeric along with inter-chromosomal pairs.  Default value: 100000."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT,LENIENT, SILENT}"
                },
                {
                    "default": true,
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/assume_sorted",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ASSUME_SORTED"
                    },
                    "doc": "If true (default), then the sort order in the header file will be ignored.  Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/stop_after",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--STOP_AFTER"
                    },
                    "doc": "Stop after processing N reads, mainly for debugging. Default value: 0. This option can be set to 'null' to clear the default value."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/create_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/create_md5_file",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_MD5_FILE"
                    },
                    "doc": "Whether to create an MD5 digest for any BAM or FASTQ files created.    Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl/gatk_collect_alignment_summary_metrics_txt",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if (inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/, '_alignment_summary_metrics.txt')\n    }\n}"
                    }
                }
            ],
            "label": "GATK-CollectAlignmentSummaryMetrics",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--java-options",
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx15G\"\n  }\n  else {\n      return \"-Xmx15G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/, '_alignment_summary_metrics.txt')\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 32000,
                    "coresMin": 1
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "gatk4",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl",
            "baseCommand": [
                "gatk",
                "MergeBamAlignment"
            ],
            "inputs": [
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/unmapped_bam",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--UNMAPPED_BAM"
                    },
                    "doc": "Original SAM or BAM file of unmapped reads, which must be in queryname order.  Reads MUST\nbe unmapped. Required.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/reference",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--REFERENCE_SEQUENCE"
                    },
                    "doc": "Reference sequence file.  Required.\n",
                    "secondaryFiles": [
                        "^.dict"
                    ]
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Merged SAM or BAM file to write to.  Required.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/add_mate_cigar",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ADD_MATE_CIGAR"
                    },
                    "doc": "Adds the mate CIGAR tag (MC) if true, does not if false.  Default value: true. Possible\nvalues: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/add_pg_tag_to_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ADD_PG_TAG_TO_READS"
                    },
                    "doc": "Add PG tag to each read in a SAM or BAM  Default value: true. Possible values: {true,\nfalse}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/aligned_bam",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File",
                            "inputBinding": {
                                "prefix": "--ALIGNED_BAM"
                            }
                        }
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "doc": "SAM or BAM file(s) with alignment data.  This argument may be specified 0 or more times.\nDefault value: null.  Cannot be used in conjunction with argument(s) READ1_ALIGNED_BAM\n(R1_ALIGNED) READ2_ALIGNED_BAM (R2_ALIGNED)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/aligned_reads_only",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ALIGNED_READS_ONLY"
                    },
                    "doc": "Whether to output only aligned reads. Default value: false. Possible values: {true,\nfalse}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/aligner_proper_pair_flags",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ALIGNER_PROPER_PAIR_FLAGS"
                    },
                    "doc": "Use the aligners idea of what a proper pair is rather than computing in this program.\nDefault value: false. Possible values: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/attributes_to_remove",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ATTRIBUTES_TO_REMOVE"
                    },
                    "doc": "Attributes from the alignment record that should be removed when merging.  This overrides\nATTRIBUTES_TO_RETAIN if they share common tags.  This argument may be specified 0 or more\ntimes. Default value: null.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/attributes_to_retain",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ATTRIBUTES_TO_RETAIN"
                    },
                    "doc": "Reserved alignment attributes (tags starting with X, Y, or Z) that should be brought over\nfrom the alignment data when merging.  This argument may be specified 0 or more times.\nDefault value: null.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/attributes_to_reverse",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ATTRIBUTES_TO_REVERSE"
                    },
                    "doc": "Attributes on negative strand reads that need to be reversed.  This argument may be\nspecified 0 or more times. Default value: [OQ, U2].\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/attributes_to_reverse_complement",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ATTRIBUTES_TO_REVERSE_COMPLEMENT"
                    },
                    "doc": "Attributes on negative strand reads that need to be reverse complemented.  This argument\nmay be specified 0 or more times. Default value: [E2, SQ].\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/clip_adapters",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CLIP_ADAPTERS"
                    },
                    "doc": "Whether to clip adapters where identified.  Default value: true. Possible values: {true,\nfalse}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/clip_overlapping_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CLIP_OVERLAPPING_READS"
                    },
                    "doc": "For paired reads, clip the 3' end of each read if necessary so that it does not extend\npast the 5' end of its mate.  Clipping will be either soft or hard clipping, depending on\nCLIP_OVERLAPPING_READS_OPERATOR setting. Hard clipped bases and their qualities will be\nstored in the XB and XQ tags respectively.  Default value: true. Possible values: {true,\nfalse}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/expected_orientations",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--EXPECTED_ORIENTATIONS"
                    },
                    "doc": "The expected orientation of proper read pairs. Replaces JUMP_SIZE  This argument may be\nspecified 0 or more times. Default value: null. Possible values: {FR, RF, TANDEM}  Cannot\nbe used in conjunction with argument(s) JUMP_SIZE (JUMP)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/hard_clip_overlapping_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--HARD_CLIP_OVERLAPPING_READS"
                    },
                    "doc": "If true, hard clipping will be applied to overlapping reads.  By default, soft clipping is\nused.  Default value: false. Possible values: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/include_secondary_alignments",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--INCLUDE_SECONDARY_ALIGNMENTS"
                    },
                    "doc": "If false, do not write secondary alignments to output.  Default value: true. Possible\nvalues: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/is_bisulfite_sequence",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--IS_BISULFITE_SEQUENCE"
                    },
                    "doc": "Whether the lane is bisulfite sequence (used when calculating the NM tag).  Default value:\nfalse. Possible values: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/jump_size",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--JUMP_SIZE"
                    },
                    "doc": "The expected jump size (required if this is a jumping library). Deprecated. Use\nEXPECTED_ORIENTATIONS instead  Default value: null.  Cannot be used in conjunction with\nargument(s) EXPECTED_ORIENTATIONS (ORIENTATIONS)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/matching_dictionary_tags",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MATCHING_DICTIONARY_TAGS"
                    },
                    "doc": "List of Sequence Records tags that must be equal (if present) in the reference dictionary\nand in the aligned file. Mismatching tags will cause an error if in this list, and a\nwarning otherwise.  This argument may be specified 0 or more times. Default value: [M5,\nLN].\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/max_insertions_or_deletions",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MAX_INSERTIONS_OR_DELETIONS"
                    },
                    "doc": "The maximum number of insertions or deletions permitted for an alignment to be included.\nAlignments with more than this many insertions or deletions will be ignored. Set to -1 to\nallow any number of insertions or deletions.  Default value: 1.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/min_unclipped_bases",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MIN_UNCLIPPED_BASES"
                    },
                    "doc": "If UNMAP_CONTAMINANT_READS is set, require this many unclipped bases or else the read will\nbe marked as contaminant.  Default value: 32.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/paired_run",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--PAIRED_RUN"
                    },
                    "doc": "DEPRECATED. This argument is ignored and will be removed.  Default value: true. Possible\nvalues: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/primary_alignment_strategy",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--PRIMARY_ALIGNMENT_STRATEGY"
                    },
                    "doc": "Strategy for selecting primary alignment when the aligner has provided more than one\nalignment for a pair or fragment, and none are marked as primary, more than one is marked\nas primary, or the primary alignment is filtered out for some reason. For all strategies,\nties are resolved arbitrarily.  Default value: BestMapq. BestMapq (Expects that multiple\nalignments will be correlated with HI tag, and prefers the pair of alignments with the\nlargest MAPQ, in the absence of a primary selected by the aligner.)\nEarliestFragment (Prefers the alignment which maps the earliest base in the read. Note\nthat EarliestFragment may not be used for paired reads.)\nBestEndMapq (Appropriate for cases in which the aligner is not pair-aware, and does not\noutput the HI tag. It simply picks the alignment for each end with the highest MAPQ, and\nmakes those alignments primary, regardless of whether the two alignments make sense\ntogether.)\nMostDistant (Appropriate for a non-pair-aware aligner. Picks the alignment pair with the\nlargest insert size. If all alignments would be chimeric, it picks the alignments for each\nend with the best MAPQ.)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/read1_aligned_bam",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File",
                            "inputBinding": {
                                "prefix": "--READ1_ALIGNED_BAM"
                            }
                        }
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "doc": "SAM or BAM file(s) with alignment data from the first read of a pair.  This argument may\nbe specified 0 or more times. Default value: null.  Cannot be used in conjunction with\nargument(s) ALIGNED_BAM (ALIGNED)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/read1_trim",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ1_TRIM"
                    },
                    "doc": "The number of bases trimmed from the beginning of read 1 prior to alignment  Default\nvalue: 0.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/read2_aligned_bam",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File",
                            "inputBinding": {
                                "prefix": "--READ2_ALIGNED_BAM"
                            }
                        }
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "doc": "SAM or BAM file(s) with alignment data from the second read of a pair.  This argument may\nbe specified 0 or more times. Default value: null.  Cannot be used in conjunction with\nargument(s) ALIGNED_BAM (ALIGNED)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/read2_trim",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ2_TRIM"
                    },
                    "doc": "The number of bases trimmed from the beginning of read 2 prior to alignment  Default\nvalue: 0.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--SORT_ORDER"
                    },
                    "doc": "The order in which the merged reads should be output.  Default value: coordinate. Possible\nvalues: {unsorted, queryname, coordinate, duplicate, unknown}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/unmap_contaminant_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--UNMAP_CONTAMINANT_READS"
                    },
                    "doc": "Detect reads originating from foreign organisms (e.g. bacterial DNA in a non-bacterial\nsample),and unmap + label those reads accordingly.  Default value: false. Possible values:\n{true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/unmapped_read_strategy",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--UNMAPPED_READ_STRATEGY"
                    },
                    "doc": "How to deal with alignment information in reads that are being unmapped (e.g. due to\ncross-species contamination.) Currently ignored unless UNMAP_CONTAMINANT_READS = true.\nNote that the DO_NOT_CHANGE strategy will actually reset the cigar and set the mapping\nquality on unmapped reads since otherwisethe result will be an invalid record. To force no\nchange use the DO_NOT_CHANGE_INVALID strategy.  Default value: DO_NOT_CHANGE. Possible\nvalues: {COPY_TO_TAG, DO_NOT_CHANGE, DO_NOT_CHANGE_INVALID, MOVE_TO_TAG}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT,LENIENT, SILENT}"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/create_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/create_md5_file",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_MD5_FILE"
                    },
                    "doc": "Whether to create an MD5 digest for any BAM or FASTQ files created.    Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl/gatk_merge_bam_alignment_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.unmapped_bam.basename.replace(/.bam|.sam/, '_merged.bam')\n    }\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "label": "GATK-MergeBamAlignment",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--java-options",
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx15G\"\n  }\n  else {\n      return \"-Xmx15G\"\n  }\n}"
                },
                {
                    "position": 1,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.unmapped_bam.basename.replace(/.bam|.sam/, '_merged.bam')\n    }\n}"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 17000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "gatk4",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#gatk_sam_to_fastq_4.1.8.0.cwl",
            "baseCommand": [
                "gatk",
                "SamToFastq"
            ],
            "inputs": [
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/fastq",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output FASTQ file (single-end fastq or, if paired, first end of the pair FASTQ)"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--INPUT"
                    },
                    "doc": "Input SAM/BAM file to extract reads from  Required."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/clipping_action",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CLIPPING_ACTION"
                    },
                    "doc": "The action that should be taken with clipped reads: 'X' means the reads and qualities should be trimmed at the clipped position; 'N' means the bases should be changed to Ns in the clipped region; and any integer means that the base qualities should be set to that value in the clipped region.  Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/clipping_attribute",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CLIPPING_ATTRIBUTE"
                    },
                    "doc": "The attribute that stores the position at which the SAM record should be clipped  Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/clipping_min_length",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CLIPPING_MIN_LENGTH"
                    },
                    "doc": "When performing clipping with the CLIPPING_ATTRIBUTE and CLIPPING_ACTION parameters, ensure that the resulting reads after clipping are at least CLIPPING_MIN_LENGTH bases long. If the original read is shorter than CLIPPING_MIN_LENGTH then the original read length will be maintained.  Default value: 0."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/compress_outputs_per_rg",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--COMPRESS_OUTPUTS_PER_RG"
                    },
                    "doc": "Compress output FASTQ files per read group using gzip and append a .gz extension to the file names.  Default value: false. Possible values: {true, false}  Cannot be used in conjunction with argument(s) FASTQ (F) SECOND_END_FASTQ (F2) UNPAIRED_FASTQ (FU)"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/compression_level",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--COMPRESSION_LEVEL"
                    },
                    "doc": "Compression level for all compressed files created (e.g. BAM and VCF).  Default value: 2."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/create_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/include_non_pf_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--INCLUDE_NON_PF_READS"
                    },
                    "doc": "Include non-PF reads from the SAM file into the output FASTQ files. PF means 'passes filtering'. Reads whose 'not passing quality controls' flag is set are non-PF reads. See GATK Dictionary for more info.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/include_non_primary_alignments",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--INCLUDE_NON_PRIMARY_ALIGNMENTS"
                    },
                    "doc": "If true, include non-primary alignments in the output.  Support of non-primary alignments in SamToFastq is not comprehensive, so there may be exceptions if this is set to true and there are paired reads with non-primary alignments.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/interleave",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--INTERLEAVE"
                    },
                    "doc": "Will generate an interleaved fastq if paired, each line will have /1 or /2 to describe which end it came from  Default value: false. Possible values: {true, false}"
                },
                {
                    "default": 50000,
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/max_records_in_ram",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MAX_RECORDS_IN_RAM"
                    },
                    "doc": "When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.  Default value: 500000."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/output_dir",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--OUTPUT_DIR"
                    },
                    "doc": "Directory in which to output the FASTQ file(s). Used only when OUTPUT_PER_RG is true. Default value: null. Cannot be used in conjunction with argument(s) FASTQ (F)."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/create_md5_file",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_MD5_FILE"
                    },
                    "doc": "Whether to create an MD5 digest for any BAM or FASTQ files created.    Default value: false. Possible values: {true, false}."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/output_per_rg",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--OUTPUT_PER_RG"
                    },
                    "doc": "Output a FASTQ file per read group (two FASTQ files per read group if the group is paired).  Default value: false. Possible values: {true, false}  Cannot be used in conjunction with argument(s) FASTQ (F) SECOND_END_FASTQ (F2) UNPAIRED_FASTQ (FU)"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/quality",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--QUALITY"
                    },
                    "doc": "End-trim reads using the phred/bwa quality trimming algorithm and this quality. Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/re_reverse",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RE_REVERSE"
                    },
                    "doc": "Re-reverse bases and qualities of reads with negative strand flag set before writing them to FASTQ  Default value: true. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/read1_max_bases_to_write",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ1_MAX_BASES_TO_WRITE"
                    },
                    "doc": "The maximum number of bases to write from read 1 after trimming. If there are fewer than this many bases left after trimming, all will be written.  If this value is null then all bases left after trimming will be written.  Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/read1_trim",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ1_TRIM"
                    },
                    "doc": "The number of bases to trim from the beginning of read 1.  Default value: 0."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/read2_max_bases_to_write",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ2_MAX_BASES_TO_WRITE"
                    },
                    "doc": "The maximum number of bases to write from read 2 after trimming. If there are fewer than this many bases left after trimming, all will be written.  If this value is null then all bases left after trimming will be written.  Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/read2_trim",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ2_TRIM"
                    },
                    "doc": "The number of bases to trim from the beginning of read 2.  Default value: 0."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/reference_sequence",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--REFERENCE_SEQUENCE"
                    },
                    "doc": "Reference sequence file. Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/rg_tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RG_TAG"
                    },
                    "doc": "The read group tag (PU or ID) to be used to output a FASTQ file per read group.  Default value: PU."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/second_end_fastq",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--SECOND_END_FASTQ"
                    },
                    "doc": "Output FASTQ file (if paired, second end of the pair FASTQ).  Default value: null.  Cannot be used in conjunction with argument(s) OUTPUT_PER_RG (OPRG) COMPRESS_OUTPUTS_PER_RG (GZOPRG)"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/unpaired_fastq",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--UNPAIRED_FASTQ"
                    },
                    "doc": "Output FASTQ file for unpaired reads; may only be provided in paired-FASTQ mode  Default value: null.  Cannot be used in conjunction with argument(s) OUTPUT_PER_RG (OPRG) COMPRESS_OUTPUTS_PER_RG (GZOPRG)"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT}"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/gatk_sam_to_fastq_fastq",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.fastq){\n      return inputs.fastq\n    } else {\n      return inputs.input.basename.replace(/.bam|.sam/, '_R1.fastq')\n    }\n}"
                    }
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/gatk_sam_to_fastq_unpaired_fastq",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "${\n    if(inputs.unpaired_fastq){\n        return inputs.unpaired_fastq\n    } else {\n      return inputs.input.basename.replace(/.bam|.sam/, '_unpaired.fastq')\n    }\n}"
                    }
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl/gatk_sam_to_fastq_second_end_fastq",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "${\n    if(inputs.second_end_fastq){\n        return inputs.second_end_fastq\n    } else {\n      return inputs.input.basename.replace(/.bam|.sam/, '_R2.fastq')\n    }\n}"
                    }
                }
            ],
            "label": "GATK-SamToFastq",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--java-options",
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx15G\"\n  }\n  else {\n      return \"-Xmx15G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "prefix": "--FASTQ",
                    "valueFrom": "${\n    if(inputs.fastq){\n        return inputs.fastq\n    } else {\n        return inputs.input.basename.replace(/.bam|.sam/, '_R1.fastq')\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 8000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "gatk4",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.0"
                }
            ]
        },
        {
            "class": "Workflow",
            "id": "#alignment.cwl",
            "label": "alignment",
            "inputs": [
                {
                    "id": "#alignment.cwl/create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 319.15625,
                    "https://www.sevenbridges.com/y": 958.8671875
                },
                {
                    "id": "#alignment.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 319.15625,
                    "https://www.sevenbridges.com/y": 852.0390625
                },
                {
                    "id": "#alignment.cwl/read_group_description",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1495.59375
                },
                {
                    "id": "#alignment.cwl/read_group_identifier",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1388.765625
                },
                {
                    "id": "#alignment.cwl/read_group_library",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1281.9375
                },
                {
                    "id": "#alignment.cwl/read_group_platform_unit",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1175.109375
                },
                {
                    "id": "#alignment.cwl/read_group_run_date",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1068.28125
                },
                {
                    "id": "#alignment.cwl/read_group_sample_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 961.453125
                },
                {
                    "id": "#alignment.cwl/read_group_sequencing_center",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 854.625
                },
                {
                    "id": "#alignment.cwl/read_group_sequencing_platform",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 747.796875
                },
                {
                    "id": "#alignment.cwl/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 427.3125
                },
                {
                    "id": "#alignment.cwl/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 106.828125
                },
                {
                    "id": "#alignment.cwl/reference",
                    "type": "File",
                    "secondaryFiles": [
                        ".amb",
                        ".fai",
                        ".sa",
                        "^.dict",
                        ".ann",
                        ".bwt",
                        ".pac"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 534.140625
                },
                {
                    "id": "#alignment.cwl/reads",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 640.96875
                },
                {
                    "id": "#alignment.cwl/output",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1709.25
                },
                {
                    "id": "#alignment.cwl/P",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1602.421875
                },
                {
                    "id": "#alignment.cwl/M",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1816.078125
                },
                {
                    "id": "#alignment.cwl/T",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 320.484375
                },
                {
                    "id": "#alignment.cwl/Y",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 0
                },
                {
                    "id": "#alignment.cwl/K",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1922.90625
                },
                {
                    "id": "#alignment.cwl/bwa_number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2029.734375
                },
                {
                    "id": "#alignment.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 213.65625
                }
            ],
            "outputs": [
                {
                    "id": "#alignment.cwl/picard_add_or_replace_read_groups_bam",
                    "outputSource": [
                        "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/picard_add_or_replace_read_groups_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 1389.239501953125,
                    "https://www.sevenbridges.com/y": 1014.8671875
                }
            ],
            "steps": [
                {
                    "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1",
                    "in": [
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/input",
                            "source": "#alignment.cwl/bwa_mem_0_7_17/bwa_mem_output_sam"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/output_file_name",
                            "source": "#alignment.cwl/output_file_name"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/sort_order",
                            "source": "#alignment.cwl/sort_order"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/read_group_identifier",
                            "source": "#alignment.cwl/read_group_identifier"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/read_group_sequencing_center",
                            "source": "#alignment.cwl/read_group_sequencing_center"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/read_group_library",
                            "source": "#alignment.cwl/read_group_library"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/read_group_platform_unit",
                            "source": "#alignment.cwl/read_group_platform_unit"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/read_group_sample_name",
                            "source": "#alignment.cwl/read_group_sample_name"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/read_group_sequencing_platform",
                            "source": "#alignment.cwl/read_group_sequencing_platform"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/read_group_description",
                            "source": "#alignment.cwl/read_group_description"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/read_group_run_date",
                            "source": "#alignment.cwl/read_group_run_date"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/validation_stringency",
                            "source": "#alignment.cwl/validation_stringency"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/create_bam_index",
                            "source": "#alignment.cwl/create_bam_index"
                        },
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/temporary_directory",
                            "source": "#alignment.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#alignment.cwl/picard_add_or_replace_read_groups_4_1_8_1/picard_add_or_replace_read_groups_bam"
                        }
                    ],
                    "run": "#picard_add_or_replace_read_groups_4.1.8.1.cwl",
                    "label": "picard_add_or_replace_read_groups_4.1.8.1",
                    "https://www.sevenbridges.com/x": 737.3328857421875,
                    "https://www.sevenbridges.com/y": 923.8671875
                },
                {
                    "id": "#alignment.cwl/bwa_mem_0_7_17",
                    "in": [
                        {
                            "id": "#alignment.cwl/bwa_mem_0_7_17/number_of_threads",
                            "source": "#alignment.cwl/bwa_number_of_threads"
                        },
                        {
                            "id": "#alignment.cwl/bwa_mem_0_7_17/reads",
                            "source": [
                                "#alignment.cwl/reads"
                            ]
                        },
                        {
                            "id": "#alignment.cwl/bwa_mem_0_7_17/reference",
                            "source": "#alignment.cwl/reference"
                        },
                        {
                            "id": "#alignment.cwl/bwa_mem_0_7_17/M",
                            "source": "#alignment.cwl/M"
                        },
                        {
                            "id": "#alignment.cwl/bwa_mem_0_7_17/P",
                            "source": "#alignment.cwl/P"
                        },
                        {
                            "id": "#alignment.cwl/bwa_mem_0_7_17/T",
                            "source": "#alignment.cwl/T"
                        },
                        {
                            "id": "#alignment.cwl/bwa_mem_0_7_17/K",
                            "source": "#alignment.cwl/K"
                        },
                        {
                            "id": "#alignment.cwl/bwa_mem_0_7_17/output",
                            "source": "#alignment.cwl/output"
                        },
                        {
                            "id": "#alignment.cwl/bwa_mem_0_7_17/Y",
                            "source": "#alignment.cwl/Y"
                        }
                    ],
                    "out": [
                        {
                            "id": "#alignment.cwl/bwa_mem_0_7_17/bwa_mem_output_sam"
                        }
                    ],
                    "run": "#bwa_mem_0.7.17.cwl",
                    "label": "bwa_mem_0.7.17",
                    "https://www.sevenbridges.com/x": 319.15625,
                    "https://www.sevenbridges.com/y": 1121.6953125
                }
            ],
            "requirements": [],
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:sumans@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Shalabh Suman"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:johnsoni@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Ian Jonhnson"
                }
            ],
            "https://schema.org/citation": "",
            "https://schema.org/codeRepository": "https://github.com/msk-access/cwl_subworkflows/alignment",
            "https://schema.org/contributor": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:johnsoni@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Ian Jonhnson"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9042-6213",
                    "https://schema.org/name": "Ronak Shah"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:sumans@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Shalabh Suman"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:murphyc4@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Charlie Murphy"
                }
            ],
            "https://schema.org/dateCreated": "2019-10-01",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0"
        },
        {
            "class": "CommandLineTool",
            "id": "#abra2_2.22.cwl",
            "baseCommand": [
                "java"
            ],
            "inputs": [
                {
                    "id": "#abra2_2.22.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#abra2_2.22.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#abra2_2.22.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#abra2_2.22.cwl/input_bam",
                    "type": [
                        "File",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--in"
                    },
                    "doc": "Required list of input sam or bam file (s) separated by comma",
                    "secondaryFiles": [
                        "^.bai"
                    ]
                },
                {
                    "id": "#abra2_2.22.cwl/working_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Set the temp directory (overrides java.io.tmpdir)"
                },
                {
                    "id": "#abra2_2.22.cwl/reference_fasta",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ref"
                    },
                    "doc": "Genome reference location",
                    "secondaryFiles": [
                        ".fai"
                    ]
                },
                {
                    "id": "#abra2_2.22.cwl/targets",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--targets"
                    }
                },
                {
                    "id": "#abra2_2.22.cwl/kmer_size",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--kmer"
                    },
                    "doc": "Optional assembly kmer size(delimit with commas if multiple sizes specified)"
                },
                {
                    "id": "#abra2_2.22.cwl/maximum_average_depth",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--mad"
                    },
                    "doc": "Regions with average depth exceeding this value will be downsampled (default: 1000)"
                },
                {
                    "id": "#abra2_2.22.cwl/soft_clip_contig",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--sc"
                    },
                    "doc": "Soft clip contig args [max_contigs,min_base_qual,frac_high_qual_bases,min_soft_clip_len] (default:16,13,80,15)"
                },
                {
                    "id": "#abra2_2.22.cwl/maximum_mixmatch_rate",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--mmr"
                    },
                    "doc": "Max allowed mismatch rate when mapping reads back to contigs (default: 0.05)"
                },
                {
                    "id": "#abra2_2.22.cwl/scoring_gap_alignments",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--sga"
                    },
                    "doc": "Scoring used for contig alignments(match, mismatch_penalty,gap_open_penalty,gap_extend_penalty) (default:8,32,48,1)"
                },
                {
                    "id": "#abra2_2.22.cwl/contig_anchor",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ca"
                    },
                    "doc": "Contig anchor [M_bases_at_contig_edge,max_mismatches_near_edge] (default:10,2)"
                },
                {
                    "id": "#abra2_2.22.cwl/window_size",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ws"
                    },
                    "doc": "Processing window size and overlap\n(size,overlap) (default: 400,200)"
                },
                {
                    "id": "#abra2_2.22.cwl/consensus_sequence",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--cons"
                    },
                    "doc": "Use positional consensus sequence when aligning high quality soft clipping"
                },
                {
                    "id": "#abra2_2.22.cwl/output_bams",
                    "type": [
                        "string",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--out"
                    },
                    "doc": "Required list of output sam or bam file (s) separated by comma"
                },
                {
                    "id": "#abra2_2.22.cwl/ignore_bad_assembly",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ignore-bad-assembly"
                    },
                    "doc": "Use this option to avoid parsing errors for corrupted assemblies"
                },
                {
                    "id": "#abra2_2.22.cwl/bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--index"
                    },
                    "doc": "Enable BAM index generation when outputting sorted alignments (may require additonal memory)"
                },
                {
                    "id": "#abra2_2.22.cwl/input_vcf",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--in-vcf"
                    },
                    "doc": "VCF containing known (or suspected) variant sites.  Very large files should be avoided."
                },
                {
                    "id": "#abra2_2.22.cwl/no_edge_complex_indel",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--no-edge-ci"
                    },
                    "doc": "Prevent output of complex indels at read start or read end"
                },
                {
                    "id": "#abra2_2.22.cwl/no_sort",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--nosort"
                    },
                    "doc": "Do not attempt to sort final output"
                }
            ],
            "outputs": [
                {
                    "id": "#abra2_2.22.cwl/abra_realigned_bam",
                    "type": [
                        "null",
                        "File",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "${\n    return inputs.output_bams\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "label": "abra2_2.22",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if (inputs.memory_per_job && inputs.memory_overhead) {\n\n    if (inputs.memory_per_job % 1000 == 0) {\n\n      return \"-Xmx\" + (inputs.memory_per_job / 1000).toString() + \"G\"\n    }\n    else {\n\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job / 1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead) {\n\n    if (inputs.memory_per_job % 1000 == 0) {\n\n      return \"-Xmx\" + (inputs.memory_per_job / 1000).toString() + \"G\"\n    }\n    else {\n\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job / 1000)).toString() + \"G\"\n    }\n  }\n  else if (!inputs.memory_per_job && inputs.memory_overhead) {\n\n    return \"-Xmx20G\"\n  }\n  else {\n\n    return \"-Xmx20G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "prefix": "-jar",
                    "valueFrom": "/usr/local/bin/abra2.jar"
                },
                {
                    "position": 0,
                    "prefix": "--threads",
                    "valueFrom": "${\n    if(inputs.number_of_threads)\n        return inputs.number_of_threads\n    return runtime.cores\n}"
                },
                {
                    "position": 0,
                    "prefix": "--tmpdir",
                    "valueFrom": "${\n    if(inputs.working_directory)\n        return inputs.working_directory;\n      return runtime.tmpdir\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 60000,
                    "coresMin": 16
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/abra2:2.22"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:kumarn1@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Nikhil Kumar"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "abra2",
                    "http://usefulinc.com/ns/doap#revision": 2.22
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#bedtools_genomecov_v2.28.0_cv2.cwl",
            "baseCommand": [
                "bedtools",
                "genomecov"
            ],
            "inputs": [
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 5,
                        "prefix": "-ibam",
                        "shellQuote": false
                    },
                    "doc": "The input file can be in BAM format (Note: BAM  must be sorted by position)",
                    "secondaryFiles": [
                        "^.bai"
                    ]
                },
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl/option_bedgraph",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-bg",
                        "separate": false
                    },
                    "doc": "option flag parameter to choose output file format. -bg refers to bedgraph format"
                }
            ],
            "outputs": [
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl/bedtools_genomecove_bedgraph",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n     if (inputs.output_file_name)\n      return inputs.output_file_name;\n    return inputs.input.basename.replace('.bam','.bedgraph');\n  }"
                    }
                }
            ],
            "label": "bedtools_genomecov",
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 20000,
                    "coresMin": 1
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/bedtools:v2.28.0_cv2"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "${\n    if (inputs.output_file_name)\n      return inputs.output_file_name;\n    return inputs.input.basename.replace('.bam','.bedgraph');\n  }",
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:sumans@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Shalabh Suman"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:sumans@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Shalabh Suman"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "bedtools",
                    "http://usefulinc.com/ns/doap#revision": "v2.28.0_cv2"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#bedtools_merge_v2.28.0_cv2.cwl",
            "baseCommand": [
                "bedtools",
                "merge"
            ],
            "inputs": [
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 5,
                        "prefix": "-i",
                        "shellQuote": false
                    },
                    "doc": "BEDgraph format file generated from Bedtools Genomecov module"
                },
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "default": 0,
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl/distance_between_features",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-d",
                        "shellQuote": false
                    },
                    "doc": "Maximum distance between features allowed for features to be merged."
                }
            ],
            "outputs": [
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl/bedtools_merge_bed",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if (inputs.output_file_name)\n      return inputs.output_file_name;\n    return inputs.input.basename.replace('.bedgraph', '.bed');\n  }"
                    }
                }
            ],
            "label": "bedtools_merge",
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 20000,
                    "coresMin": 1
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/bedtools:v2.28.0_cv2"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "${\n    if (inputs.output_file_name)\n      return inputs.output_file_name;\n    return inputs.input.basename.replace('.bedgraph', '.bed');\n  }",
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:sumans@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Shalabh Suman"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:sumans@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Shalabh Suman"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "bedtools",
                    "http://usefulinc.com/ns/doap#revision": "v2.28.0_cv2"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "baseCommand": [
                "bwa",
                "mem"
            ],
            "inputs": [
                {
                    "id": "#bwa_mem_0.7.17.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/reads",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 3
                    }
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/reference",
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    },
                    "secondaryFiles": [
                        ".amb",
                        ".ann",
                        ".bwt",
                        ".pac",
                        ".sa",
                        ".fai"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/A",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-A"
                    },
                    "doc": "score for a sequence match, which scales options -TdBOELU unless overridden [1]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/B",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-B"
                    },
                    "doc": "penalty for a mismatch [4]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/C",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-C"
                    },
                    "doc": "append FASTA/FASTQ comment to SAM output"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/E",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-E",
                        "itemSeparator": ","
                    },
                    "doc": "gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/L",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-L",
                        "itemSeparator": ","
                    },
                    "doc": "penalty for 5'- and 3'-end clipping [5,5]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/M",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-M"
                    }
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/O",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-O",
                        "itemSeparator": ","
                    },
                    "doc": "gap open penalties for deletions and insertions [6,6]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/P",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-P"
                    },
                    "doc": "skip pairing; mate rescue performed unless -S also in use"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/S",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-S"
                    },
                    "doc": "skip mate rescue"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/T",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-T"
                    },
                    "doc": "minimum score to output [30]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/U",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-U"
                    },
                    "doc": "penalty for an unpaired read pair [17]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/a",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-a"
                    },
                    "doc": "output all alignments for SE or unpaired PE"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/c",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-c"
                    },
                    "doc": "skip seeds with more than INT occurrences [500]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/d",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-d"
                    },
                    "doc": "off-diagonal X-dropoff [100]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/k",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-k"
                    },
                    "doc": "minimum seed length [19]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/K",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-K"
                    },
                    "doc": "process INT input bases in each batch regardless of nThreads (for reproducibility) []"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/output",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/p",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-p"
                    },
                    "doc": "smart pairing (ignoring in2.fq)"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/r",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-r"
                    },
                    "doc": "look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/v",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-v"
                    },
                    "doc": "verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/w",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-w"
                    },
                    "doc": "band width for banded alignment [100]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/y",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-y"
                    },
                    "doc": "seed occurrence for the 3rd round seeding [20]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/D",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-D"
                    },
                    "doc": "drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/W",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-W"
                    },
                    "doc": "discard a chain if seeded bases shorter than INT [0]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/m",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-m"
                    },
                    "doc": "perform at most INT rounds of mate rescues for each read [50]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/e",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-e"
                    }
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/x",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-x"
                    },
                    "doc": "read type. Setting -x changes multiple parameters unless overridden [null] pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref) ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref) intractg: -B9 -O16 -L5  (intra-species contigs to ref)"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/H",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-H"
                    },
                    "doc": "Use hard clipping \u2019H\u2019 in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/j",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-j"
                    },
                    "doc": "treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/he",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-h",
                        "itemSeparator": ","
                    },
                    "doc": "if there are <INT hits with score >80% of the max score, output all in XA [5,200]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/V",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-V"
                    },
                    "doc": "output the reference FASTA header in the XR tag"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/Y",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-Y"
                    },
                    "doc": "use soft clipping for supplementary alignments"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/I",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-M"
                    }
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/R",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "STR read group header line such as '@RG\\tID -foo\\tSM -bar' [null]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/sample_id",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/lane_id",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/platform",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/platform_unit",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/center_name",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl/library_id",
                    "type": [
                        "null",
                        "string"
                    ]
                }
            ],
            "outputs": [
                {
                    "id": "#bwa_mem_0.7.17.cwl/bwa_mem_output_sam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n  if (inputs.output)\n    return inputs.output;\n  return inputs.reads[0].basename.replace(/(fastq.gz)|(fq.gz)/, 'sam');\n}"
                    }
                }
            ],
            "doc": "bwa mem [-aCHMpP] [-t nThreads] [-k minSeedLen] [-w bandWidth] [-d zDropoff] [-r seedSplitRatio] [-c maxOcc] [-A matchScore] [-B mmPenalty] [-O gapOpenPen] [-E gapExtPen] [-L clipPen] [-U unpairPen] [-R RGline] [-v verboseLevel] db.prefix reads.fq [mates.fq]\nAlign 70bp-1Mbp query sequences with the BWA-MEM algorithm. Briefly, the algorithm works by seeding alignments with maximal exact matches (MEMs) and then extending seeds with the affine-gap Smith-Waterman algorithm (SW).\n\nIf mates.fq file is absent and option -p is not set, this command regards input reads are single-end. If mates.fq is present, this command assumes the i-th read in reads.fq and the i-th read in mates.fq constitute a read pair. If -p is used, the command assumes the 2i-th and the (2i+1)-th read in reads.fq constitute a read pair (such input file is said to be interleaved). In this case, mates.fq is ignored. In the paired-end mode, the mem command will infer the read orientation and the insert size distribution from a batch of reads.\n\nThe BWA-MEM algorithm performs local alignment. It may produce multiple primary alignments for different part of a query sequence. This is a crucial feature for long sequences. However, some tools such as Picard\u2019s markDuplicates does not work with split alignments. One may consider to use option -M to flag shorter split hits as secondary.",
            "label": "bwa_mem_0.7.17",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "-t",
                    "valueFrom": "$(runtime.cores)"
                },
                {
                    "position": 0,
                    "prefix": "-R",
                    "valueFrom": "${\n    if (inputs.sample_id) {\n        var rg_id = \"@RG\\\\tID:\" + inputs.sample_id + \"\\\\tSM:\" + inputs.sample_id;\n        if (inputs.library_id) {\n            rg_id += \"\\\\tLB:\" + inputs.library_id;\n        } if (inputs.platform) {\n            rg_id += \"\\\\tPL:\" + inputs.platform;\n        } if (inputs.platform_unit) {\n            rg_id += \"\\\\tPU:\" + inputs.platform_unit;\n        } if (inputs.center_name) {\n            rg_id += \"\\\\tCN:\" + inputs.center_name;\n        }\n        return rg_id\n    } else {\n        return inputs.R\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 34000,
                    "coresMin": 16
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/bwa:0.7.17"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "${\n  if (inputs.output)\n    return inputs.output;\n  return inputs.reads[0].basename.replace(/(fastq.gz)|(fq.gz)/, 'sam');\n}",
            "id": "#bwa_mem_0.7.17.cwl",
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:johnsoni@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ian Johnson"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "bwa",
                    "http://usefulinc.com/ns/doap#revision": "0.7.17"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#fgbio_filter_consensus_reads_1.2.0.cwl",
            "baseCommand": [
                "fgbio"
            ],
            "inputs": [
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--input",
                        "shellQuote": false
                    },
                    "doc": "The input SAM or BAM file."
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output SAM or BAM file to write consensus reads."
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/reference_fasta",
                    "type": "File",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--ref"
                    },
                    "doc": "Reference fasta file.",
                    "secondaryFiles": [
                        ".fai",
                        "^.dict"
                    ]
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/reverse_per_base_tags",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--reverse-per-base-tags"
                    },
                    "doc": "Reverse [complement] per base tags on reverse strand reads."
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/min_reads",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--min-reads",
                        "itemSeparator": " ",
                        "shellQuote": false
                    },
                    "doc": "The minimum number of reads supporting a consensus base/read. (Max 3 values)"
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/max_read_error_rate",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "float"
                        }
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--max-read-error-rate",
                        "itemSeparator": " "
                    },
                    "doc": "The maximum raw-read error rate across the entire consensus read. (Max 3 values)"
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/max_base_error_rate",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "float"
                        }
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--max-base-error-rate",
                        "itemSeparator": " "
                    },
                    "doc": "The maximum error rate for a single consensus base. (Max 3 values)"
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/min_base_quality",
                    "type": "int",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--min-base-quality"
                    },
                    "doc": "Mask (make N) consensus bases with quality less than this threshold."
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/max_no_call_fraction",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--max-no-call-fraction"
                    },
                    "doc": "Maximum fraction of no-calls in the read after filtering"
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/min_mean_base_quality",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--min-mean-base-quality"
                    },
                    "doc": "The minimum mean base quality across the consensus read"
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/require_single_strand_agreement",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--require-single-strand-agreement"
                    },
                    "doc": "Mask (make N) consensus bases where the AB and BA consensus reads disagree (for duplex-sequencing only)."
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null."
                },
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/async_io",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "separate": false,
                        "prefix": "--async-io="
                    },
                    "doc": "'Use asynchronous I/O where possible, e.g. for SAM and BAM files [=true|false].'"
                }
            ],
            "outputs": [
                {
                    "id": "#fgbio_filter_consensus_reads_1.2.0.cwl/fgbio_filter_consensus_reads_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n    return  inputs.input.basename.replace(/.bam/,'_filtered.bam');\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "doc": "Filters consensus reads generated by CallMolecularConsensusReads or CallDuplexConsensusReads. Two kinds of filtering are performed:\n\n1. Masking/filtering of individual bases in reads\n2. Filtering out of reads (i.e. not writing them to the output file)\n\nBase-level filtering/masking is only applied if per-base tags are present (see CallDuplexConsensusReads and CallMolecularConsensusReads for descriptions of these tags). Read-level filtering is always applied. When filtering reads, secondary alignments and supplementary records may be removed independently if they fail one or more filters; if either R1 or R2 primary alignments fail a filter then all records for the template will be filtered out.\n\nThe filters applied are as follows:\n\n1. Reads with fewer than min-reads contributing reads are filtered out\n2. Reads with an average consensus error rate higher than max-read-error-rate are filtered out\n3. Reads with mean base quality of the consensus read, prior to any masking, less than min-mean-base-quality are filtered out (if specified)\n4. Bases with quality scores below min-base-quality are masked to Ns\n5. Bases with fewer than min-reads contributing raw reads are masked to Ns\n6. Bases with a consensus error rate (defined as the fraction of contributing reads that voted for a different base than the consensus call) higher than max-base-error-rate are masked to Ns\n7. For duplex reads, if require-single-strand-agreement is provided, masks to Ns any bases where the base was observed in both single-strand consensus reads and the two reads did not agree\n8. Reads with a proportion of Ns higher than max-no-call-fraction after per-base filtering are filtered out",
            "label": "fgbio_filter_consensus_reads_1.2.0",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx12G\"\n  }\n  else {\n      return \"-Xmx12G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "valueFrom": "-XX:-UseGCOverheadLimit"
                },
                {
                    "position": 1,
                    "valueFrom": "FilterConsensusReads"
                },
                {
                    "position": 0,
                    "prefix": "--tmp-dir=",
                    "separate": false,
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 2,
                    "prefix": "--output",
                    "shellQuote": false,
                    "valueFrom": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n      return  inputs.input.basename.replace(/.bam/,'_filtered.bam');\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 16000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/fgbio:1.2.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "fgbio FilterConsensusReads",
                    "http://usefulinc.com/ns/doap#revision": "1.2.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#fgbio_postprocessing_simplex_filter_0.1.8.cwl",
            "baseCommand": [
                "simplex_filter"
            ],
            "inputs": [
                {
                    "id": "#fgbio_postprocessing_simplex_filter_0.1.8.cwl/input_bam",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--input_bam"
                    },
                    "doc": "Input file (bam or sam).  Required.",
                    "secondaryFiles": [
                        "^.bai"
                    ]
                },
                {
                    "id": "#fgbio_postprocessing_simplex_filter_0.1.8.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--output_filename"
                    },
                    "doc": "Output file (bam or sam)."
                },
                {
                    "id": "#fgbio_postprocessing_simplex_filter_0.1.8.cwl/min_simplex_reads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--min_simplex_reads"
                    },
                    "doc": "Minimum number of simplex reads to pass filter for consensus reads"
                }
            ],
            "outputs": [
                {
                    "id": "#fgbio_postprocessing_simplex_filter_0.1.8.cwl/fgbio_postprocessing_simplex_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if (inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input_bam.basename.replace(/.bam$/,'_simplex.bam')\n    }\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "label": "fgbio_postprocessing_simplex_filter_0.1.8",
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 16000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/fgbio_postprocessing:0.2.1"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:johnsoni@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ian Johnson"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:johnsoni@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ian Johnson"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "fgbio_postprocessing",
                    "http://usefulinc.com/ns/doap#revision": "0.1.8"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2",
            "baseCommand": [
                "gatk",
                "CollectAlignmentSummaryMetrics"
            ],
            "inputs": [
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-I"
                    },
                    "doc": "Input file (bam or sam).  Required."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "File to write the output to.  Required."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/reference",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-R"
                    },
                    "doc": "Reference sequence file. Note that while this argument is not required, without it only a small subset of the metrics will be calculated. Note also that if a reference sequence is provided, it must be accompanied by a sequence dictionary.  Default value: null.",
                    "secondaryFiles": [
                        "^.fasta.fai",
                        "^.dict"
                    ]
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/adaptor_sequence",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ADAPTER_SEQUENCE"
                    },
                    "doc": "List of adapter sequences to use when processing the alignment metrics.  This argument may be specified 0 or more times. Default value: [AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG]."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/metrics_acciumulation_level",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--METRIC_ACCUMULATION_LEVEL"
                    },
                    "doc": "The level(s) at which to accumulate metrics. Default value: [ALL_READS]. This option can be set to 'null' to clear the default value. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} This option may be specified 0 or more times. This option can be set to 'null' to clear the default list."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/expected_pair_orientations",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--EXPECTED_PAIR_ORIENTATIONS"
                    },
                    "doc": "Paired-end reads that do not have this expected orientation will be considered chimeric. This argument may be specified 0 or more times. Default value: [FR]. Possible values: {FR, RF, TANDEM}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/is_bisulfite_sequenced",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--IS_BISULFITE_SEQUENCED"
                    },
                    "doc": "Whether the SAM or BAM file consists of bisulfite sequenced reads.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/max_insert_size",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MAX_INSERT_SIZE"
                    },
                    "doc": "Paired-end reads above this insert size will be considered chimeric along with inter-chromosomal pairs.  Default value: 100000."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT,LENIENT, SILENT}"
                },
                {
                    "default": true,
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/assume_sorted",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ASSUME_SORTED"
                    },
                    "doc": "If true (default), then the sort order in the header file will be ignored.  Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/stop_after",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--STOP_AFTER"
                    },
                    "doc": "Stop after processing N reads, mainly for debugging. Default value: 0. This option can be set to 'null' to clear the default value."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/create_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/create_md5_file",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_MD5_FILE"
                    },
                    "doc": "Whether to create an MD5 digest for any BAM or FASTQ files created.    Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2/gatk_collect_alignment_summary_metrics_txt",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if (inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/, '_alignment_summary_metrics.txt')\n    }\n}"
                    }
                }
            ],
            "label": "GATK-CollectAlignmentSummaryMetrics",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--java-options",
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx15G\"\n  }\n  else {\n      return \"-Xmx15G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/, '_alignment_summary_metrics.txt')\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 32000,
                    "coresMin": 1
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "gatk4",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl",
            "baseCommand": [
                "java"
            ],
            "inputs": [
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-I"
                    },
                    "doc": "Input file ( sam).  Required."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output file name (bam or sam). Not Required"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-SO"
                    },
                    "doc": "Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.Default value: null. Possible values: {unsorted, queryname, coordinate}"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/read_group_identifier",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGID"
                    },
                    "doc": "Read Group ID  Default value: 1. This option can be set to 'null' to clear the default value  Required"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/read_group_sequencing_center",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGCN"
                    },
                    "doc": "Read Group sequencing center name  Default value: null. Required"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/read_group_library",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGLB"
                    },
                    "doc": "Read Group Library.  Required"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/read_group_platform_unit",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGPU"
                    },
                    "doc": "Read Group platform unit (eg. run barcode)  Required."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/read_group_sample_name",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGSM"
                    },
                    "doc": "Read Group sample name.  Required"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/read_group_sequencing_platform",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGPL"
                    },
                    "doc": "Read Group platform (e.g. illumina, solid)  Required."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/read_group_description",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGDS"
                    },
                    "doc": "Read Group description  Default value: null."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/read_group_run_date",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGDT"
                    },
                    "doc": "Read Group run date  Default value: null."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT,LENIENT, SILENT}"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/bam_compression_level",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--COMPRESSION_LEVEL"
                    },
                    "doc": "Compression level for all compressed files created (e.g. BAM and GELI). Default value:5. This option can be set to 'null' to clear the default value."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "default": true,
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value:false. This option can be set to 'null' to clear the default value. Possible values:{true, false}"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl/picard_add_or_replace_read_groups_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n    return inputs.input.basename.replace(/.sam$/, '_srt.bam');\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "label": "picard_add_or_replace_read_groups_4.1.8.1",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx15G\"\n  }\n  else {\n      return \"-Xmx15G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "prefix": "-Djava.io.tmpdir=",
                    "separate": false,
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "shellQuote": false,
                    "valueFrom": "-XX:-UseGCOverheadLimit"
                },
                {
                    "position": 0,
                    "prefix": "-jar",
                    "valueFrom": "/gatk/gatk-package-4.1.8.1-local.jar"
                },
                {
                    "position": 0,
                    "valueFrom": "AddOrReplaceReadGroups"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n      return inputs.input.basename.replace(/.sam$/, '_srt.bam');\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 17000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.1"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:kumarn1@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Nikhil Kumar"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "picard",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.1"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#picard_fix_mate_information_4.1.8.1.cwl",
            "baseCommand": [
                "java"
            ],
            "inputs": [
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-I"
                    },
                    "doc": "The input file to fix.  This option may be specified 0 or more times"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output file name (bam or sam). Not Required"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-SO"
                    },
                    "doc": "Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.Default value: null. Possible values: {unsorted, queryname, coordinate}"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT,LENIENT, SILENT}"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/bam_compression_level",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--COMPRESSION_LEVEL"
                    },
                    "doc": "Compression level for all compressed files created (e.g. BAM and GELI). Default value:5. This option can be set to 'null' to clear the default value."
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "default": true,
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value:false. This option can be set to 'null' to clear the default value. Possible values:{true, false}"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl/picard_fix_mate_information_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/,'_fm.bam')\n    }\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "label": "picard_fix_mate_information_4.1.8.1",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx20G\"\n  }\n  else {\n      return \"-Xmx20G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "shellQuote": false,
                    "valueFrom": "-XX:-UseGCOverheadLimit"
                },
                {
                    "position": 0,
                    "prefix": "-jar",
                    "valueFrom": "/gatk/gatk-package-4.1.8.1-local.jar"
                },
                {
                    "position": 0,
                    "valueFrom": "FixMateInformation"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/,'_fm.bam')\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 30000,
                    "coresMin": 12
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.1"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:kumarn1@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Nikhil Kumar"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "picard",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.1"
                }
            ]
        },
        {
            "class": "Workflow",
            "id": "#fgbio_separate_bams.cwl",
            "label": "fgbio_separate_bams",
            "inputs": [
                {
                    "id": "#fgbio_separate_bams.cwl/reference_fasta",
                    "type": "File",
                    "secondaryFiles": [
                        ".fai",
                        "^.dict"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 853.8671875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/input",
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3201.7734375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/reverse_per_base_tags_simplex_duplex",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 426.9375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/require_single_strand_agreement_simplex_duplex",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 640.40625
                },
                {
                    "id": "#fgbio_separate_bams.cwl/output_file_name_simplex_duplex",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 960.5859375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1494.1796875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/min_reads_simplex_duplex",
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1707.6171875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/min_mean_base_quality_simplex_duplex",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1921.0625
                },
                {
                    "id": "#fgbio_separate_bams.cwl/max_base_error_rate_simplex_duplex",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "float"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2988.3359375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/max_no_call_fraction_simplex_duplex",
                    "type": [
                        "null",
                        "float"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2774.8984375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/min_base_quality_simplex_duplex",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2134.53125
                },
                {
                    "id": "#fgbio_separate_bams.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2348
                },
                {
                    "id": "#fgbio_separate_bams.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2454.734375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/max_read_error_rate_simplex_duplex",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "float"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2561.4609375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/reverse_per_base_tags_duplex",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 533.671875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/require_single_strand_agreement_duplex",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 747.140625
                },
                {
                    "id": "#fgbio_separate_bams.cwl/output_file_name_duplex",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1387.4609375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/min_reads_duplex",
                    "type": {
                        "type": "array",
                        "items": "int"
                    },
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1814.3359375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/min_mean_base_quality_duplex",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2027.796875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/min_base_quality_duplex",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2241.265625
                },
                {
                    "id": "#fgbio_separate_bams.cwl/max_read_error_rate_duplex",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "float"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2668.1796875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/max_no_call_fraction_duplex",
                    "type": [
                        "null",
                        "float"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2881.6171875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/max_base_error_rate_duplex",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "float"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3095.0546875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 0
                },
                {
                    "id": "#fgbio_separate_bams.cwl/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 106.734375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 213.46875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/output_file_name_duplex_aln_metrics",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1280.7421875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/create_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 454.671875,
                    "https://www.sevenbridges.com/y": 1805.625
                },
                {
                    "id": "#fgbio_separate_bams.cwl/assume_sorted",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 454.671875,
                    "https://www.sevenbridges.com/y": 1912.34375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/output_file_name_simplex_aln_metrics",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1067.3046875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/output_file_name_simpex",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1174.0234375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/min_simplex_reads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1600.8984375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 320.203125
                },
                {
                    "id": "#fgbio_separate_bams.cwl/async_io",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3308.5
                }
            ],
            "outputs": [
                {
                    "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_duplex_bam",
                    "outputSource": [
                        "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/fgbio_filter_consensus_reads_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 1072.9705810546875,
                    "https://www.sevenbridges.com/y": 1828.3515625
                },
                {
                    "id": "#fgbio_separate_bams.cwl/fgbio_postprocessing_simplex_bam",
                    "outputSource": [
                        "#fgbio_separate_bams.cwl/fgbio_postprocessing_simplex_filter_0_1_8/fgbio_postprocessing_simplex_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 1616.9268798828125,
                    "https://www.sevenbridges.com/y": 1809.984375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_txt_duplex",
                    "outputSource": [
                        "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/gatk_collect_alignment_summary_metrics_txt"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1616.9268798828125,
                    "https://www.sevenbridges.com/y": 1498.515625
                },
                {
                    "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_simplex_duplex_bam",
                    "outputSource": [
                        "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/fgbio_filter_consensus_reads_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 1072.9705810546875,
                    "https://www.sevenbridges.com/y": 1721.6171875
                },
                {
                    "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_txt_simplex",
                    "outputSource": [
                        "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex/gatk_collect_alignment_summary_metrics_txt"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 2134.5888671875,
                    "https://www.sevenbridges.com/y": 1654.25
                }
            ],
            "steps": [
                {
                    "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex",
                    "in": [
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/memory_overhead",
                            "source": "#fgbio_separate_bams.cwl/memory_overhead"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/number_of_threads",
                            "source": "#fgbio_separate_bams.cwl/number_of_threads"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/input",
                            "source": "#fgbio_separate_bams.cwl/input"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/output_file_name",
                            "source": "#fgbio_separate_bams.cwl/output_file_name_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/reference_fasta",
                            "source": "#fgbio_separate_bams.cwl/reference_fasta"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/reverse_per_base_tags",
                            "source": "#fgbio_separate_bams.cwl/reverse_per_base_tags_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/min_reads",
                            "source": [
                                "#fgbio_separate_bams.cwl/min_reads_duplex"
                            ]
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/max_read_error_rate",
                            "source": [
                                "#fgbio_separate_bams.cwl/max_read_error_rate_duplex"
                            ]
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/max_base_error_rate",
                            "source": [
                                "#fgbio_separate_bams.cwl/max_base_error_rate_duplex"
                            ]
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/min_base_quality",
                            "source": "#fgbio_separate_bams.cwl/min_base_quality_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/max_no_call_fraction",
                            "source": "#fgbio_separate_bams.cwl/max_no_call_fraction_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/min_mean_base_quality",
                            "source": "#fgbio_separate_bams.cwl/min_mean_base_quality_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/require_single_strand_agreement",
                            "source": "#fgbio_separate_bams.cwl/require_single_strand_agreement_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/temporary_directory",
                            "source": "#fgbio_separate_bams.cwl/temporary_directory"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/async_io",
                            "source": "#fgbio_separate_bams.cwl/async_io"
                        }
                    ],
                    "out": [
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/fgbio_filter_consensus_reads_bam"
                        }
                    ],
                    "run": "#fgbio_filter_consensus_reads_1.2.0.cwl",
                    "label": "fgbio_filter_consensus_reads_1.2.0_duplex",
                    "https://www.sevenbridges.com/x": 454.671875,
                    "https://www.sevenbridges.com/y": 1600.8984375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex",
                    "in": [
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/memory_per_job",
                            "source": "#fgbio_separate_bams.cwl/memory_per_job"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/memory_overhead",
                            "source": "#fgbio_separate_bams.cwl/memory_overhead"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/number_of_threads",
                            "source": "#fgbio_separate_bams.cwl/number_of_threads"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/input",
                            "source": "#fgbio_separate_bams.cwl/input"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/output_file_name",
                            "source": "#fgbio_separate_bams.cwl/output_file_name_simplex_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/reference_fasta",
                            "source": "#fgbio_separate_bams.cwl/reference_fasta"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/reverse_per_base_tags",
                            "source": "#fgbio_separate_bams.cwl/reverse_per_base_tags_simplex_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/min_reads",
                            "source": [
                                "#fgbio_separate_bams.cwl/min_reads_simplex_duplex"
                            ]
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/max_read_error_rate",
                            "source": [
                                "#fgbio_separate_bams.cwl/max_read_error_rate_simplex_duplex"
                            ]
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/max_base_error_rate",
                            "source": [
                                "#fgbio_separate_bams.cwl/max_base_error_rate_simplex_duplex"
                            ]
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/min_base_quality",
                            "source": "#fgbio_separate_bams.cwl/min_base_quality_simplex_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/max_no_call_fraction",
                            "source": "#fgbio_separate_bams.cwl/max_no_call_fraction_simplex_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/min_mean_base_quality",
                            "source": "#fgbio_separate_bams.cwl/min_mean_base_quality_simplex_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/require_single_strand_agreement",
                            "source": "#fgbio_separate_bams.cwl/require_single_strand_agreement_simplex_duplex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/temporary_directory",
                            "source": "#fgbio_separate_bams.cwl/temporary_directory"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/async_io",
                            "source": "#fgbio_separate_bams.cwl/async_io"
                        }
                    ],
                    "out": [
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/fgbio_filter_consensus_reads_bam"
                        }
                    ],
                    "run": "#fgbio_filter_consensus_reads_1.2.0.cwl",
                    "label": "fgbio_filter_consensus_reads_1.2.0_simplex_duplex",
                    "https://www.sevenbridges.com/x": 454.671875,
                    "https://www.sevenbridges.com/y": 1291.1640625
                },
                {
                    "id": "#fgbio_separate_bams.cwl/fgbio_postprocessing_simplex_filter_0_1_8",
                    "in": [
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_postprocessing_simplex_filter_0_1_8/input_bam",
                            "source": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_1_simplex_duplex/fgbio_filter_consensus_reads_bam"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_postprocessing_simplex_filter_0_1_8/output_file_name",
                            "source": "#fgbio_separate_bams.cwl/output_file_name_simpex"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_postprocessing_simplex_filter_0_1_8/min_simplex_reads",
                            "source": "#fgbio_separate_bams.cwl/min_simplex_reads"
                        }
                    ],
                    "out": [
                        {
                            "id": "#fgbio_separate_bams.cwl/fgbio_postprocessing_simplex_filter_0_1_8/fgbio_postprocessing_simplex_bam"
                        }
                    ],
                    "run": "#fgbio_postprocessing_simplex_filter_0.1.8.cwl",
                    "label": "fgbio_postprocessing_simplex_filter_0.1.8",
                    "https://www.sevenbridges.com/x": 1072.9705810546875,
                    "https://www.sevenbridges.com/y": 1600.8828125
                },
                {
                    "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex",
                    "in": [
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/input",
                            "source": "#fgbio_separate_bams.cwl/fgbio_filter_consensus_reads_1_2_0_duplex/fgbio_filter_consensus_reads_bam"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/output_file_name",
                            "source": "#fgbio_separate_bams.cwl/output_file_name_duplex_aln_metrics"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/reference",
                            "source": "#fgbio_separate_bams.cwl/reference_fasta"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/validation_stringency",
                            "source": "#fgbio_separate_bams.cwl/validation_stringency"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/assume_sorted",
                            "source": "#fgbio_separate_bams.cwl/assume_sorted"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/create_index",
                            "source": "#fgbio_separate_bams.cwl/create_index"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/use_jdk_deflater",
                            "source": "#fgbio_separate_bams.cwl/use_jdk_deflater"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/use_jdk_inflater",
                            "source": "#fgbio_separate_bams.cwl/use_jdk_inflater"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/temporary_directory",
                            "source": "#fgbio_separate_bams.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_duplex/gatk_collect_alignment_summary_metrics_txt"
                        }
                    ],
                    "run": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2",
                    "label": "GATK-CollectAlignmentSummaryMetrics",
                    "https://www.sevenbridges.com/x": 1072.9705810546875,
                    "https://www.sevenbridges.com/y": 1424.1484375
                },
                {
                    "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex",
                    "in": [
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex/input",
                            "source": "#fgbio_separate_bams.cwl/fgbio_postprocessing_simplex_filter_0_1_8/fgbio_postprocessing_simplex_bam"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex/output_file_name",
                            "source": "#fgbio_separate_bams.cwl/output_file_name_simplex_aln_metrics"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex/reference",
                            "source": "#fgbio_separate_bams.cwl/reference_fasta"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex/validation_stringency",
                            "source": "#fgbio_separate_bams.cwl/validation_stringency"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex/assume_sorted",
                            "source": "#fgbio_separate_bams.cwl/assume_sorted"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex/create_index",
                            "source": "#fgbio_separate_bams.cwl/create_index"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex/use_jdk_deflater",
                            "source": "#fgbio_separate_bams.cwl/use_jdk_deflater"
                        },
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex/use_jdk_inflater",
                            "source": "#fgbio_separate_bams.cwl/use_jdk_inflater"
                        }
                    ],
                    "out": [
                        {
                            "id": "#fgbio_separate_bams.cwl/gatk_collect_alignment_summary_metrics_4.1.8.0_simplex/gatk_collect_alignment_summary_metrics_txt"
                        }
                    ],
                    "run": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_2",
                    "label": "GATK-CollectAlignmentSummaryMetrics",
                    "https://www.sevenbridges.com/x": 1616.9268798828125,
                    "https://www.sevenbridges.com/y": 1654.25
                }
            ],
            "requirements": [],
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9042-6213",
                    "https://schema.org/name": "Ronak Shah"
                }
            ],
            "https://schema.org/citation": "",
            "https://schema.org/codeRepository": "https://github.com/msk-access/cwl_subworkflows/fgbio_separate_bams",
            "https://schema.org/contributor": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9042-6213",
                    "https://schema.org/name": "Ronak Shah"
                }
            ],
            "https://schema.org/dateCreated": "2020-06-09",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0"
        },
        {
            "class": "Workflow",
            "id": "#indel_realignment.cwl",
            "label": "indel_realignment",
            "inputs": [
                {
                    "id": "#indel_realignment.cwl/window_size",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 0
                },
                {
                    "id": "#indel_realignment.cwl/soft_clip_contig",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 426.796875
                },
                {
                    "id": "#indel_realignment.cwl/scoring_gap_alignments",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 533.53125
                },
                {
                    "id": "#indel_realignment.cwl/reference_fasta",
                    "type": "File",
                    "secondaryFiles": [
                        ".fai"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 640.21875
                },
                {
                    "id": "#indel_realignment.cwl/no_sort",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1066.875
                },
                {
                    "id": "#indel_realignment.cwl/maximum_mixmatch_rate",
                    "type": [
                        "null",
                        "float"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1280.25
                },
                {
                    "id": "#indel_realignment.cwl/maximum_average_depth",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1386.9375
                },
                {
                    "id": "#indel_realignment.cwl/input_bam",
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1493.625
                },
                {
                    "id": "#indel_realignment.cwl/ignore_bad_assembly",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1600.3125
                },
                {
                    "id": "#indel_realignment.cwl/contig_anchor",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1813.6875
                },
                {
                    "id": "#indel_realignment.cwl/consensus_sequence",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1920.375
                },
                {
                    "id": "#indel_realignment.cwl/bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2027.015625
                },
                {
                    "id": "#indel_realignment.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 960.234375
                },
                {
                    "id": "#indel_realignment.cwl/option_bedgraph",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 853.546875
                },
                {
                    "id": "#indel_realignment.cwl/no_edge_complex_indel",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1173.5625
                },
                {
                    "id": "#indel_realignment.cwl/distance_between_features",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1707
                },
                {
                    "id": "#indel_realignment.cwl/output_bams",
                    "type": [
                        "string",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 746.859375
                },
                {
                    "id": "#indel_realignment.cwl/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 106.6875
                },
                {
                    "id": "#indel_realignment.cwl/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 320.109375
                },
                {
                    "id": "#indel_realignment.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 992.927978515625,
                    "https://www.sevenbridges.com/y": 794.8671875
                },
                {
                    "id": "#indel_realignment.cwl/create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 992.927978515625,
                    "https://www.sevenbridges.com/y": 901.5078125
                },
                {
                    "id": "#indel_realignment.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 213.421875
                }
            ],
            "outputs": [
                {
                    "id": "#indel_realignment.cwl/indel_realignment_bam",
                    "outputSource": [
                        "#indel_realignment.cwl/picard_fix_mate_information_4_1_8_1/picard_fix_mate_information_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 1981.323974609375,
                    "https://www.sevenbridges.com/y": 1013.4609375
                }
            ],
            "steps": [
                {
                    "id": "#indel_realignment.cwl/abra2_2_22",
                    "in": [
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/number_of_threads",
                            "source": "#indel_realignment.cwl/number_of_threads"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/input_bam",
                            "source": [
                                "#indel_realignment.cwl/input_bam"
                            ]
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/working_directory",
                            "source": "#indel_realignment.cwl/temporary_directory"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/reference_fasta",
                            "source": "#indel_realignment.cwl/reference_fasta"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/targets",
                            "source": "#indel_realignment.cwl/bedtools_merge/bedtools_merge_bed"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/maximum_average_depth",
                            "source": "#indel_realignment.cwl/maximum_average_depth"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/soft_clip_contig",
                            "source": "#indel_realignment.cwl/soft_clip_contig"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/maximum_mixmatch_rate",
                            "source": "#indel_realignment.cwl/maximum_mixmatch_rate"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/scoring_gap_alignments",
                            "source": "#indel_realignment.cwl/scoring_gap_alignments"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/contig_anchor",
                            "source": "#indel_realignment.cwl/contig_anchor"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/window_size",
                            "source": "#indel_realignment.cwl/window_size"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/consensus_sequence",
                            "source": "#indel_realignment.cwl/consensus_sequence"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/output_bams",
                            "source": [
                                "#indel_realignment.cwl/output_bams"
                            ]
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/ignore_bad_assembly",
                            "source": "#indel_realignment.cwl/ignore_bad_assembly"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/bam_index",
                            "source": "#indel_realignment.cwl/bam_index"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/no_edge_complex_indel",
                            "source": "#indel_realignment.cwl/no_edge_complex_indel"
                        },
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/no_sort",
                            "source": "#indel_realignment.cwl/no_sort"
                        }
                    ],
                    "out": [
                        {
                            "id": "#indel_realignment.cwl/abra2_2_22/abra_realigned_bam"
                        }
                    ],
                    "run": "#abra2_2.22.cwl",
                    "label": "abra2_2.22",
                    "https://www.sevenbridges.com/x": 992.927978515625,
                    "https://www.sevenbridges.com/y": 1120.1484375
                },
                {
                    "id": "#indel_realignment.cwl/bedtools_genomecov",
                    "in": [
                        {
                            "id": "#indel_realignment.cwl/bedtools_genomecov/input",
                            "source": "#indel_realignment.cwl/input_bam"
                        },
                        {
                            "id": "#indel_realignment.cwl/bedtools_genomecov/option_bedgraph",
                            "source": "#indel_realignment.cwl/option_bedgraph"
                        }
                    ],
                    "out": [
                        {
                            "id": "#indel_realignment.cwl/bedtools_genomecov/bedtools_genomecove_bedgraph"
                        }
                    ],
                    "run": "#bedtools_genomecov_v2.28.0_cv2.cwl",
                    "label": "bedtools_genomecov",
                    "https://www.sevenbridges.com/x": 269.59375,
                    "https://www.sevenbridges.com/y": 1006.4609375
                },
                {
                    "id": "#indel_realignment.cwl/bedtools_merge",
                    "in": [
                        {
                            "id": "#indel_realignment.cwl/bedtools_merge/input",
                            "source": "#indel_realignment.cwl/bedtools_genomecov/bedtools_genomecove_bedgraph"
                        },
                        {
                            "id": "#indel_realignment.cwl/bedtools_merge/distance_between_features",
                            "source": "#indel_realignment.cwl/distance_between_features"
                        }
                    ],
                    "out": [
                        {
                            "id": "#indel_realignment.cwl/bedtools_merge/bedtools_merge_bed"
                        }
                    ],
                    "run": "#bedtools_merge_v2.28.0_cv2.cwl",
                    "label": "bedtools_merge",
                    "https://www.sevenbridges.com/x": 635.5108642578125,
                    "https://www.sevenbridges.com/y": 1006.4609375
                },
                {
                    "id": "#indel_realignment.cwl/picard_fix_mate_information_4_1_8_1",
                    "in": [
                        {
                            "id": "#indel_realignment.cwl/picard_fix_mate_information_4_1_8_1/input",
                            "source": "#indel_realignment.cwl/abra2_2_22/abra_realigned_bam"
                        },
                        {
                            "id": "#indel_realignment.cwl/picard_fix_mate_information_4_1_8_1/output_file_name",
                            "source": "#indel_realignment.cwl/output_file_name"
                        },
                        {
                            "id": "#indel_realignment.cwl/picard_fix_mate_information_4_1_8_1/sort_order",
                            "source": "#indel_realignment.cwl/sort_order"
                        },
                        {
                            "id": "#indel_realignment.cwl/picard_fix_mate_information_4_1_8_1/validation_stringency",
                            "source": "#indel_realignment.cwl/validation_stringency"
                        },
                        {
                            "id": "#indel_realignment.cwl/picard_fix_mate_information_4_1_8_1/create_bam_index",
                            "source": "#indel_realignment.cwl/create_bam_index"
                        },
                        {
                            "id": "#indel_realignment.cwl/picard_fix_mate_information_4_1_8_1/temporary_directory",
                            "source": "#indel_realignment.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#indel_realignment.cwl/picard_fix_mate_information_4_1_8_1/picard_fix_mate_information_bam"
                        }
                    ],
                    "run": "#picard_fix_mate_information_4.1.8.1.cwl",
                    "label": "picard_fix_mate_information_4.1.8.1",
                    "https://www.sevenbridges.com/x": 1546.70458984375,
                    "https://www.sevenbridges.com/y": 978.328125
                }
            ],
            "requirements": [],
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:murphyc4@mskcc.org",
                    "https://schema.org/name": "Charlie Murphy"
                }
            ],
            "https://schema.org/citation": "",
            "https://schema.org/codeRepository": "https://github.com/msk-access/cwl_subworkflows/indel_realignment",
            "https://schema.org/contributor": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:murphyc4@mskcc.org",
                    "https://schema.org/name": "Charlie Murphy"
                }
            ],
            "https://schema.org/dateCreated": "2020-09-14",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0"
        },
        {
            "class": "CommandLineTool",
            "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3",
            "baseCommand": [
                "gatk",
                "CollectAlignmentSummaryMetrics"
            ],
            "inputs": [
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-I"
                    },
                    "doc": "Input file (bam or sam).  Required."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "File to write the output to.  Required."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/reference",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-R"
                    },
                    "doc": "Reference sequence file. Note that while this argument is not required, without it only a small subset of the metrics will be calculated. Note also that if a reference sequence is provided, it must be accompanied by a sequence dictionary.  Default value: null.",
                    "secondaryFiles": [
                        "^.fasta.fai",
                        "^.dict"
                    ]
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/adaptor_sequence",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ADAPTER_SEQUENCE"
                    },
                    "doc": "List of adapter sequences to use when processing the alignment metrics.  This argument may be specified 0 or more times. Default value: [AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG]."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/metrics_acciumulation_level",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--METRIC_ACCUMULATION_LEVEL"
                    },
                    "doc": "The level(s) at which to accumulate metrics. Default value: [ALL_READS]. This option can be set to 'null' to clear the default value. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} This option may be specified 0 or more times. This option can be set to 'null' to clear the default list."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/expected_pair_orientations",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--EXPECTED_PAIR_ORIENTATIONS"
                    },
                    "doc": "Paired-end reads that do not have this expected orientation will be considered chimeric. This argument may be specified 0 or more times. Default value: [FR]. Possible values: {FR, RF, TANDEM}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/is_bisulfite_sequenced",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--IS_BISULFITE_SEQUENCED"
                    },
                    "doc": "Whether the SAM or BAM file consists of bisulfite sequenced reads.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/max_insert_size",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MAX_INSERT_SIZE"
                    },
                    "doc": "Paired-end reads above this insert size will be considered chimeric along with inter-chromosomal pairs.  Default value: 100000."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT,LENIENT, SILENT}"
                },
                {
                    "default": true,
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/assume_sorted",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ASSUME_SORTED"
                    },
                    "doc": "If true (default), then the sort order in the header file will be ignored.  Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/stop_after",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--STOP_AFTER"
                    },
                    "doc": "Stop after processing N reads, mainly for debugging. Default value: 0. This option can be set to 'null' to clear the default value."
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/create_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/create_md5_file",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_MD5_FILE"
                    },
                    "doc": "Whether to create an MD5 digest for any BAM or FASTQ files created.    Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3/gatk_collect_alignment_summary_metrics_txt",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if (inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/, '_alignment_summary_metrics.txt')\n    }\n}"
                    }
                }
            ],
            "label": "GATK-CollectAlignmentSummaryMetrics",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--java-options",
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx15G\"\n  }\n  else {\n      return \"-Xmx15G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/, '_alignment_summary_metrics.txt')\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 32000,
                    "coresMin": 1
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "gatk4",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.0"
                }
            ]
        },
        {
            "class": "Workflow",
            "id": "#main",
            "doc": "This workflow takes a READ1 and READ2 fastq.gz file generated for MSK-ACCESS assay and generated four different Binary Alignment Map file along with alignment metrics for each.",
            "label": "nucleo",
            "inputs": [
                {
                    "id": "#reference_sequence",
                    "type": "File",
                    "doc": "Reference sequence file.  Please include \".fai\", \"^.dict\", \".amb\" , \".sa\", \".bwt\", \".pac\", \".ann\" as  secondary files if they are not present in the same location as the  \".fasta\" file",
                    "secondaryFiles": [
                        ".amb",
                        ".fai",
                        ".sa",
                        "^.dict",
                        ".ann",
                        ".bwt",
                        ".pac"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1388.765625
                },
                {
                    "id": "#gatk_base_recalibrator_known_sites",
                    "type": {
                        "type": "array",
                        "items": "File",
                        "inputBinding": {
                            "prefix": "--known-sites"
                        }
                    },
                    "secondaryFiles": [
                        ".idx"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2991.1875
                },
                {
                    "id": "#sequencing-center",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1068.28125
                },
                {
                    "id": "#run-date",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1281.9375
                },
                {
                    "id": "#sample",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1175.109375
                },
                {
                    "id": "#read-structures",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1495.59375
                },
                {
                    "id": "#read-group-id",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1602.421875
                },
                {
                    "id": "#platform-unit",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1709.25
                },
                {
                    "id": "#platform-model",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1816.078125
                },
                {
                    "id": "#platform",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1922.90625
                },
                {
                    "id": "#library",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2350.21875
                },
                {
                    "id": "#validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 0
                },
                {
                    "id": "#UBG_picard_SamToFastq_R1_output_fastq",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 213.65625
                },
                {
                    "id": "#UBG_picard_SamToFastq_R2_output_fastq",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 106.828125
                },
                {
                    "id": "#fastp_read2_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5341.40625
                },
                {
                    "id": "#fastp_read2_adapter_sequence",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5448.234375
                },
                {
                    "id": "#fastp_read1_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5555.0625
                },
                {
                    "id": "#fastp_read1_adapter_sequence",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5661.890625
                },
                {
                    "id": "#fastp_minimum_read_length",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5768.71875
                },
                {
                    "id": "#fastp_html_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5982.375
                },
                {
                    "id": "#fastp_json_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5875.546875
                },
                {
                    "id": "#sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 961.453125
                },
                {
                    "id": "#bwa_mem_T",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 6516.515625
                },
                {
                    "id": "#bwa_mem_Y",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 6409.6875
                },
                {
                    "id": "#UBG_picard_addRG_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 427.3125
                },
                {
                    "id": "#create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 6302.859375
                },
                {
                    "id": "#UBG_bwa_mem_output",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 640.96875
                },
                {
                    "id": "#bwa_mem_K",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 6623.34375
                },
                {
                    "id": "#UBG_gatk_merge_bam_alignment_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 534.140625
                },
                {
                    "id": "#optical_duplicate_pixel_distance",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2136.5625
                },
                {
                    "id": "#gatk_mark_duplicates_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2563.875
                },
                {
                    "id": "#gatk_mark_duplicates_duplication_metrics_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2670.703125
                },
                {
                    "id": "#bedtools_merge_distance_between_features",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 6730.171875
                },
                {
                    "id": "#bedtools_genomecov_option_bedgraph",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 6837
                },
                {
                    "id": "#apply_bqsr_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 611.546875,
                    "https://www.sevenbridges.com/y": 4433.3671875
                },
                {
                    "id": "#abra2_window_size",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 7798.453125
                },
                {
                    "id": "#abra2_soft_clip_contig",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 7905.28125
                },
                {
                    "id": "#abra2_scoring_gap_alignments",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 8012.109375
                },
                {
                    "id": "#UBG_abra2_output_bams",
                    "type": [
                        "string",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 747.796875
                },
                {
                    "id": "#abra2_no_edge_complex_indel",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 8225.765625
                },
                {
                    "id": "#abra2_maximum_mixmatch_rate",
                    "type": [
                        "null",
                        "float"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 8332.59375
                },
                {
                    "id": "#abra2_maximum_average_depth",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 8439.421875
                },
                {
                    "id": "#abra2_consensus_sequence",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 8653.078125
                },
                {
                    "id": "#abra2_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 8759.90625
                },
                {
                    "id": "#fgbio_fastq_to_bam_input",
                    "type": {
                        "type": "array",
                        "items": {
                            "items": "File",
                            "type": "array"
                        }
                    },
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4593.609375
                },
                {
                    "id": "#UBG_picard_fixmateinformation_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 320.484375
                },
                {
                    "id": "#merge_sam_files_sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2243.390625
                },
                {
                    "id": "#gatk_merge_sam_files_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2457.046875
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_intervals",
                    "type": [
                        "null",
                        "File"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4807.265625
                },
                {
                    "id": "#fgbio_group_reads_by_umi_strategy",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3311.671875
                },
                {
                    "id": "#fgbio_group_reads_by_umi_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3418.5
                },
                {
                    "id": "#fgbio_group_reads_by_umi_family_size_histogram",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3525.328125
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_output_prefix",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4700.4375
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_duplex_umi_counts",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4914.09375
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5020.921875
                },
                {
                    "id": "#fgbio_call_duplex_consensus_reads_min_reads",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5127.75
                },
                {
                    "id": "#BC_gatk_sam_to_fastq_output_name_R2",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 7157.484375
                },
                {
                    "id": "#BC_gatk_sam_to_fastq_output_name_R1",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 7264.3125
                },
                {
                    "id": "#BC_picard_fixmate_information_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 6943.828125
                },
                {
                    "id": "#BC_picard_addRG_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 7050.65625
                },
                {
                    "id": "#BC_gatk_merge_bam_alignment_output_file_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 7371.140625
                },
                {
                    "id": "#fgbio_postprocessing_output_file_name_simplex",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3204.84375
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2884.359375
                },
                {
                    "id": "#abra2_contig_anchor",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 8546.25
                },
                {
                    "id": "#BC_abra2_output_bams",
                    "type": [
                        "string",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 7584.796875
                },
                {
                    "id": "#fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3632.15625
                },
                {
                    "id": "#fgbio_filter_consensus_read_output_file_name_simplex_duplex",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3738.984375
                },
                {
                    "id": "#fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3845.8125
                },
                {
                    "id": "#fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3952.640625
                },
                {
                    "id": "#fgbio_filter_consensus_read_output_file_name_duplex",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4059.46875
                },
                {
                    "id": "#fgbio_filter_consensus_read_min_reads_duplex",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4273.125
                },
                {
                    "id": "#fgbio_filter_consensus_read_min_base_quality_simplex_duplex",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4379.953125
                },
                {
                    "id": "#fgbio_filter_consensus_read_min_base_quality_duplex",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4486.78125
                },
                {
                    "id": "#BC_bwa_mem_output",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 7477.96875
                },
                {
                    "id": "#fgbio_filter_consensus_read_min_reads_simplex_duplex",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4166.296875
                },
                {
                    "id": "#picard_addRG_sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2029.734375
                },
                {
                    "id": "#disable_trim_poly_g",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 6089.203125
                },
                {
                    "id": "#disable_quality_filtering",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 6196.03125
                },
                {
                    "id": "#gatk_base_recalibrator_add_output_sam_program_record",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3098.015625
                },
                {
                    "id": "#gatk_collect_aln_summary_metrics_bqsr_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2777.53125
                },
                {
                    "id": "#abra2_no_sort",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 8118.9375
                },
                {
                    "id": "#base_recalibrator_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 7691.625
                },
                {
                    "id": "#temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 854.625
                },
                {
                    "id": "#fgbio_async_io",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5234.578125
                }
            ],
            "outputs": [
                {
                    "id": "#fastp_html_output",
                    "outputSource": [
                        "#uncollapsed_bam_generation/fastp_html_output"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1871.76708984375,
                    "https://www.sevenbridges.com/y": 4006.3671875
                },
                {
                    "id": "#fastp_json_output",
                    "outputSource": [
                        "#uncollapsed_bam_generation/fastp_json_output"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1871.76708984375,
                    "https://www.sevenbridges.com/y": 3899.5390625
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_txt_uncollapsed",
                    "outputSource": [
                        "#gatk_collect_alignment_summary_metrics_4_1_8_0/gatk_collect_alignment_summary_metrics_txt"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3862.5732421875,
                    "https://www.sevenbridges.com/y": 4379.953125
                },
                {
                    "id": "#indel_realignment_bam",
                    "outputSource": [
                        "#uncollapsed_bam_generation/indel_realignment_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 1871.76708984375,
                    "https://www.sevenbridges.com/y": 3792.7109375
                },
                {
                    "id": "#picard_mark_duplicates_metrics",
                    "outputSource": [
                        "#uncollapsed_bam_generation/picard_mark_duplicates_metrics"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 1871.76708984375,
                    "https://www.sevenbridges.com/y": 3685.8828125
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_txt_simplex",
                    "outputSource": [
                        "#bam_collapsing/gatk_collect_alignment_summary_metrics_txt_simplex"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 3664.5703125
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_txt_duplex",
                    "outputSource": [
                        "#bam_collapsing/gatk_collect_alignment_summary_metrics_txt_duplex"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 3771.3984375
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_txt_collapsed",
                    "outputSource": [
                        "#bam_collapsing/gatk_collect_alignment_summary_metrics_txt_collapsed"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 3878.2265625
                },
                {
                    "id": "#fgbio_postprocessing_simplex_bam",
                    "outputSource": [
                        "#bam_collapsing/fgbio_postprocessing_simplex_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 4133.8828125
                },
                {
                    "id": "#fgbio_group_reads_by_umi_histogram",
                    "outputSource": [
                        "#bam_collapsing/fgbio_group_reads_by_umi_histogram"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 4240.7109375
                },
                {
                    "id": "#fgbio_group_reads_by_umi_bam",
                    "outputSource": [
                        "#bam_collapsing/fgbio_group_reads_by_umi_bam"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 4347.5390625
                },
                {
                    "id": "#fgbio_filter_consensus_reads_duplex_bam",
                    "outputSource": [
                        "#bam_collapsing/fgbio_filter_consensus_reads_duplex_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 4454.3671875
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_umi_counts",
                    "outputSource": [
                        "#bam_collapsing/fgbio_collect_duplex_seq_metrics_umi_counts"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 4561.1953125
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_family_size",
                    "outputSource": [
                        "#bam_collapsing/fgbio_collect_duplex_seq_metrics_family_size"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 4668.0234375
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_duplex_yield_metrics",
                    "outputSource": [
                        "#bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_yield_metrics"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 4774.8515625
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt",
                    "outputSource": [
                        "#bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 4881.6796875
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_duplex_qc",
                    "outputSource": [
                        "#bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_qc"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 4988.5078125
                },
                {
                    "id": "#fgbio_collect_duplex_seq_metrics_duplex_family_size",
                    "outputSource": [
                        "#bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_family_size"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 5095.3359375
                },
                {
                    "id": "#fgbio_collapsed_bam",
                    "outputSource": [
                        "#bam_collapsing/fgbio_collapsed_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 5202.1640625
                },
                {
                    "id": "#uncollapsed_bam",
                    "outputSource": [
                        "#base_quality_recalibration/gatk_apply_bqsr_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 3557.7421875
                }
            ],
            "steps": [
                {
                    "id": "#uncollapsed_bam_generation",
                    "in": [
                        {
                            "id": "#uncollapsed_bam_generation/sequencing-center",
                            "default": "MSKCC",
                            "source": "#sequencing-center"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/sample",
                            "source": "#sample"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/run-date",
                            "source": "#run-date"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/read-structures",
                            "default": [
                                "3M2S+T",
                                "3M2S+T"
                            ],
                            "source": [
                                "#read-structures"
                            ]
                        },
                        {
                            "id": "#uncollapsed_bam_generation/read-group-id",
                            "source": "#read-group-id"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/platform-unit",
                            "source": "#platform-unit"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/platform-model",
                            "default": "novaseq",
                            "source": "#platform-model"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/platform",
                            "default": "ILLUMINA",
                            "source": "#platform"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/library",
                            "source": "#library"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/validation_stringency",
                            "default": "LENIENT",
                            "source": "#validation_stringency"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/R1_output_fastq",
                            "source": "#UBG_picard_SamToFastq_R1_output_fastq"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/R2_output_fastq",
                            "source": "#UBG_picard_SamToFastq_R2_output_fastq"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/reference_sequence",
                            "source": "#reference_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_read2_output_file_name",
                            "source": "#fastp_read2_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_read2_adapter_sequence",
                            "default": "AGATCGGAAGAGC",
                            "source": "#fastp_read2_adapter_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_read1_output_file_name",
                            "source": "#fastp_read1_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_read1_adapter_sequence",
                            "default": "GATCGGAAGAGC",
                            "source": "#fastp_read1_adapter_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_minimum_read_length",
                            "default": 25,
                            "source": "#fastp_minimum_read_length"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_json_output_file_name",
                            "source": "#fastp_json_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_html_output_file_name",
                            "source": "#fastp_html_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/bwa_mem_Y",
                            "default": true,
                            "source": "#bwa_mem_Y"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/bwa_mem_T",
                            "source": "#bwa_mem_T"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/sort_order",
                            "default": "coordinate",
                            "source": "#sort_order"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/picard_addRG_output_file_name",
                            "source": "#UBG_picard_addRG_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/bwa_mem_output",
                            "source": "#UBG_bwa_mem_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/bwa_mem_K",
                            "source": "#bwa_mem_K"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/create_bam_index",
                            "default": true,
                            "source": "#create_bam_index"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/gatk_merge_bam_alignment_output_file_name",
                            "source": "#UBG_gatk_merge_bam_alignment_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/optical_duplicate_pixel_distance",
                            "default": 2500,
                            "source": "#optical_duplicate_pixel_distance"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/gatk_mark_duplicates_output_file_name",
                            "source": "#gatk_mark_duplicates_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/gatk_mark_duplicates_duplication_metrics_file_name",
                            "source": "#gatk_mark_duplicates_duplication_metrics_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_window_size",
                            "default": "800,700",
                            "source": "#abra2_window_size"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_soft_clip_contig",
                            "default": "100,30,80,15",
                            "source": "#abra2_soft_clip_contig"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_scoring_gap_alignments",
                            "default": "8,32,48,1",
                            "source": "#abra2_scoring_gap_alignments"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_output_bams",
                            "source": [
                                "#UBG_abra2_output_bams"
                            ]
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_maximum_average_depth",
                            "default": 1000,
                            "source": "#abra2_maximum_average_depth"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_bam_index",
                            "source": "#abra2_bam_index"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_contig_anchor",
                            "default": "10,1",
                            "source": "#abra2_contig_anchor"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_consensus_sequence",
                            "default": true,
                            "source": "#abra2_consensus_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/bedtools_merge_distance_between_features",
                            "default": 10,
                            "source": "#bedtools_merge_distance_between_features"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_maximum_mixmatch_rate",
                            "default": 0.1,
                            "source": "#abra2_maximum_mixmatch_rate"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/bedtools_genomecov_option_bedgraph",
                            "default": true,
                            "source": "#bedtools_genomecov_option_bedgraph"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/picard_fixmateinformation_output_file_name",
                            "source": "#UBG_picard_fixmateinformation_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_no_sort",
                            "default": true,
                            "source": "#abra2_no_sort"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/abra2_no_edge_complex_indel",
                            "default": true,
                            "source": "#abra2_no_edge_complex_indel"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/merge_sam_files_sort_order",
                            "default": "queryname",
                            "source": "#merge_sam_files_sort_order"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/gatk_merge_sam_files_output_file_name",
                            "source": "#gatk_merge_sam_files_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fgbio_fastq_to_bam_input",
                            "source": [
                                "#fgbio_fastq_to_bam_input"
                            ]
                        },
                        {
                            "id": "#uncollapsed_bam_generation/picard_addRG_sort_order",
                            "default": "queryname",
                            "source": "#picard_addRG_sort_order"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/disable_trim_poly_g",
                            "default": true,
                            "source": "#disable_trim_poly_g"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/disable_quality_filtering",
                            "default": true,
                            "source": "#disable_quality_filtering"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/temporary_directory",
                            "source": "#temporary_directory"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fgbio_async_io",
                            "default": "true",
                            "source": "#fgbio_async_io"
                        }
                    ],
                    "out": [
                        {
                            "id": "#uncollapsed_bam_generation/gatk_sam_to_fastq_unpaired_fastq"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_unpaired2_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_unpaired1_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_json_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/fastp_html_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/picard_mark_duplicates_metrics"
                        },
                        {
                            "id": "#uncollapsed_bam_generation/indel_realignment_bam"
                        }
                    ],
                    "run": "#uncollapsed_bam_generation.cwl",
                    "label": "Uncollapsed BAM Generation",
                    "https://www.sevenbridges.com/x": 611.546875,
                    "https://www.sevenbridges.com/y": 3962.5390625
                },
                {
                    "id": "#bam_collapsing",
                    "in": [
                        {
                            "id": "#bam_collapsing/fgbio_group_reads_by_umi_input",
                            "source": "#uncollapsed_bam_generation/indel_realignment_bam"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_group_reads_by_umi_strategy",
                            "default": "paired",
                            "source": "#fgbio_group_reads_by_umi_strategy"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_group_reads_by_umi_output_file_name",
                            "source": "#fgbio_group_reads_by_umi_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_group_reads_by_umi_family_size_histogram",
                            "source": "#fgbio_group_reads_by_umi_family_size_histogram"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_collect_duplex_seq_metrics_intervals",
                            "source": "#fgbio_collect_duplex_seq_metrics_intervals"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_collect_duplex_seq_metrics_output_prefix",
                            "source": "#fgbio_collect_duplex_seq_metrics_output_prefix"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_umi_counts",
                            "default": true,
                            "source": "#fgbio_collect_duplex_seq_metrics_duplex_umi_counts"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_call_duplex_consensus_reads_read_group_id",
                            "source": "#read-group-id"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_call_duplex_consensus_reads_output_file_name",
                            "source": "#fgbio_call_duplex_consensus_reads_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_call_duplex_consensus_reads_min_reads",
                            "default": [
                                1,
                                1,
                                0
                            ],
                            "source": [
                                "#fgbio_call_duplex_consensus_reads_min_reads"
                            ]
                        },
                        {
                            "id": "#bam_collapsing/reference_sequence",
                            "source": "#reference_sequence"
                        },
                        {
                            "id": "#bam_collapsing/validation_stringency",
                            "default": "LENIENT",
                            "source": "#validation_stringency"
                        },
                        {
                            "id": "#bam_collapsing/gatk_sam_to_fastq_output_name_R2",
                            "source": "#BC_gatk_sam_to_fastq_output_name_R2"
                        },
                        {
                            "id": "#bam_collapsing/gatk_sam_to_fastq_output_name_R1",
                            "source": "#BC_gatk_sam_to_fastq_output_name_R1"
                        },
                        {
                            "id": "#bam_collapsing/bwa_mem_Y",
                            "default": true,
                            "source": "#bwa_mem_Y"
                        },
                        {
                            "id": "#bam_collapsing/bwa_mem_T",
                            "source": "#bwa_mem_T"
                        },
                        {
                            "id": "#bam_collapsing/sort_order",
                            "default": "coordinate",
                            "source": "#sort_order"
                        },
                        {
                            "id": "#bam_collapsing/picard_addRG_read_group_sequencing_platform",
                            "default": "ILLUMINA",
                            "source": "#platform"
                        },
                        {
                            "id": "#bam_collapsing/picard_addRG_read_group_sequencing_center",
                            "default": "MSKCC",
                            "source": "#sequencing-center"
                        },
                        {
                            "id": "#bam_collapsing/picard_addRG_read_group_run_date",
                            "source": "#run-date"
                        },
                        {
                            "id": "#bam_collapsing/picard_addRG_read_group_platform_unit",
                            "source": "#platform-unit"
                        },
                        {
                            "id": "#bam_collapsing/picard_addRG_read_group_library",
                            "source": "#library"
                        },
                        {
                            "id": "#bam_collapsing/picard_addRG_read_group_identifier",
                            "source": "#read-group-id"
                        },
                        {
                            "id": "#bam_collapsing/picard_addRG_output_file_name",
                            "source": "#BC_picard_addRG_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing/bwa_mem_output",
                            "source": "#BC_bwa_mem_output"
                        },
                        {
                            "id": "#bam_collapsing/create_bam_index",
                            "default": true,
                            "source": "#create_bam_index"
                        },
                        {
                            "id": "#bam_collapsing/bwa_mem_K",
                            "source": "#bwa_mem_K"
                        },
                        {
                            "id": "#bam_collapsing/abra2_window_size",
                            "default": "800,700",
                            "source": "#abra2_window_size"
                        },
                        {
                            "id": "#bam_collapsing/abra2_soft_clip_contig",
                            "default": "100,30,80,15",
                            "source": "#abra2_soft_clip_contig"
                        },
                        {
                            "id": "#bam_collapsing/abra2_scoring_gap_alignments",
                            "default": "8,32,48,1",
                            "source": "#abra2_scoring_gap_alignments"
                        },
                        {
                            "id": "#bam_collapsing/picard_fixmate_information_output_file_name",
                            "source": "#BC_picard_fixmate_information_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing/abra2_output_bams",
                            "source": [
                                "#BC_abra2_output_bams"
                            ]
                        },
                        {
                            "id": "#bam_collapsing/bedtools_genomecov_option_bedgraph",
                            "default": true,
                            "source": "#bedtools_genomecov_option_bedgraph"
                        },
                        {
                            "id": "#bam_collapsing/abra2_no_sort",
                            "default": true,
                            "source": "#abra2_no_sort"
                        },
                        {
                            "id": "#bam_collapsing/abra2_no_edge_complex_indel",
                            "default": true,
                            "source": "#abra2_no_edge_complex_indel"
                        },
                        {
                            "id": "#bam_collapsing/abra2_maximum_mixmatch_rate",
                            "default": 0.1,
                            "source": "#abra2_maximum_mixmatch_rate"
                        },
                        {
                            "id": "#bam_collapsing/abra2_maximum_average_depth",
                            "default": 1000,
                            "source": "#abra2_maximum_average_depth"
                        },
                        {
                            "id": "#bam_collapsing/bedtools_merge_distance_between_features",
                            "default": 10,
                            "source": "#bedtools_merge_distance_between_features"
                        },
                        {
                            "id": "#bam_collapsing/abra2_contig_anchor",
                            "default": "10,1",
                            "source": "#abra2_contig_anchor"
                        },
                        {
                            "id": "#bam_collapsing/abra2_consensus_sequence",
                            "default": true,
                            "source": "#abra2_consensus_sequence"
                        },
                        {
                            "id": "#bam_collapsing/gatk_merge_bam_alignment_output_file_name",
                            "source": "#BC_gatk_merge_bam_alignment_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex",
                            "default": true,
                            "source": "#fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_read_reverse_per_base_tags_duplex",
                            "default": true
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_read_min_base_quality_duplex",
                            "default": 30,
                            "source": "#fgbio_filter_consensus_read_min_base_quality_duplex"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_read_min_base_quality_simplex_duplex",
                            "default": 30,
                            "source": "#fgbio_filter_consensus_read_min_base_quality_simplex_duplex"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_read_min_reads_duplex",
                            "default": [
                                2,
                                1,
                                1
                            ],
                            "source": [
                                "#fgbio_filter_consensus_read_min_reads_duplex"
                            ]
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_read_min_reads_simplex_duplex",
                            "default": [
                                3,
                                3,
                                0
                            ],
                            "source": [
                                "#fgbio_filter_consensus_read_min_reads_simplex_duplex"
                            ]
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_read_output_file_name_simplex_duplex",
                            "source": "#fgbio_filter_consensus_read_output_file_name_simplex_duplex"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics",
                            "source": "#fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_postprocessing_output_file_name_simplex",
                            "source": "#fgbio_postprocessing_output_file_name_simplex"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics",
                            "source": "#fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_read_output_file_name_duplex",
                            "source": "#fgbio_filter_consensus_read_output_file_name_duplex"
                        },
                        {
                            "id": "#bam_collapsing/picard_addRG_read_group_sample_name",
                            "source": "#sample"
                        },
                        {
                            "id": "#bam_collapsing/gatk_collect_alignment_summary_metrics_output_file_name",
                            "source": "#gatk_collect_alignment_summary_metrics_output_file_name"
                        },
                        {
                            "id": "#bam_collapsing/picard_addRG_sort_order",
                            "default": "queryname",
                            "source": "#picard_addRG_sort_order"
                        },
                        {
                            "id": "#bam_collapsing/temporary_directory",
                            "source": "#temporary_directory"
                        },
                        {
                            "id": "#bam_collapsing/async_io",
                            "default": "true",
                            "source": "#fgbio_async_io"
                        }
                    ],
                    "out": [
                        {
                            "id": "#bam_collapsing/fgbio_group_reads_by_umi_histogram"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_collect_duplex_seq_metrics_umi_counts"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_collect_duplex_seq_metrics_family_size"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_yield_metrics"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_umi_counts_txt"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_qc"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_collect_duplex_seq_metrics_duplex_family_size"
                        },
                        {
                            "id": "#bam_collapsing/gatk_sam_to_fastq_unpaired_fastq"
                        },
                        {
                            "id": "#bam_collapsing/gatk_sam_to_fastq_second_end_fastq"
                        },
                        {
                            "id": "#bam_collapsing/gatk_sam_to_fastq_fastq"
                        },
                        {
                            "id": "#bam_collapsing/gatk_collect_alignment_summary_metrics_txt_simplex"
                        },
                        {
                            "id": "#bam_collapsing/gatk_collect_alignment_summary_metrics_txt_duplex"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_postprocessing_simplex_bam"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_reads_duplex_bam"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_collapsed_bam"
                        },
                        {
                            "id": "#bam_collapsing/gatk_collect_alignment_summary_metrics_txt_collapsed"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_group_reads_by_umi_bam"
                        },
                        {
                            "id": "#bam_collapsing/fgbio_filter_consensus_reads_simplex_duplex_bam"
                        }
                    ],
                    "run": "#bam_collapsing.cwl",
                    "label": "bam_collapsing",
                    "https://www.sevenbridges.com/x": 1871.76708984375,
                    "https://www.sevenbridges.com/y": 4696.0234375
                },
                {
                    "id": "#base_quality_recalibration",
                    "in": [
                        {
                            "id": "#base_quality_recalibration/input",
                            "source": "#uncollapsed_bam_generation/indel_realignment_bam"
                        },
                        {
                            "id": "#base_quality_recalibration/reference",
                            "source": "#reference_sequence"
                        },
                        {
                            "id": "#base_quality_recalibration/known_sites",
                            "source": [
                                "#gatk_base_recalibrator_known_sites"
                            ]
                        },
                        {
                            "id": "#base_quality_recalibration/base_recalibrator_output_file_name",
                            "source": "#base_recalibrator_output_file_name"
                        },
                        {
                            "id": "#base_quality_recalibration/add_output_sam_program_record",
                            "default": true,
                            "source": "#gatk_base_recalibrator_add_output_sam_program_record"
                        },
                        {
                            "id": "#base_quality_recalibration/apply_bqsr_create_output_bam_index",
                            "source": "#create_bam_index"
                        },
                        {
                            "id": "#base_quality_recalibration/apply_bqsr_output_file_name",
                            "source": "#apply_bqsr_output_file_name"
                        },
                        {
                            "id": "#base_quality_recalibration/temporary_directory",
                            "source": "#temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#base_quality_recalibration/gatk_apply_bqsr_bam"
                        }
                    ],
                    "run": "#base_quality_recalibration.cwl",
                    "label": "base_quality_recalibration",
                    "https://www.sevenbridges.com/x": 1871.76708984375,
                    "https://www.sevenbridges.com/y": 4162.1953125
                },
                {
                    "id": "#gatk_collect_alignment_summary_metrics_4_1_8_0",
                    "in": [
                        {
                            "id": "#gatk_collect_alignment_summary_metrics_4_1_8_0/input",
                            "source": "#base_quality_recalibration/gatk_apply_bqsr_bam"
                        },
                        {
                            "id": "#gatk_collect_alignment_summary_metrics_4_1_8_0/output_file_name",
                            "source": "#gatk_collect_aln_summary_metrics_bqsr_output_file_name"
                        },
                        {
                            "id": "#gatk_collect_alignment_summary_metrics_4_1_8_0/reference",
                            "source": "#reference_sequence"
                        },
                        {
                            "id": "#gatk_collect_alignment_summary_metrics_4_1_8_0/temporary_directory",
                            "source": "#temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#gatk_collect_alignment_summary_metrics_4_1_8_0/gatk_collect_alignment_summary_metrics_txt"
                        }
                    ],
                    "run": "#gatk_collect_alignment_summary_metrics_4.1.8.0.cwl_3",
                    "label": "GATK-CollectAlignmentSummaryMetrics",
                    "https://www.sevenbridges.com/x": 3346.1357421875,
                    "https://www.sevenbridges.com/y": 4006.0546875
                }
            ],
            "requirements": [
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Ronak Shah"
                }
            ],
            "https://schema.org/citation": "",
            "https://schema.org/codeRepository": "https://github.com/msk-access/nucleo",
            "https://schema.org/contributor": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9042-6213",
                    "https://schema.org/name": "Ronak Shah"
                }
            ],
            "https://schema.org/dateCreated": "2020-11-23",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0"
        },
        {
            "class": "Workflow",
            "id": "#base_quality_recalibration.cwl",
            "label": "base_quality_recalibration",
            "inputs": [
                {
                    "id": "#base_quality_recalibration.cwl/input",
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 533.390625
                },
                {
                    "id": "#base_quality_recalibration.cwl/reference",
                    "type": "File",
                    "secondaryFiles": [
                        ".fai",
                        "^.dict"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 106.703125
                },
                {
                    "id": "#base_quality_recalibration.cwl/read_filter",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 213.375
                },
                {
                    "id": "#base_quality_recalibration.cwl/known_sites",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "secondaryFiles": [
                        ".idx"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 426.71875
                },
                {
                    "id": "#base_quality_recalibration.cwl/base_recalibrator_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 746.734375
                },
                {
                    "id": "#base_quality_recalibration.cwl/add_output_sam_program_record",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 853.4375
                },
                {
                    "id": "#base_quality_recalibration.cwl/disable_read_filter",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 640.0625
                },
                {
                    "id": "#base_quality_recalibration.cwl/lenient",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 320.046875
                },
                {
                    "id": "#base_quality_recalibration.cwl/apply_bqsr_create_output_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 337.34375,
                    "https://www.sevenbridges.com/y": 533.453125
                },
                {
                    "id": "#base_quality_recalibration.cwl/apply_bqsr_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 337.34375,
                    "https://www.sevenbridges.com/y": 426.71875
                },
                {
                    "id": "#base_quality_recalibration.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 0
                }
            ],
            "outputs": [
                {
                    "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_bam",
                    "outputSource": [
                        "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/gatk_apply_bqsr_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 1269.836181640625,
                    "https://www.sevenbridges.com/y": 426.71875
                }
            ],
            "steps": [
                {
                    "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1",
                    "in": [
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/input",
                            "source": "#base_quality_recalibration.cwl/input"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/known_sites",
                            "source": [
                                "#base_quality_recalibration.cwl/known_sites"
                            ]
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/reference",
                            "source": "#base_quality_recalibration.cwl/reference"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/output_file_name",
                            "source": "#base_quality_recalibration.cwl/base_recalibrator_output_file_name"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/add_output_sam_program_record",
                            "source": "#base_quality_recalibration.cwl/add_output_sam_program_record"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/disable_read_filter",
                            "source": [
                                "#base_quality_recalibration.cwl/disable_read_filter"
                            ]
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/lenient",
                            "source": "#base_quality_recalibration.cwl/lenient"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/read_filter",
                            "source": [
                                "#base_quality_recalibration.cwl/read_filter"
                            ]
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/temporary_directory",
                            "source": "#base_quality_recalibration.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/gatk_base_recalibrator_output"
                        }
                    ],
                    "run": "#gatk_base_recalibrator_4.1.8.1.cwl",
                    "label": "gatk_base_recalibrator_4.1.8.1",
                    "https://www.sevenbridges.com/x": 337.34375,
                    "https://www.sevenbridges.com/y": 263.8515625
                },
                {
                    "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1",
                    "in": [
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/reference",
                            "source": "#base_quality_recalibration.cwl/reference"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/create_output_bam_index",
                            "source": "#base_quality_recalibration.cwl/apply_bqsr_create_output_bam_index"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/bqsr_recal_file",
                            "source": "#base_quality_recalibration.cwl/gatk_base_recalibrator_4_1_8_1/gatk_base_recalibrator_output"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/input",
                            "source": "#base_quality_recalibration.cwl/input"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/output_file_name",
                            "source": "#base_quality_recalibration.cwl/apply_bqsr_output_file_name"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/disable_read_filter",
                            "source": [
                                "#base_quality_recalibration.cwl/disable_read_filter"
                            ]
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/lenient",
                            "source": "#base_quality_recalibration.cwl/lenient"
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/read_filter",
                            "source": [
                                "#base_quality_recalibration.cwl/read_filter"
                            ]
                        },
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/temporary_directory",
                            "source": "#base_quality_recalibration.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#base_quality_recalibration.cwl/gatk_apply_bqsr_4_1_8_1/gatk_apply_bqsr_bam"
                        }
                    ],
                    "run": "#gatk_apply_bqsr_4.1.8.1.cwl",
                    "label": "gatk_apply_bqsr_4.1.8.1",
                    "https://www.sevenbridges.com/x": 837.3018188476562,
                    "https://www.sevenbridges.com/y": 370.5859375
                }
            ],
            "requirements": [],
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9042-6213",
                    "https://schema.org/name": "Ronak Shah"
                }
            ],
            "https://schema.org/citation": "",
            "https://schema.org/codeRepository": "https://github.com/msk-access/cwl_subworkflows/base_quality_recalibration",
            "https://schema.org/contributor": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9042-6213",
                    "https://schema.org/name": "Ronak Shah"
                }
            ],
            "https://schema.org/dateCreated": "2020-06-09",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0"
        },
        {
            "class": "CommandLineTool",
            "id": "#gatk_apply_bqsr_4.1.8.1.cwl",
            "baseCommand": [
                "gatk",
                "ApplyBQSR"
            ],
            "inputs": [
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/reference",
                    "type": "File",
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--reference"
                    },
                    "doc": "Reference sequence",
                    "secondaryFiles": [
                        ".fai",
                        "^.dict"
                    ]
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/create_output_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--create-output-bam-index"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/bqsr_recal_file",
                    "type": "File",
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--bqsr-recal-file"
                    },
                    "doc": "Input recalibration table for BQSR. Only run ApplyBQSR with the covariates table created from the input BAM"
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--input"
                    },
                    "doc": "A BAM file containing input read data",
                    "secondaryFiles": [
                        "^.bai"
                    ]
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output file name. Not Required"
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/add_output_sam_program_record",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--add-output-sam-program-record"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/add_output_vcf_command_line",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--add-output-vcf-command-line"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/arguments_file",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--arguments_file"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/cloud_index_prefetch_buffer",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--cloud-index-prefetch-buffer"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/cloud_prefetch_buffer",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--cloud-prefetch-buffer"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/create_output_bam_md5",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--create-output-bam-md5"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/create_output_variant_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--create-output-variant-index"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/create_output_variant_md5",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--create-output-variant-md5"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/disable_bam_index_caching",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--disable-bam-index-caching"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/disable_read_filter",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string",
                            "inputBinding": {
                                "prefix": "--disable-read-filter"
                            }
                        }
                    ],
                    "inputBinding": {
                        "position": 6
                    },
                    "doc": "Read filters to be disabled before analysis"
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/disable_sequence_dictionary_validation",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--disable-sequence-dictionary-validation"
                    }
                },
                {
                    "default": true,
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/emit_original_quals",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--emit-original-quals"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/exclude_intervals",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--exclude-intervals"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/gatk_config_file",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--gatk-config-file"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/gcs_max_retries",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--gcs-max-retries"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/gcs_project_for_requester_pays",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--gcs-project-for-requester-pays"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/global_qscore_prior",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--global-qscore-prior"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/interval_exclusion_padding",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--interval-exclusion-padding"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/interval_merging_rule",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--interval-merging-rule"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/interval_padding",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--interval-padding"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/interval_set_rule",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--interval-set-rule"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/intervals",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--intervals"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/lenient",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--lenient"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/preserve_qscores_less_than",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--preserve-qscores-less-than"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/quantize_quals",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--quantize-quals"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/quiet",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--QUIET"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/read_filter",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string",
                            "inputBinding": {
                                "prefix": "--read-filter"
                            }
                        }
                    ],
                    "inputBinding": {
                        "position": 6
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/read_index",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--read-index"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/read_validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--read-validation-stringency"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/seconds_between_progress_updates",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--seconds-between-progress-updates"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/sequence_dictionary",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--sequence-dictionary"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/sites_only_vcf_output",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--sites-only-vcf-output"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--use-jdk-deflater"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--use-jdk-inflater"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/use_original_qualities",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--use-original-qualities"
                    }
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#gatk_apply_bqsr_4.1.8.1.cwl/gatk_apply_bqsr_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/, '_bqsr.bam')\n    }\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "label": "gatk_apply_bqsr_4.1.8.1",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--java-options",
                    "valueFrom": "${\n     if(inputs.memory_per_job && inputs.memory_overhead){\n        if(inputs.memory_per_job % 1000 == 0){\n            return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n        } else {\n            return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n        }\n     } else if (inputs.memory_per_job && !inputs.memory_overhead){\n        if(inputs.memory_per_job % 1000 == 0) {\n            return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n        } else {\n            return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n        }\n     } else if(!inputs.memory_per_job && inputs.memory_overhead){\n        return \"-Xmx12G\"\n     } else {\n        return \"-Xmx12G\"\n     }\n}"
                },
                {
                    "position": 2,
                    "prefix": "--tmp-dir",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 2,
                    "prefix": "--output",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/, '_bqsr.bam')\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 16000,
                    "coresMin": 4
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.1"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                },
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:sumans@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Shalabh Suman"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "gatk4",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.1"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#gatk_base_recalibrator_4.1.8.1.cwl",
            "baseCommand": [
                "gatk",
                "BaseRecalibrator"
            ],
            "inputs": [
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--input"
                    },
                    "doc": "BAM/SAM file containing reads",
                    "secondaryFiles": [
                        "^.bai"
                    ]
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/known_sites",
                    "type": {
                        "type": "array",
                        "items": "File",
                        "inputBinding": {
                            "prefix": "--known-sites"
                        }
                    },
                    "inputBinding": {
                        "position": 3
                    },
                    "doc": "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis",
                    "secondaryFiles": [
                        ".idx"
                    ]
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/reference",
                    "type": "File",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--reference"
                    },
                    "doc": "Reference sequence file",
                    "secondaryFiles": [
                        ".fai",
                        "^.dict"
                    ]
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output file name. Not Required"
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/add_output_sam_program_record",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--add-output-sam-program-record"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/add_output_vcf_command_line",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--add-output-vcf-command-line"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/arguments_file",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File",
                            "inputBinding": {
                                "position": 0,
                                "prefix": "--arguments_file"
                            }
                        }
                    ]
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/binary_tag_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--binary-tag-name"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/bqsr_baq_gap_open_penalty",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--bqsr-baq-gap-open-penalty"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/cloud-index-prefetch-buffer",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--cloud-index-prefetch-buffer"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/cloud_prefetch_buffer",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--cloud-prefetch-buffer"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/create_output_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--create-output-bam-index"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/create_output_bam_md5",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--create-output-bam-md5"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/create_output_variant_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--create-output-variant-index"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/create_output_variant_md5",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--create-output-variant-md5"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/default_base_qualities",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--default-base-qualities"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/deletions_default_quality",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--deletions-default-quality"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/disable_bam_index_caching",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--disable-bam-index-caching"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/disable_read_filter",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string",
                            "inputBinding": {
                                "prefix": "--disable-read-filter"
                            }
                        }
                    ],
                    "inputBinding": {
                        "position": 10
                    },
                    "doc": "Read filters to be disabled before analysis"
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/disable_sequence_dictionary_validation",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--disable-sequence-dictionary-validation"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/exclude_intervals",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--exclude-intervals"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/gatk_config_file",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--gatk-config-file"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/gcs_max_retries",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--gcs-max-retries"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/gcs_project_for_requester_pays",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--gcs-project-for-requester-pays"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/indels_context_size",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--indels-context-size"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/insertions_default_quality",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--insertions-default-quality"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/interval_exclusion_padding",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--interval-exclusion-padding"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/interval_merging_rule",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--interval-merging-rule"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/interval_padding",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--interval-padding"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/interval_set_rule",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--interval-set-rule"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/intervals",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--intervals"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/lenient",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--lenient"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/low_quality_tail",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--low-quality-tail"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/maximum_cycle_value",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--maximum-cycle-value"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/mismatches_context_size",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--mismatches-context-size"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/mismatches_default_quality",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--mismatches-default-quality"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/preserve_qscores_less_than",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--preserve-qscores-less-than"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/quantizing_levels",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--quantizing-levels"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/QUIET",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--QUIET"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/read_filter",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string",
                            "inputBinding": {
                                "prefix": "--read-filter"
                            }
                        }
                    ],
                    "inputBinding": {
                        "position": 10
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/read_index",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--read-index"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/seconds_between_progress_updates",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--seconds-between-progress-updates"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/sequence_dictionary",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--sequence-dictionary"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/sites_only_vcf_output",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--sites-only-vcf-output"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/use_original_qualities",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--use-original-qualities"
                    }
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#gatk_base_recalibrator_4.1.8.1.cwl/gatk_base_recalibrator_output",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/, '_bqsr.table')\n    }\n}"
                    }
                }
            ],
            "label": "gatk_base_recalibrator_4.1.8.1",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--java-options",
                    "valueFrom": "${\n     if(inputs.memory_per_job && inputs.memory_overhead){\n        if(inputs.memory_per_job % 1000 == 0){\n            return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n        } else {\n            return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n        }\n     } else if (inputs.memory_per_job && !inputs.memory_overhead){\n        if(inputs.memory_per_job % 1000 == 0) {\n            return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n        } else {\n            return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n        }\n     } else if(!inputs.memory_per_job && inputs.memory_overhead){\n        return \"-Xmx12G\"\n     } else {\n        return \"-Xmx12G\"\n     }\n}"
                },
                {
                    "position": 2,
                    "prefix": "--tmp-dir",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 2,
                    "prefix": "--output",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/, '_bqsr.table')\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 32000,
                    "coresMin": 8
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.1"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                },
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:sumans@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Shalabh Suman"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#fastp_0.20.1.cwl",
            "baseCommand": [
                "fastp"
            ],
            "inputs": [
                {
                    "id": "#fastp_0.20.1.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#fastp_0.20.1.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#fastp_0.20.1.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "worker thread number, default is 2 (int [=2])"
                },
                {
                    "id": "#fastp_0.20.1.cwl/read1_input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--in1"
                    },
                    "doc": "read1 input file name\n"
                },
                {
                    "id": "#fastp_0.20.1.cwl/read1_output_path",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--out1"
                    },
                    "doc": "read1 output file name\n"
                },
                {
                    "id": "#fastp_0.20.1.cwl/read2_input",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--in2"
                    },
                    "doc": "read2 input file name, for PE data\n"
                },
                {
                    "id": "#fastp_0.20.1.cwl/read2_output_path",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--out2"
                    },
                    "doc": "read2 output file name\n"
                },
                {
                    "id": "#fastp_0.20.1.cwl/unpaired1_path",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--unpaired1"
                    },
                    "doc": "for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. \n"
                },
                {
                    "id": "#fastp_0.20.1.cwl/unpaired2_path",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--unpaired2"
                    },
                    "doc": "for PE input, if read2 passed QC but read1 not, it will be written to unpaired2.\n"
                },
                {
                    "id": "#fastp_0.20.1.cwl/failed_reads_path",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--failed_out"
                    },
                    "doc": "specify the file to store reads that cannot pass the filters.\n"
                },
                {
                    "id": "#fastp_0.20.1.cwl/read1_adapter_sequence",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--adapter_sequence"
                    },
                    "doc": "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped.\n"
                },
                {
                    "id": "#fastp_0.20.1.cwl/read2_adapter_sequence",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--adapter_sequence_r2"
                    },
                    "doc": "the adapter for read2. For PE data, this is used if R1/R2 are found not overlapped.\n"
                },
                {
                    "id": "#fastp_0.20.1.cwl/minimum_read_length",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--length_required"
                    },
                    "doc": "reads shorter than length_required will be discarded, default is 15.\n"
                },
                {
                    "default": "fastp.json",
                    "id": "#fastp_0.20.1.cwl/json_output_path",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--json"
                    },
                    "doc": "the json format report file name\n"
                },
                {
                    "default": "fastp.html",
                    "id": "#fastp_0.20.1.cwl/html_output_path",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--html"
                    },
                    "doc": "the html format report file name\n"
                },
                {
                    "id": "#fastp_0.20.1.cwl/disable_quality_filtering",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--disable_quality_filtering"
                    },
                    "doc": "quality filtering is enabled by default. If this option is specified, quality filtering is disabled"
                },
                {
                    "id": "#fastp_0.20.1.cwl/disable_trim_poly_g",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--disable_trim_poly_g"
                    },
                    "doc": "disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data"
                },
                {
                    "id": "#fastp_0.20.1.cwl/verbose",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--verbose"
                    },
                    "doc": "output verbose log information (i.e. when every 1M reads are processed)"
                }
            ],
            "outputs": [
                {
                    "id": "#fastp_0.20.1.cwl/fastp_json_output",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.json_output_path)"
                    }
                },
                {
                    "id": "#fastp_0.20.1.cwl/fastp_html_output",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.html_output_path)"
                    }
                },
                {
                    "id": "#fastp_0.20.1.cwl/fastp_read1_output",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.read1_output_path)"
                    }
                },
                {
                    "id": "#fastp_0.20.1.cwl/fastp_read2_output",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.read2_output_path)"
                    }
                },
                {
                    "id": "#fastp_0.20.1.cwl/fastp_unpaired1_output",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.unpaired1_path)"
                    }
                },
                {
                    "id": "#fastp_0.20.1.cwl/fastp_unpaired2_output",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.unpaired2_path)"
                    }
                }
            ],
            "doc": "Setup and execute Fastp",
            "label": "fastp_0.20.1",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--thread",
                    "valueFrom": "${\n    if(inputs.number_of_threads)\n        return inputs.number_of_threads\n    return runtime.cores\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 17000,
                    "coresMin": 4
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/fastp:0.20.1--h8b12597_0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charlie Murphy"
                        },
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:fraihaa@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Adrian Fraiha"
                        },
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "fastp",
                    "http://usefulinc.com/ns/doap#revision": "0.20.1"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#fgbio_fastq_to_bam_1.2.0.cwl",
            "baseCommand": [
                "fgbio"
            ],
            "inputs": [
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/input",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--input",
                        "itemSeparator": " ",
                        "shellQuote": false
                    },
                    "label": "PathToFastq",
                    "doc": "Fastq files corresponding to each sequencing read (e.g. R1, I1, etc.)."
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "The output SAM or BAM file to be written."
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/read-structures",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--read-structures",
                        "itemSeparator": " ",
                        "shellQuote": false
                    },
                    "doc": "Read structures, one for each of the FASTQs. https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/sort",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--sort",
                        "shellQuote": false
                    },
                    "doc": "If true, queryname sort the BAM file, otherwise preserve input order."
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/umi-tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--umi-tag",
                        "shellQuote": false
                    },
                    "doc": "Tag in which to store molecular barcodes/UMIs"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/read-group-id",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--read-group-id",
                        "shellQuote": false
                    },
                    "doc": "Read group ID to use in the file header."
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/sample",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--sample",
                        "shellQuote": false
                    },
                    "doc": "The name of the sequenced sample."
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/library",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--library",
                        "shellQuote": false
                    },
                    "doc": "The name/ID of the sequenced library."
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/platform",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--platform",
                        "shellQuote": false
                    },
                    "doc": "Sequencing Platform"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/platform-unit",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--platform-unit",
                        "shellQuote": false
                    },
                    "doc": "Platform unit (e.g. \u2018..')"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/platform-model",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--platform-model",
                        "shellQuote": false
                    },
                    "doc": "Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/sequencing-center",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--sequencing-center",
                        "shellQuote": false
                    },
                    "doc": "The sequencing center from which the data originated"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/predicted-insert-size",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--predicted-insert-size",
                        "shellQuote": false
                    },
                    "doc": "Predicted median insert size, to insert into the read group header"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/description",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--description"
                    },
                    "doc": "Description of the read group."
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/comment",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--comment"
                    },
                    "doc": "Comment(s) to include in the output file\u2019s header"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/run-date",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--run-date",
                        "shellQuote": false
                    },
                    "doc": "Date the run was produced, to insert into the read group header"
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null."
                },
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/async_io",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "separate": false,
                        "prefix": "--async-io="
                    },
                    "doc": "'Use asynchronous I/O where possible, e.g. for SAM and BAM files [=true|false].'"
                }
            ],
            "outputs": [
                {
                    "id": "#fgbio_fastq_to_bam_1.2.0.cwl/fgbio_fastq_to_bam_ubam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n    return  inputs.input[0].basename.replace(/.fastq.gz/,'_ubam.bam');\n}"
                    }
                }
            ],
            "doc": "Generates an unmapped BAM (or SAM or CRAM) file from fastq files. Takes in one or more fastq files (optionally gzipped), each representing a different sequencing read (e.g. R1, R2, I1 or I2) and can use a set of read structures to allocate bases in those reads to template reads, sample indices, unique molecular indices, or to designate bases to be skipped over.\n\nRead structures are made up of <number><operator> pairs much like the CIGAR string in BAM files. Four kinds of operators are recognized:\n\n1. T identifies a template read\n2. B identifies a sample barcode read\n3. M identifies a unique molecular index read\n4. S identifies a set of bases that should be skipped or ignored\n\nThe last <number><operator> pair may be specified using a + sign instead of number to denote \u201call remaining bases\u201d. This is useful if, e.g., fastqs have been trimmed and contain reads of varying length. For example to convert a paired-end run with an index read and where the first 5 bases of R1 are a UMI and the second five bases are monotemplate you might specify:",
            "label": "fgbio_fastq_to_bam_1.2.0",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx12G\"\n  }\n  else {\n      return \"-Xmx12G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "valueFrom": "-XX:-UseGCOverheadLimit"
                },
                {
                    "position": 1,
                    "valueFrom": "FastqToBam"
                },
                {
                    "position": 0,
                    "prefix": "--tmp-dir=",
                    "separate": false,
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 2,
                    "prefix": "--output",
                    "shellQuote": false,
                    "valueFrom": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n      return  inputs.input[0].basename.replace(/.fastq.gz/,'_ubam.bam');\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 16000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/fgbio:1.2.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "fgbio FastqToBam",
                    "http://usefulinc.com/ns/doap#revision": "1.2.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2",
            "baseCommand": [
                "gatk",
                "MergeBamAlignment"
            ],
            "inputs": [
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/unmapped_bam",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--UNMAPPED_BAM"
                    },
                    "doc": "Original SAM or BAM file of unmapped reads, which must be in queryname order.  Reads MUST\nbe unmapped. Required.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/reference",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--REFERENCE_SEQUENCE"
                    },
                    "doc": "Reference sequence file.  Required.\n",
                    "secondaryFiles": [
                        "^.dict"
                    ]
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Merged SAM or BAM file to write to.  Required.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/add_mate_cigar",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ADD_MATE_CIGAR"
                    },
                    "doc": "Adds the mate CIGAR tag (MC) if true, does not if false.  Default value: true. Possible\nvalues: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/add_pg_tag_to_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ADD_PG_TAG_TO_READS"
                    },
                    "doc": "Add PG tag to each read in a SAM or BAM  Default value: true. Possible values: {true,\nfalse}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/aligned_bam",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File",
                            "inputBinding": {
                                "prefix": "--ALIGNED_BAM"
                            }
                        }
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "doc": "SAM or BAM file(s) with alignment data.  This argument may be specified 0 or more times.\nDefault value: null.  Cannot be used in conjunction with argument(s) READ1_ALIGNED_BAM\n(R1_ALIGNED) READ2_ALIGNED_BAM (R2_ALIGNED)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/aligned_reads_only",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ALIGNED_READS_ONLY"
                    },
                    "doc": "Whether to output only aligned reads. Default value: false. Possible values: {true,\nfalse}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/aligner_proper_pair_flags",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ALIGNER_PROPER_PAIR_FLAGS"
                    },
                    "doc": "Use the aligners idea of what a proper pair is rather than computing in this program.\nDefault value: false. Possible values: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/attributes_to_remove",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ATTRIBUTES_TO_REMOVE"
                    },
                    "doc": "Attributes from the alignment record that should be removed when merging.  This overrides\nATTRIBUTES_TO_RETAIN if they share common tags.  This argument may be specified 0 or more\ntimes. Default value: null.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/attributes_to_retain",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ATTRIBUTES_TO_RETAIN"
                    },
                    "doc": "Reserved alignment attributes (tags starting with X, Y, or Z) that should be brought over\nfrom the alignment data when merging.  This argument may be specified 0 or more times.\nDefault value: null.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/attributes_to_reverse",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ATTRIBUTES_TO_REVERSE"
                    },
                    "doc": "Attributes on negative strand reads that need to be reversed.  This argument may be\nspecified 0 or more times. Default value: [OQ, U2].\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/attributes_to_reverse_complement",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ATTRIBUTES_TO_REVERSE_COMPLEMENT"
                    },
                    "doc": "Attributes on negative strand reads that need to be reverse complemented.  This argument\nmay be specified 0 or more times. Default value: [E2, SQ].\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/clip_adapters",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CLIP_ADAPTERS"
                    },
                    "doc": "Whether to clip adapters where identified.  Default value: true. Possible values: {true,\nfalse}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/clip_overlapping_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CLIP_OVERLAPPING_READS"
                    },
                    "doc": "For paired reads, clip the 3' end of each read if necessary so that it does not extend\npast the 5' end of its mate.  Clipping will be either soft or hard clipping, depending on\nCLIP_OVERLAPPING_READS_OPERATOR setting. Hard clipped bases and their qualities will be\nstored in the XB and XQ tags respectively.  Default value: true. Possible values: {true,\nfalse}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/expected_orientations",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--EXPECTED_ORIENTATIONS"
                    },
                    "doc": "The expected orientation of proper read pairs. Replaces JUMP_SIZE  This argument may be\nspecified 0 or more times. Default value: null. Possible values: {FR, RF, TANDEM}  Cannot\nbe used in conjunction with argument(s) JUMP_SIZE (JUMP)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/hard_clip_overlapping_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--HARD_CLIP_OVERLAPPING_READS"
                    },
                    "doc": "If true, hard clipping will be applied to overlapping reads.  By default, soft clipping is\nused.  Default value: false. Possible values: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/include_secondary_alignments",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--INCLUDE_SECONDARY_ALIGNMENTS"
                    },
                    "doc": "If false, do not write secondary alignments to output.  Default value: true. Possible\nvalues: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/is_bisulfite_sequence",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--IS_BISULFITE_SEQUENCE"
                    },
                    "doc": "Whether the lane is bisulfite sequence (used when calculating the NM tag).  Default value:\nfalse. Possible values: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/jump_size",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--JUMP_SIZE"
                    },
                    "doc": "The expected jump size (required if this is a jumping library). Deprecated. Use\nEXPECTED_ORIENTATIONS instead  Default value: null.  Cannot be used in conjunction with\nargument(s) EXPECTED_ORIENTATIONS (ORIENTATIONS)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/matching_dictionary_tags",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MATCHING_DICTIONARY_TAGS"
                    },
                    "doc": "List of Sequence Records tags that must be equal (if present) in the reference dictionary\nand in the aligned file. Mismatching tags will cause an error if in this list, and a\nwarning otherwise.  This argument may be specified 0 or more times. Default value: [M5,\nLN].\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/max_insertions_or_deletions",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MAX_INSERTIONS_OR_DELETIONS"
                    },
                    "doc": "The maximum number of insertions or deletions permitted for an alignment to be included.\nAlignments with more than this many insertions or deletions will be ignored. Set to -1 to\nallow any number of insertions or deletions.  Default value: 1.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/min_unclipped_bases",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MIN_UNCLIPPED_BASES"
                    },
                    "doc": "If UNMAP_CONTAMINANT_READS is set, require this many unclipped bases or else the read will\nbe marked as contaminant.  Default value: 32.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/paired_run",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--PAIRED_RUN"
                    },
                    "doc": "DEPRECATED. This argument is ignored and will be removed.  Default value: true. Possible\nvalues: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/primary_alignment_strategy",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--PRIMARY_ALIGNMENT_STRATEGY"
                    },
                    "doc": "Strategy for selecting primary alignment when the aligner has provided more than one\nalignment for a pair or fragment, and none are marked as primary, more than one is marked\nas primary, or the primary alignment is filtered out for some reason. For all strategies,\nties are resolved arbitrarily.  Default value: BestMapq. BestMapq (Expects that multiple\nalignments will be correlated with HI tag, and prefers the pair of alignments with the\nlargest MAPQ, in the absence of a primary selected by the aligner.)\nEarliestFragment (Prefers the alignment which maps the earliest base in the read. Note\nthat EarliestFragment may not be used for paired reads.)\nBestEndMapq (Appropriate for cases in which the aligner is not pair-aware, and does not\noutput the HI tag. It simply picks the alignment for each end with the highest MAPQ, and\nmakes those alignments primary, regardless of whether the two alignments make sense\ntogether.)\nMostDistant (Appropriate for a non-pair-aware aligner. Picks the alignment pair with the\nlargest insert size. If all alignments would be chimeric, it picks the alignments for each\nend with the best MAPQ.)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/read1_aligned_bam",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File",
                            "inputBinding": {
                                "prefix": "--READ1_ALIGNED_BAM"
                            }
                        }
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "doc": "SAM or BAM file(s) with alignment data from the first read of a pair.  This argument may\nbe specified 0 or more times. Default value: null.  Cannot be used in conjunction with\nargument(s) ALIGNED_BAM (ALIGNED)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/read1_trim",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ1_TRIM"
                    },
                    "doc": "The number of bases trimmed from the beginning of read 1 prior to alignment  Default\nvalue: 0.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/read2_aligned_bam",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File",
                            "inputBinding": {
                                "prefix": "--READ2_ALIGNED_BAM"
                            }
                        }
                    ],
                    "inputBinding": {
                        "position": 1
                    },
                    "doc": "SAM or BAM file(s) with alignment data from the second read of a pair.  This argument may\nbe specified 0 or more times. Default value: null.  Cannot be used in conjunction with\nargument(s) ALIGNED_BAM (ALIGNED)\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/read2_trim",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ2_TRIM"
                    },
                    "doc": "The number of bases trimmed from the beginning of read 2 prior to alignment  Default\nvalue: 0.\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--SORT_ORDER"
                    },
                    "doc": "The order in which the merged reads should be output.  Default value: coordinate. Possible\nvalues: {unsorted, queryname, coordinate, duplicate, unknown}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/unmap_contaminant_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--UNMAP_CONTAMINANT_READS"
                    },
                    "doc": "Detect reads originating from foreign organisms (e.g. bacterial DNA in a non-bacterial\nsample),and unmap + label those reads accordingly.  Default value: false. Possible values:\n{true, false}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/unmapped_read_strategy",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--UNMAPPED_READ_STRATEGY"
                    },
                    "doc": "How to deal with alignment information in reads that are being unmapped (e.g. due to\ncross-species contamination.) Currently ignored unless UNMAP_CONTAMINANT_READS = true.\nNote that the DO_NOT_CHANGE strategy will actually reset the cigar and set the mapping\nquality on unmapped reads since otherwisethe result will be an invalid record. To force no\nchange use the DO_NOT_CHANGE_INVALID strategy.  Default value: DO_NOT_CHANGE. Possible\nvalues: {COPY_TO_TAG, DO_NOT_CHANGE, DO_NOT_CHANGE_INVALID, MOVE_TO_TAG}\n"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT,LENIENT, SILENT}"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/create_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/create_md5_file",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_MD5_FILE"
                    },
                    "doc": "Whether to create an MD5 digest for any BAM or FASTQ files created.    Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2/gatk_merge_bam_alignment_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.unmapped_bam.basename.replace(/.bam|.sam/, '_merged.bam')\n    }\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "label": "GATK-MergeBamAlignment",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--java-options",
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx15G\"\n  }\n  else {\n      return \"-Xmx15G\"\n  }\n}"
                },
                {
                    "position": 1,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.unmapped_bam.basename.replace(/.bam|.sam/, '_merged.bam')\n    }\n}"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 17000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "gatk4",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#gatk_merge_sam_files_4.1.8.0.cwl",
            "baseCommand": [
                "gatk",
                "MergeSamFiles"
            ],
            "inputs": [
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/input",
                    "type": {
                        "type": "array",
                        "items": "File",
                        "inputBinding": {
                            "prefix": "-I"
                        }
                    },
                    "inputBinding": {
                        "position": 1
                    },
                    "doc": "SAM or BAM input file  This argument must be specified at least once. Required.\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "SAM or BAM file to write merged result to  Required."
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/assume_sorted",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--ASSUME_SORTED"
                    },
                    "doc": "If true, assume that the input files are in the same sort order as the requested output\nsort order, even if their headers say otherwise.  Default value: false. Possible values:\n{true, false}\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/comment",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--COMMENT"
                    },
                    "doc": "Comment(s) to include in the merged output files header.  This argument may be specified\n0 or more times. Default value: null.\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/create_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value:\nfalse. Possible values: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/create_md5_file",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--CREATE_MD5_FILE"
                    },
                    "doc": "Whether to create an MD5 digest for any BAM or FASTQ files created.    Default value:\nfalse. Possible values: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/intervals",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--INTERVALS"
                    },
                    "doc": "An interval list file that contains the locations of the positions to merge. Assume bam\nare sorted and indexed. The resulting file will contain alignments that may overlap with\ngenomic regions outside the requested region. Unmapped reads are discarded.  Default\nvalue: null.\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/merge_sequence_dictionaries",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--MERGE_SEQUENCE_DICTIONARIES"
                    },
                    "doc": "Merge the sequence dictionaries  Default value: false. Possible values: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/reference_sequence",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--REFERENCE_SEQUENCE"
                    },
                    "doc": "Reference sequence file.  Default value: null.\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--SORT_ORDER"
                    },
                    "doc": "Sort order of output file  Default value: coordinate. Possible values: {unsorted,\nqueryname, coordinate, duplicate, unknown}\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/use_threading",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--USE_THREADING"
                    },
                    "doc": "Option to create a background thread to encode, compress and write to disk the output\nfile. The threaded version uses about 20% more CPU and decreases runtime by ~20% when\nwriting out a compressed BAM file.  Default value: false. Possible values: {true, false}\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to\nSILENT can improve performance when processing a BAM file in which variable-length data\n(read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT.\nPossible values: {STRICT, LENIENT, SILENT}\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/verbosity",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--VERBOSITY"
                    },
                    "doc": "Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING,\nINFO, DEBUG}\n"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#gatk_merge_sam_files_4.1.8.0.cwl/gatk_merge_sam_files_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return 'merged.bam'\n    }\n}"
                    }
                }
            ],
            "label": "GATK-MergeSamFiles",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--java-options",
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx15G\"\n  }\n  else {\n      return \"-Xmx15G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 2,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return 'merged.bam'\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 17000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "gatk4",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2",
            "baseCommand": [
                "gatk",
                "SamToFastq"
            ],
            "inputs": [
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/fastq",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output FASTQ file (single-end fastq or, if paired, first end of the pair FASTQ)"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--INPUT"
                    },
                    "doc": "Input SAM/BAM file to extract reads from  Required."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/clipping_action",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CLIPPING_ACTION"
                    },
                    "doc": "The action that should be taken with clipped reads: 'X' means the reads and qualities should be trimmed at the clipped position; 'N' means the bases should be changed to Ns in the clipped region; and any integer means that the base qualities should be set to that value in the clipped region.  Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/clipping_attribute",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CLIPPING_ATTRIBUTE"
                    },
                    "doc": "The attribute that stores the position at which the SAM record should be clipped  Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/clipping_min_length",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CLIPPING_MIN_LENGTH"
                    },
                    "doc": "When performing clipping with the CLIPPING_ATTRIBUTE and CLIPPING_ACTION parameters, ensure that the resulting reads after clipping are at least CLIPPING_MIN_LENGTH bases long. If the original read is shorter than CLIPPING_MIN_LENGTH then the original read length will be maintained.  Default value: 0."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/compress_outputs_per_rg",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--COMPRESS_OUTPUTS_PER_RG"
                    },
                    "doc": "Compress output FASTQ files per read group using gzip and append a .gz extension to the file names.  Default value: false. Possible values: {true, false}  Cannot be used in conjunction with argument(s) FASTQ (F) SECOND_END_FASTQ (F2) UNPAIRED_FASTQ (FU)"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/compression_level",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--COMPRESSION_LEVEL"
                    },
                    "doc": "Compression level for all compressed files created (e.g. BAM and VCF).  Default value: 2."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/create_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/include_non_pf_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--INCLUDE_NON_PF_READS"
                    },
                    "doc": "Include non-PF reads from the SAM file into the output FASTQ files. PF means 'passes filtering'. Reads whose 'not passing quality controls' flag is set are non-PF reads. See GATK Dictionary for more info.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/include_non_primary_alignments",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--INCLUDE_NON_PRIMARY_ALIGNMENTS"
                    },
                    "doc": "If true, include non-primary alignments in the output.  Support of non-primary alignments in SamToFastq is not comprehensive, so there may be exceptions if this is set to true and there are paired reads with non-primary alignments.  Default value: false. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/interleave",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--INTERLEAVE"
                    },
                    "doc": "Will generate an interleaved fastq if paired, each line will have /1 or /2 to describe which end it came from  Default value: false. Possible values: {true, false}"
                },
                {
                    "default": 50000,
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/max_records_in_ram",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--MAX_RECORDS_IN_RAM"
                    },
                    "doc": "When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.  Default value: 500000."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/output_dir",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--OUTPUT_DIR"
                    },
                    "doc": "Directory in which to output the FASTQ file(s). Used only when OUTPUT_PER_RG is true. Default value: null. Cannot be used in conjunction with argument(s) FASTQ (F)."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/create_md5_file",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_MD5_FILE"
                    },
                    "doc": "Whether to create an MD5 digest for any BAM or FASTQ files created.    Default value: false. Possible values: {true, false}."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/output_per_rg",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--OUTPUT_PER_RG"
                    },
                    "doc": "Output a FASTQ file per read group (two FASTQ files per read group if the group is paired).  Default value: false. Possible values: {true, false}  Cannot be used in conjunction with argument(s) FASTQ (F) SECOND_END_FASTQ (F2) UNPAIRED_FASTQ (FU)"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/quality",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--QUALITY"
                    },
                    "doc": "End-trim reads using the phred/bwa quality trimming algorithm and this quality. Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/re_reverse",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RE_REVERSE"
                    },
                    "doc": "Re-reverse bases and qualities of reads with negative strand flag set before writing them to FASTQ  Default value: true. Possible values: {true, false}"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/read1_max_bases_to_write",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ1_MAX_BASES_TO_WRITE"
                    },
                    "doc": "The maximum number of bases to write from read 1 after trimming. If there are fewer than this many bases left after trimming, all will be written.  If this value is null then all bases left after trimming will be written.  Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/read1_trim",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ1_TRIM"
                    },
                    "doc": "The number of bases to trim from the beginning of read 1.  Default value: 0."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/read2_max_bases_to_write",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ2_MAX_BASES_TO_WRITE"
                    },
                    "doc": "The maximum number of bases to write from read 2 after trimming. If there are fewer than this many bases left after trimming, all will be written.  If this value is null then all bases left after trimming will be written.  Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/read2_trim",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ2_TRIM"
                    },
                    "doc": "The number of bases to trim from the beginning of read 2.  Default value: 0."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/reference_sequence",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--REFERENCE_SEQUENCE"
                    },
                    "doc": "Reference sequence file. Default value: null."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/rg_tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RG_TAG"
                    },
                    "doc": "The read group tag (PU or ID) to be used to output a FASTQ file per read group.  Default value: PU."
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/second_end_fastq",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--SECOND_END_FASTQ"
                    },
                    "doc": "Output FASTQ file (if paired, second end of the pair FASTQ).  Default value: null.  Cannot be used in conjunction with argument(s) OUTPUT_PER_RG (OPRG) COMPRESS_OUTPUTS_PER_RG (GZOPRG)"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/unpaired_fastq",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--UNPAIRED_FASTQ"
                    },
                    "doc": "Output FASTQ file for unpaired reads; may only be provided in paired-FASTQ mode  Default value: null.  Cannot be used in conjunction with argument(s) OUTPUT_PER_RG (OPRG) COMPRESS_OUTPUTS_PER_RG (GZOPRG)"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT}"
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/gatk_sam_to_fastq_fastq",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.fastq){\n      return inputs.fastq\n    } else {\n      return inputs.input.basename.replace(/.bam|.sam/, '_R1.fastq')\n    }\n}"
                    }
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/gatk_sam_to_fastq_unpaired_fastq",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "${\n    if(inputs.unpaired_fastq){\n        return inputs.unpaired_fastq\n    } else {\n      return inputs.input.basename.replace(/.bam|.sam/, '_unpaired.fastq')\n    }\n}"
                    }
                },
                {
                    "id": "#gatk_sam_to_fastq_4.1.8.0.cwl_2/gatk_sam_to_fastq_second_end_fastq",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "${\n    if(inputs.second_end_fastq){\n        return inputs.second_end_fastq\n    } else {\n      return inputs.input.basename.replace(/.bam|.sam/, '_R2.fastq')\n    }\n}"
                    }
                }
            ],
            "label": "GATK-SamToFastq",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "--java-options",
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx15G\"\n  }\n  else {\n      return \"-Xmx15G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "prefix": "--FASTQ",
                    "valueFrom": "${\n    if(inputs.fastq){\n        return inputs.fastq\n    } else {\n        return inputs.input.basename.replace(/.bam|.sam/, '_R1.fastq')\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 8000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.0"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:murphyc4@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Charles Murphy"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "gatk4",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.0"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#picard_mark_duplicates_4.1.8.1.cwl",
            "baseCommand": [
                "java"
            ],
            "inputs": [
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-I"
                    },
                    "doc": "Input file (bam or sam).  Required."
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output file (bam or sam)."
                },
                {
                    "default": "$( inputs.input.basename.replace(/.bam/, '_md.metrics') )",
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/duplication_metrics",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-M"
                    },
                    "doc": "File to write duplication metrics to Required."
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/assume_sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-ASO"
                    },
                    "doc": "Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.Default value: null. Possible values: {unsorted, queryname, coordinate}"
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/tmp_dir",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--TMP_DIR"
                    },
                    "doc": "This option may be specified 0 or more times"
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT,LENIENT, SILENT}"
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/bam_compression_level",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--COMPRESSION_LEVEL"
                    },
                    "doc": "Compression level for all compressed files created (e.g. BAM and GELI). Default value:5. This option can be set to 'null' to clear the default value."
                },
                {
                    "default": true,
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value:false. This option can be set to 'null' to clear the default value. Possible values:{true, false}"
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/read_name_regex",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--READ_NAME_REGEX"
                    },
                    "doc": "MarkDuplicates can use the tile and cluster positions to estimate the rate of optical duplication in addition to the dominant source of duplication, PCR, to provide a more accurate estimation of library size. By default (with no READ_NAME_REGEX specified), MarkDuplicates will attempt to extract coordinates using a split on ':' (see Note below). Set READ_NAME_REGEX to 'null' to disable optical duplicate detection. Note that without optical duplicate counts, library size estimation will be less accurate. If the read name does not follow a standard Illumina colon-separation convention, but does contain tile and x,y coordinates, a regular expression can be specified to extract three variables: tile/region, x coordinate and y coordinate from a read name. The regular expression must contain three capture groups for the three variables, in order. It must match the entire read name. e.g. if field names were separated by semi-colon (';') this example regex could be specified (?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$ Note that if no READ_NAME_REGEX is specified, the read name is split on ':'. For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values."
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/sorting_collection_size_ratio",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--SORTING_COLLECTION_SIZE_RATIO"
                    },
                    "doc": "This number, plus the maximum RAM available to the JVM, determine the memory footprint used by some of the sorting collections. If you are running out of memory, try reducing this number."
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/duplicate_scoring_strategy",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--DUPLICATE_SCORING_STRATEGY"
                    },
                    "doc": "The scoring strategy for choosing the non-duplicate among candidates. Default value:SUM_OF_BASE_QUALITIES. This option can be set to 'null' to clear the default value.Possible values: {SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM}"
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/optical_duplicate_pixel_distance",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--OPTICAL_DUPLICATE_PIXEL_DISTANCE"
                    },
                    "doc": "The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is moreappropriate. For other platforms and models, users should experiment to find what works best.  Default value: 100. This option can be set to 'null' to clear the default value."
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/picard_mark_duplicates_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/,'_md.bam')\n    }\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                },
                {
                    "id": "#picard_mark_duplicates_4.1.8.1.cwl/picard_mark_duplicates_metrics",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.duplication_metrics){\n        return inputs.duplication_metrics\n    } else {\n        return inputs.input.basename.replace(/.bam/,'_md.metrics')\n    }\n}"
                    }
                }
            ],
            "label": "picard_mark_duplicates_4.1.8.1",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${ if(inputs.memory_per_job && inputs.memory_overhead) { if(inputs.memory_per_job % 1000 == 0) { return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\" } else { return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\" } } else if (inputs.memory_per_job && !inputs.memory_overhead){ if(inputs.memory_per_job % 1000 == 0) { return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\" } else { return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\" } } else if(!inputs.memory_per_job && inputs.memory_overhead){ return \"-Xmx15G\" } else { return \"-Xmx15G\" } }"
                },
                {
                    "position": 0,
                    "prefix": "-jar",
                    "valueFrom": "/gatk/gatk-package-4.1.8.1-local.jar"
                },
                {
                    "position": 0,
                    "valueFrom": "MarkDuplicates"
                },
                {
                    "position": 0,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/,'_md.bam')\n    }\n}"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 17000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.1"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:kumarn1@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Nikhil Kumar"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "picard",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.1"
                }
            ]
        },
        {
            "class": "Workflow",
            "id": "#alignment.cwl_2",
            "label": "alignment",
            "inputs": [
                {
                    "id": "#alignment.cwl_2/create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 319.15625,
                    "https://www.sevenbridges.com/y": 958.8671875
                },
                {
                    "id": "#alignment.cwl_2/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 319.15625,
                    "https://www.sevenbridges.com/y": 852.0390625
                },
                {
                    "id": "#alignment.cwl_2/read_group_description",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1495.59375
                },
                {
                    "id": "#alignment.cwl_2/read_group_identifier",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1388.765625
                },
                {
                    "id": "#alignment.cwl_2/read_group_library",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1281.9375
                },
                {
                    "id": "#alignment.cwl_2/read_group_platform_unit",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1175.109375
                },
                {
                    "id": "#alignment.cwl_2/read_group_run_date",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1068.28125
                },
                {
                    "id": "#alignment.cwl_2/read_group_sample_name",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 961.453125
                },
                {
                    "id": "#alignment.cwl_2/read_group_sequencing_center",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 854.625
                },
                {
                    "id": "#alignment.cwl_2/read_group_sequencing_platform",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 747.796875
                },
                {
                    "id": "#alignment.cwl_2/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 427.3125
                },
                {
                    "id": "#alignment.cwl_2/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 106.828125
                },
                {
                    "id": "#alignment.cwl_2/reference",
                    "type": "File",
                    "secondaryFiles": [
                        ".amb",
                        ".fai",
                        ".sa",
                        "^.dict",
                        ".ann",
                        ".bwt",
                        ".pac"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 534.140625
                },
                {
                    "id": "#alignment.cwl_2/reads",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 640.96875
                },
                {
                    "id": "#alignment.cwl_2/output",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1709.25
                },
                {
                    "id": "#alignment.cwl_2/P",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1602.421875
                },
                {
                    "id": "#alignment.cwl_2/M",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1816.078125
                },
                {
                    "id": "#alignment.cwl_2/T",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 320.484375
                },
                {
                    "id": "#alignment.cwl_2/Y",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 0
                },
                {
                    "id": "#alignment.cwl_2/K",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1922.90625
                },
                {
                    "id": "#alignment.cwl_2/bwa_number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2029.734375
                },
                {
                    "id": "#alignment.cwl_2/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 213.65625
                }
            ],
            "outputs": [
                {
                    "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_bam",
                    "outputSource": [
                        "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/picard_add_or_replace_read_groups_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 1389.239501953125,
                    "https://www.sevenbridges.com/y": 1014.8671875
                }
            ],
            "steps": [
                {
                    "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1",
                    "in": [
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/input",
                            "source": "#alignment.cwl_2/bwa_mem_0_7_17/bwa_mem_output_sam"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/output_file_name",
                            "source": "#alignment.cwl_2/output_file_name"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/sort_order",
                            "source": "#alignment.cwl_2/sort_order"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/read_group_identifier",
                            "source": "#alignment.cwl_2/read_group_identifier"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/read_group_sequencing_center",
                            "source": "#alignment.cwl_2/read_group_sequencing_center"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/read_group_library",
                            "source": "#alignment.cwl_2/read_group_library"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/read_group_platform_unit",
                            "source": "#alignment.cwl_2/read_group_platform_unit"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/read_group_sample_name",
                            "source": "#alignment.cwl_2/read_group_sample_name"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/read_group_sequencing_platform",
                            "source": "#alignment.cwl_2/read_group_sequencing_platform"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/read_group_description",
                            "source": "#alignment.cwl_2/read_group_description"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/read_group_run_date",
                            "source": "#alignment.cwl_2/read_group_run_date"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/validation_stringency",
                            "source": "#alignment.cwl_2/validation_stringency"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/create_bam_index",
                            "source": "#alignment.cwl_2/create_bam_index"
                        },
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/temporary_directory",
                            "source": "#alignment.cwl_2/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#alignment.cwl_2/picard_add_or_replace_read_groups_4_1_8_1/picard_add_or_replace_read_groups_bam"
                        }
                    ],
                    "run": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2",
                    "label": "picard_add_or_replace_read_groups_4.1.8.1",
                    "https://www.sevenbridges.com/x": 737.3328857421875,
                    "https://www.sevenbridges.com/y": 923.8671875
                },
                {
                    "id": "#alignment.cwl_2/bwa_mem_0_7_17",
                    "in": [
                        {
                            "id": "#alignment.cwl_2/bwa_mem_0_7_17/number_of_threads",
                            "source": "#alignment.cwl_2/bwa_number_of_threads"
                        },
                        {
                            "id": "#alignment.cwl_2/bwa_mem_0_7_17/reads",
                            "source": [
                                "#alignment.cwl_2/reads"
                            ]
                        },
                        {
                            "id": "#alignment.cwl_2/bwa_mem_0_7_17/reference",
                            "source": "#alignment.cwl_2/reference"
                        },
                        {
                            "id": "#alignment.cwl_2/bwa_mem_0_7_17/M",
                            "source": "#alignment.cwl_2/M"
                        },
                        {
                            "id": "#alignment.cwl_2/bwa_mem_0_7_17/P",
                            "source": "#alignment.cwl_2/P"
                        },
                        {
                            "id": "#alignment.cwl_2/bwa_mem_0_7_17/T",
                            "source": "#alignment.cwl_2/T"
                        },
                        {
                            "id": "#alignment.cwl_2/bwa_mem_0_7_17/K",
                            "source": "#alignment.cwl_2/K"
                        },
                        {
                            "id": "#alignment.cwl_2/bwa_mem_0_7_17/output",
                            "source": "#alignment.cwl_2/output"
                        },
                        {
                            "id": "#alignment.cwl_2/bwa_mem_0_7_17/Y",
                            "source": "#alignment.cwl_2/Y"
                        }
                    ],
                    "out": [
                        {
                            "id": "#alignment.cwl_2/bwa_mem_0_7_17/bwa_mem_output_sam"
                        }
                    ],
                    "run": "#bwa_mem_0.7.17.cwl_2",
                    "label": "bwa_mem_0.7.17",
                    "https://www.sevenbridges.com/x": 319.15625,
                    "https://www.sevenbridges.com/y": 1121.6953125
                }
            ],
            "requirements": [],
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:sumans@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Shalabh Suman"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:johnsoni@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Ian Jonhnson"
                }
            ],
            "https://schema.org/citation": "",
            "https://schema.org/codeRepository": "https://github.com/msk-access/cwl_subworkflows/alignment",
            "https://schema.org/contributor": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:johnsoni@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Ian Jonhnson"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9042-6213",
                    "https://schema.org/name": "Ronak Shah"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:sumans@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Shalabh Suman"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:murphyc4@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Charlie Murphy"
                }
            ],
            "https://schema.org/dateCreated": "2019-10-01",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0"
        },
        {
            "class": "CommandLineTool",
            "id": "#abra2_2.22.cwl_2",
            "baseCommand": [
                "java"
            ],
            "inputs": [
                {
                    "id": "#abra2_2.22.cwl_2/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#abra2_2.22.cwl_2/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#abra2_2.22.cwl_2/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#abra2_2.22.cwl_2/input_bam",
                    "type": [
                        "File",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--in"
                    },
                    "doc": "Required list of input sam or bam file (s) separated by comma",
                    "secondaryFiles": [
                        "^.bai"
                    ]
                },
                {
                    "id": "#abra2_2.22.cwl_2/working_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Set the temp directory (overrides java.io.tmpdir)"
                },
                {
                    "id": "#abra2_2.22.cwl_2/reference_fasta",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ref"
                    },
                    "doc": "Genome reference location",
                    "secondaryFiles": [
                        ".fai"
                    ]
                },
                {
                    "id": "#abra2_2.22.cwl_2/targets",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--targets"
                    }
                },
                {
                    "id": "#abra2_2.22.cwl_2/kmer_size",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--kmer"
                    },
                    "doc": "Optional assembly kmer size(delimit with commas if multiple sizes specified)"
                },
                {
                    "id": "#abra2_2.22.cwl_2/maximum_average_depth",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--mad"
                    },
                    "doc": "Regions with average depth exceeding this value will be downsampled (default: 1000)"
                },
                {
                    "id": "#abra2_2.22.cwl_2/soft_clip_contig",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--sc"
                    },
                    "doc": "Soft clip contig args [max_contigs,min_base_qual,frac_high_qual_bases,min_soft_clip_len] (default:16,13,80,15)"
                },
                {
                    "id": "#abra2_2.22.cwl_2/maximum_mixmatch_rate",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--mmr"
                    },
                    "doc": "Max allowed mismatch rate when mapping reads back to contigs (default: 0.05)"
                },
                {
                    "id": "#abra2_2.22.cwl_2/scoring_gap_alignments",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--sga"
                    },
                    "doc": "Scoring used for contig alignments(match, mismatch_penalty,gap_open_penalty,gap_extend_penalty) (default:8,32,48,1)"
                },
                {
                    "id": "#abra2_2.22.cwl_2/contig_anchor",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ca"
                    },
                    "doc": "Contig anchor [M_bases_at_contig_edge,max_mismatches_near_edge] (default:10,2)"
                },
                {
                    "id": "#abra2_2.22.cwl_2/window_size",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ws"
                    },
                    "doc": "Processing window size and overlap\n(size,overlap) (default: 400,200)"
                },
                {
                    "id": "#abra2_2.22.cwl_2/consensus_sequence",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--cons"
                    },
                    "doc": "Use positional consensus sequence when aligning high quality soft clipping"
                },
                {
                    "id": "#abra2_2.22.cwl_2/output_bams",
                    "type": [
                        "string",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--out"
                    },
                    "doc": "Required list of output sam or bam file (s) separated by comma"
                },
                {
                    "id": "#abra2_2.22.cwl_2/ignore_bad_assembly",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--ignore-bad-assembly"
                    },
                    "doc": "Use this option to avoid parsing errors for corrupted assemblies"
                },
                {
                    "id": "#abra2_2.22.cwl_2/bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--index"
                    },
                    "doc": "Enable BAM index generation when outputting sorted alignments (may require additonal memory)"
                },
                {
                    "id": "#abra2_2.22.cwl_2/input_vcf",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--in-vcf"
                    },
                    "doc": "VCF containing known (or suspected) variant sites.  Very large files should be avoided."
                },
                {
                    "id": "#abra2_2.22.cwl_2/no_edge_complex_indel",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--no-edge-ci"
                    },
                    "doc": "Prevent output of complex indels at read start or read end"
                },
                {
                    "id": "#abra2_2.22.cwl_2/no_sort",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--nosort"
                    },
                    "doc": "Do not attempt to sort final output"
                }
            ],
            "outputs": [
                {
                    "id": "#abra2_2.22.cwl_2/abra_realigned_bam",
                    "type": [
                        "null",
                        "File",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "${\n    return inputs.output_bams\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "label": "abra2_2.22",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if (inputs.memory_per_job && inputs.memory_overhead) {\n\n    if (inputs.memory_per_job % 1000 == 0) {\n\n      return \"-Xmx\" + (inputs.memory_per_job / 1000).toString() + \"G\"\n    }\n    else {\n\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job / 1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead) {\n\n    if (inputs.memory_per_job % 1000 == 0) {\n\n      return \"-Xmx\" + (inputs.memory_per_job / 1000).toString() + \"G\"\n    }\n    else {\n\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job / 1000)).toString() + \"G\"\n    }\n  }\n  else if (!inputs.memory_per_job && inputs.memory_overhead) {\n\n    return \"-Xmx20G\"\n  }\n  else {\n\n    return \"-Xmx20G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "prefix": "-jar",
                    "valueFrom": "/usr/local/bin/abra2.jar"
                },
                {
                    "position": 0,
                    "prefix": "--threads",
                    "valueFrom": "${\n    if(inputs.number_of_threads)\n        return inputs.number_of_threads\n    return runtime.cores\n}"
                },
                {
                    "position": 0,
                    "prefix": "--tmpdir",
                    "valueFrom": "${\n    if(inputs.working_directory)\n        return inputs.working_directory;\n      return runtime.tmpdir\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 60000,
                    "coresMin": 16
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/abra2:2.22"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:kumarn1@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Nikhil Kumar"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "abra2",
                    "http://usefulinc.com/ns/doap#revision": 2.22
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#bedtools_genomecov_v2.28.0_cv2.cwl_2",
            "baseCommand": [
                "bedtools",
                "genomecov"
            ],
            "inputs": [
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl_2/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 5,
                        "prefix": "-ibam",
                        "shellQuote": false
                    },
                    "doc": "The input file can be in BAM format (Note: BAM  must be sorted by position)",
                    "secondaryFiles": [
                        "^.bai"
                    ]
                },
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl_2/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl_2/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl_2/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl_2/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl_2/option_bedgraph",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-bg",
                        "separate": false
                    },
                    "doc": "option flag parameter to choose output file format. -bg refers to bedgraph format"
                }
            ],
            "outputs": [
                {
                    "id": "#bedtools_genomecov_v2.28.0_cv2.cwl_2/bedtools_genomecove_bedgraph",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n     if (inputs.output_file_name)\n      return inputs.output_file_name;\n    return inputs.input.basename.replace('.bam','.bedgraph');\n  }"
                    }
                }
            ],
            "label": "bedtools_genomecov",
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 20000,
                    "coresMin": 1
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/bedtools:v2.28.0_cv2"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "${\n    if (inputs.output_file_name)\n      return inputs.output_file_name;\n    return inputs.input.basename.replace('.bam','.bedgraph');\n  }",
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:sumans@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Shalabh Suman"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:sumans@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Shalabh Suman"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "bedtools",
                    "http://usefulinc.com/ns/doap#revision": "v2.28.0_cv2"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#bedtools_merge_v2.28.0_cv2.cwl_2",
            "baseCommand": [
                "bedtools",
                "merge"
            ],
            "inputs": [
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl_2/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 5,
                        "prefix": "-i",
                        "shellQuote": false
                    },
                    "doc": "BEDgraph format file generated from Bedtools Genomecov module"
                },
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl_2/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl_2/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl_2/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl_2/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "default": 0,
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl_2/distance_between_features",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-d",
                        "shellQuote": false
                    },
                    "doc": "Maximum distance between features allowed for features to be merged."
                }
            ],
            "outputs": [
                {
                    "id": "#bedtools_merge_v2.28.0_cv2.cwl_2/bedtools_merge_bed",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if (inputs.output_file_name)\n      return inputs.output_file_name;\n    return inputs.input.basename.replace('.bedgraph', '.bed');\n  }"
                    }
                }
            ],
            "label": "bedtools_merge",
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 20000,
                    "coresMin": 1
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/bedtools:v2.28.0_cv2"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "${\n    if (inputs.output_file_name)\n      return inputs.output_file_name;\n    return inputs.input.basename.replace('.bedgraph', '.bed');\n  }",
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:sumans@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Shalabh Suman"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:sumans@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Shalabh Suman"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "bedtools",
                    "http://usefulinc.com/ns/doap#revision": "v2.28.0_cv2"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "baseCommand": [
                "bwa",
                "mem"
            ],
            "inputs": [
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/reads",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 3
                    }
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/reference",
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    },
                    "secondaryFiles": [
                        ".amb",
                        ".ann",
                        ".bwt",
                        ".pac",
                        ".sa",
                        ".fai"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/A",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-A"
                    },
                    "doc": "score for a sequence match, which scales options -TdBOELU unless overridden [1]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/B",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-B"
                    },
                    "doc": "penalty for a mismatch [4]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/C",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-C"
                    },
                    "doc": "append FASTA/FASTQ comment to SAM output"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/E",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-E",
                        "itemSeparator": ","
                    },
                    "doc": "gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/L",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-L",
                        "itemSeparator": ","
                    },
                    "doc": "penalty for 5'- and 3'-end clipping [5,5]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/M",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-M"
                    }
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/O",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-O",
                        "itemSeparator": ","
                    },
                    "doc": "gap open penalties for deletions and insertions [6,6]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/P",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-P"
                    },
                    "doc": "skip pairing; mate rescue performed unless -S also in use"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/S",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-S"
                    },
                    "doc": "skip mate rescue"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/T",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-T"
                    },
                    "doc": "minimum score to output [30]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/U",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-U"
                    },
                    "doc": "penalty for an unpaired read pair [17]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/a",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-a"
                    },
                    "doc": "output all alignments for SE or unpaired PE"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/c",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-c"
                    },
                    "doc": "skip seeds with more than INT occurrences [500]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/d",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-d"
                    },
                    "doc": "off-diagonal X-dropoff [100]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/k",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-k"
                    },
                    "doc": "minimum seed length [19]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/K",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-K"
                    },
                    "doc": "process INT input bases in each batch regardless of nThreads (for reproducibility) []"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/output",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/p",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-p"
                    },
                    "doc": "smart pairing (ignoring in2.fq)"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/r",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-r"
                    },
                    "doc": "look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/v",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-v"
                    },
                    "doc": "verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/w",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-w"
                    },
                    "doc": "band width for banded alignment [100]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/y",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-y"
                    },
                    "doc": "seed occurrence for the 3rd round seeding [20]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/D",
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-D"
                    },
                    "doc": "drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/W",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-W"
                    },
                    "doc": "discard a chain if seeded bases shorter than INT [0]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/m",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-m"
                    },
                    "doc": "perform at most INT rounds of mate rescues for each read [50]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/e",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-e"
                    }
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/x",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-x"
                    },
                    "doc": "read type. Setting -x changes multiple parameters unless overridden [null] pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref) ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref) intractg: -B9 -O16 -L5  (intra-species contigs to ref)"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/H",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-H"
                    },
                    "doc": "Use hard clipping \u2019H\u2019 in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/j",
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-j"
                    },
                    "doc": "treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/he",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "int"
                        }
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-h",
                        "itemSeparator": ","
                    },
                    "doc": "if there are <INT hits with score >80% of the max score, output all in XA [5,200]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/V",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-V"
                    },
                    "doc": "output the reference FASTA header in the XR tag"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/Y",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-Y"
                    },
                    "doc": "use soft clipping for supplementary alignments"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/I",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-M"
                    }
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/R",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "STR read group header line such as '@RG\\tID -foo\\tSM -bar' [null]"
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/sample_id",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/lane_id",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/platform",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/platform_unit",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/center_name",
                    "type": [
                        "null",
                        "string"
                    ]
                },
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/library_id",
                    "type": [
                        "null",
                        "string"
                    ]
                }
            ],
            "outputs": [
                {
                    "id": "#bwa_mem_0.7.17.cwl_2/bwa_mem_output_sam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n  if (inputs.output)\n    return inputs.output;\n  return inputs.reads[0].basename.replace(/(fastq.gz)|(fq.gz)/, 'sam');\n}"
                    }
                }
            ],
            "doc": "bwa mem [-aCHMpP] [-t nThreads] [-k minSeedLen] [-w bandWidth] [-d zDropoff] [-r seedSplitRatio] [-c maxOcc] [-A matchScore] [-B mmPenalty] [-O gapOpenPen] [-E gapExtPen] [-L clipPen] [-U unpairPen] [-R RGline] [-v verboseLevel] db.prefix reads.fq [mates.fq]\nAlign 70bp-1Mbp query sequences with the BWA-MEM algorithm. Briefly, the algorithm works by seeding alignments with maximal exact matches (MEMs) and then extending seeds with the affine-gap Smith-Waterman algorithm (SW).\n\nIf mates.fq file is absent and option -p is not set, this command regards input reads are single-end. If mates.fq is present, this command assumes the i-th read in reads.fq and the i-th read in mates.fq constitute a read pair. If -p is used, the command assumes the 2i-th and the (2i+1)-th read in reads.fq constitute a read pair (such input file is said to be interleaved). In this case, mates.fq is ignored. In the paired-end mode, the mem command will infer the read orientation and the insert size distribution from a batch of reads.\n\nThe BWA-MEM algorithm performs local alignment. It may produce multiple primary alignments for different part of a query sequence. This is a crucial feature for long sequences. However, some tools such as Picard\u2019s markDuplicates does not work with split alignments. One may consider to use option -M to flag shorter split hits as secondary.",
            "label": "bwa_mem_0.7.17",
            "arguments": [
                {
                    "position": 0,
                    "prefix": "-t",
                    "valueFrom": "$(runtime.cores)"
                },
                {
                    "position": 0,
                    "prefix": "-R",
                    "valueFrom": "${\n    if (inputs.sample_id) {\n        var rg_id = \"@RG\\\\tID:\" + inputs.sample_id + \"\\\\tSM:\" + inputs.sample_id;\n        if (inputs.library_id) {\n            rg_id += \"\\\\tLB:\" + inputs.library_id;\n        } if (inputs.platform) {\n            rg_id += \"\\\\tPL:\" + inputs.platform;\n        } if (inputs.platform_unit) {\n            rg_id += \"\\\\tPU:\" + inputs.platform_unit;\n        } if (inputs.center_name) {\n            rg_id += \"\\\\tCN:\" + inputs.center_name;\n        }\n        return rg_id\n    } else {\n        return inputs.R\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ResourceRequirement",
                    "ramMin": 34000,
                    "coresMin": 16
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/bwa:0.7.17"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "stdout": "${\n  if (inputs.output)\n    return inputs.output;\n  return inputs.reads[0].basename.replace(/(fastq.gz)|(fq.gz)/, 'sam');\n}",
            "id": "#bwa_mem_0.7.17.cwl_2",
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:johnsoni@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ian Johnson"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "bwa",
                    "http://usefulinc.com/ns/doap#revision": "0.7.17"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2",
            "baseCommand": [
                "java"
            ],
            "inputs": [
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-I"
                    },
                    "doc": "Input file ( sam).  Required."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output file name (bam or sam). Not Required"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-SO"
                    },
                    "doc": "Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.Default value: null. Possible values: {unsorted, queryname, coordinate}"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/read_group_identifier",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGID"
                    },
                    "doc": "Read Group ID  Default value: 1. This option can be set to 'null' to clear the default value  Required"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/read_group_sequencing_center",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGCN"
                    },
                    "doc": "Read Group sequencing center name  Default value: null. Required"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/read_group_library",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGLB"
                    },
                    "doc": "Read Group Library.  Required"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/read_group_platform_unit",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGPU"
                    },
                    "doc": "Read Group platform unit (eg. run barcode)  Required."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/read_group_sample_name",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGSM"
                    },
                    "doc": "Read Group sample name.  Required"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/read_group_sequencing_platform",
                    "type": "string",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGPL"
                    },
                    "doc": "Read Group platform (e.g. illumina, solid)  Required."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/read_group_description",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGDS"
                    },
                    "doc": "Read Group description  Default value: null."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/read_group_run_date",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--RGDT"
                    },
                    "doc": "Read Group run date  Default value: null."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT,LENIENT, SILENT}"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/bam_compression_level",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--COMPRESSION_LEVEL"
                    },
                    "doc": "Compression level for all compressed files created (e.g. BAM and GELI). Default value:5. This option can be set to 'null' to clear the default value."
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "default": true,
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value:false. This option can be set to 'null' to clear the default value. Possible values:{true, false}"
                },
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#picard_add_or_replace_read_groups_4.1.8.1.cwl_2/picard_add_or_replace_read_groups_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n    return inputs.input.basename.replace(/.sam$/, '_srt.bam');\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "label": "picard_add_or_replace_read_groups_4.1.8.1",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx15G\"\n  }\n  else {\n      return \"-Xmx15G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "prefix": "-Djava.io.tmpdir=",
                    "separate": false,
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "shellQuote": false,
                    "valueFrom": "-XX:-UseGCOverheadLimit"
                },
                {
                    "position": 0,
                    "prefix": "-jar",
                    "valueFrom": "/gatk/gatk-package-4.1.8.1-local.jar"
                },
                {
                    "position": 0,
                    "valueFrom": "AddOrReplaceReadGroups"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name)\n        return inputs.output_file_name;\n      return inputs.input.basename.replace(/.sam$/, '_srt.bam');\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 17000,
                    "coresMin": 2
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.1"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:kumarn1@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Nikhil Kumar"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "picard",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.1"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#picard_fix_mate_information_4.1.8.1.cwl_2",
            "baseCommand": [
                "java"
            ],
            "inputs": [
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/memory_per_job",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory per job in megabytes"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/memory_overhead",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory overhead per job in megabytes"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ]
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/input",
                    "type": "File",
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-I"
                    },
                    "doc": "The input file to fix.  This option may be specified 0 or more times"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output file name (bam or sam). Not Required"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-SO"
                    },
                    "doc": "Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT.Default value: null. Possible values: {unsorted, queryname, coordinate}"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--VALIDATION_STRINGENCY"
                    },
                    "doc": "Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT,LENIENT, SILENT}"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/bam_compression_level",
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--COMPRESSION_LEVEL"
                    },
                    "doc": "Compression level for all compressed files created (e.g. BAM and GELI). Default value:5. This option can be set to 'null' to clear the default value."
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/use_jdk_deflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_DEFLATER"
                    },
                    "doc": "Use the JDK Deflater instead of the Intel Deflater for writing compressed output"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/use_jdk_inflater",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--USE_JDK_INFLATER"
                    },
                    "doc": "Use the JDK Inflater instead of the Intel Inflater for reading compressed input"
                },
                {
                    "default": true,
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--CREATE_INDEX"
                    },
                    "doc": "Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value:false. This option can be set to 'null' to clear the default value. Possible values:{true, false}"
                },
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Default value: null. This option may be specified 0 or more times."
                }
            ],
            "outputs": [
                {
                    "id": "#picard_fix_mate_information_4.1.8.1.cwl_2/picard_fix_mate_information_bam",
                    "type": "File",
                    "outputBinding": {
                        "glob": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/,'_fm.bam')\n    }\n}"
                    },
                    "secondaryFiles": [
                        "^.bai"
                    ]
                }
            ],
            "label": "picard_fix_mate_information_4.1.8.1",
            "arguments": [
                {
                    "position": 0,
                    "valueFrom": "${\n  if(inputs.memory_per_job && inputs.memory_overhead) {\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if (inputs.memory_per_job && !inputs.memory_overhead){\n    if(inputs.memory_per_job % 1000 == 0) {\n      return \"-Xmx\" + (inputs.memory_per_job/1000).toString() + \"G\"\n    }\n    else {\n      return \"-Xmx\" + Math.floor((inputs.memory_per_job/1000)).toString() + \"G\"\n    }\n  }\n  else if(!inputs.memory_per_job && inputs.memory_overhead){\n    return \"-Xmx20G\"\n  }\n  else {\n      return \"-Xmx20G\"\n  }\n}"
                },
                {
                    "position": 0,
                    "shellQuote": false,
                    "valueFrom": "-XX:-UseGCOverheadLimit"
                },
                {
                    "position": 0,
                    "prefix": "-jar",
                    "valueFrom": "/gatk/gatk-package-4.1.8.1-local.jar"
                },
                {
                    "position": 0,
                    "valueFrom": "FixMateInformation"
                },
                {
                    "position": 0,
                    "prefix": "--TMP_DIR",
                    "valueFrom": "${\n    if(inputs.temporary_directory)\n        return inputs.temporary_directory;\n      return runtime.tmpdir\n}"
                },
                {
                    "position": 0,
                    "prefix": "-O",
                    "valueFrom": "${\n    if(inputs.output_file_name){\n        return inputs.output_file_name\n    } else {\n        return inputs.input.basename.replace(/.bam/,'_fm.bam')\n    }\n}"
                }
            ],
            "requirements": [
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "ResourceRequirement",
                    "ramMin": 30000,
                    "coresMin": 12
                },
                {
                    "class": "DockerRequirement",
                    "dockerPull": "ghcr.io/msk-access/gatk:4.1.8.1"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "http://purl.org/dc/terms/contributor": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:kumarn1@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Nikhil Kumar"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://purl.org/dc/terms/creator": [
                {
                    "class": "http://xmlns.com/foaf/0.1/Organization",
                    "http://xmlns.com/foaf/0.1/member": [
                        {
                            "class": "http://xmlns.com/foaf/0.1/Person",
                            "http://xmlns.com/foaf/0.1/mbox": "mailto:shahr2@mskcc.org",
                            "http://xmlns.com/foaf/0.1/name": "Ronak Shah"
                        }
                    ],
                    "http://xmlns.com/foaf/0.1/name": "Memorial Sloan Kettering Cancer Center"
                }
            ],
            "http://usefulinc.com/ns/doap#release": [
                {
                    "class": "http://usefulinc.com/ns/doap#Version",
                    "http://usefulinc.com/ns/doap#name": "picard",
                    "http://usefulinc.com/ns/doap#revision": "4.1.8.1"
                }
            ]
        },
        {
            "class": "Workflow",
            "id": "#indel_realignment.cwl_2",
            "label": "indel_realignment",
            "inputs": [
                {
                    "id": "#indel_realignment.cwl_2/window_size",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 0
                },
                {
                    "id": "#indel_realignment.cwl_2/soft_clip_contig",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 426.796875
                },
                {
                    "id": "#indel_realignment.cwl_2/scoring_gap_alignments",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 533.53125
                },
                {
                    "id": "#indel_realignment.cwl_2/reference_fasta",
                    "type": "File",
                    "secondaryFiles": [
                        ".fai"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 640.21875
                },
                {
                    "id": "#indel_realignment.cwl_2/no_sort",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1066.875
                },
                {
                    "id": "#indel_realignment.cwl_2/maximum_mixmatch_rate",
                    "type": [
                        "null",
                        "float"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1280.25
                },
                {
                    "id": "#indel_realignment.cwl_2/maximum_average_depth",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1386.9375
                },
                {
                    "id": "#indel_realignment.cwl_2/input_bam",
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1493.625
                },
                {
                    "id": "#indel_realignment.cwl_2/ignore_bad_assembly",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1600.3125
                },
                {
                    "id": "#indel_realignment.cwl_2/contig_anchor",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1813.6875
                },
                {
                    "id": "#indel_realignment.cwl_2/consensus_sequence",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1920.375
                },
                {
                    "id": "#indel_realignment.cwl_2/bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2027.015625
                },
                {
                    "id": "#indel_realignment.cwl_2/number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 960.234375
                },
                {
                    "id": "#indel_realignment.cwl_2/option_bedgraph",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 853.546875
                },
                {
                    "id": "#indel_realignment.cwl_2/no_edge_complex_indel",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1173.5625
                },
                {
                    "id": "#indel_realignment.cwl_2/distance_between_features",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1707
                },
                {
                    "id": "#indel_realignment.cwl_2/output_bams",
                    "type": [
                        "string",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 746.859375
                },
                {
                    "id": "#indel_realignment.cwl_2/validation_stringency",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 106.6875
                },
                {
                    "id": "#indel_realignment.cwl_2/sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 320.109375
                },
                {
                    "id": "#indel_realignment.cwl_2/output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 992.927978515625,
                    "https://www.sevenbridges.com/y": 794.8671875
                },
                {
                    "id": "#indel_realignment.cwl_2/create_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 992.927978515625,
                    "https://www.sevenbridges.com/y": 901.5078125
                },
                {
                    "id": "#indel_realignment.cwl_2/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 213.421875
                }
            ],
            "outputs": [
                {
                    "id": "#indel_realignment.cwl_2/indel_realignment_bam",
                    "outputSource": [
                        "#indel_realignment.cwl_2/picard_fix_mate_information_4_1_8_1/picard_fix_mate_information_bam"
                    ],
                    "type": "File",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 1981.323974609375,
                    "https://www.sevenbridges.com/y": 1013.4609375
                }
            ],
            "steps": [
                {
                    "id": "#indel_realignment.cwl_2/abra2_2_22",
                    "in": [
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/number_of_threads",
                            "source": "#indel_realignment.cwl_2/number_of_threads"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/input_bam",
                            "source": [
                                "#indel_realignment.cwl_2/input_bam"
                            ]
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/working_directory",
                            "source": "#indel_realignment.cwl_2/temporary_directory"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/reference_fasta",
                            "source": "#indel_realignment.cwl_2/reference_fasta"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/targets",
                            "source": "#indel_realignment.cwl_2/bedtools_merge/bedtools_merge_bed"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/maximum_average_depth",
                            "source": "#indel_realignment.cwl_2/maximum_average_depth"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/soft_clip_contig",
                            "source": "#indel_realignment.cwl_2/soft_clip_contig"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/maximum_mixmatch_rate",
                            "source": "#indel_realignment.cwl_2/maximum_mixmatch_rate"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/scoring_gap_alignments",
                            "source": "#indel_realignment.cwl_2/scoring_gap_alignments"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/contig_anchor",
                            "source": "#indel_realignment.cwl_2/contig_anchor"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/window_size",
                            "source": "#indel_realignment.cwl_2/window_size"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/consensus_sequence",
                            "source": "#indel_realignment.cwl_2/consensus_sequence"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/output_bams",
                            "source": [
                                "#indel_realignment.cwl_2/output_bams"
                            ]
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/ignore_bad_assembly",
                            "source": "#indel_realignment.cwl_2/ignore_bad_assembly"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/bam_index",
                            "source": "#indel_realignment.cwl_2/bam_index"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/no_edge_complex_indel",
                            "source": "#indel_realignment.cwl_2/no_edge_complex_indel"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/no_sort",
                            "source": "#indel_realignment.cwl_2/no_sort"
                        }
                    ],
                    "out": [
                        {
                            "id": "#indel_realignment.cwl_2/abra2_2_22/abra_realigned_bam"
                        }
                    ],
                    "run": "#abra2_2.22.cwl_2",
                    "label": "abra2_2.22",
                    "https://www.sevenbridges.com/x": 992.927978515625,
                    "https://www.sevenbridges.com/y": 1120.1484375
                },
                {
                    "id": "#indel_realignment.cwl_2/bedtools_genomecov",
                    "in": [
                        {
                            "id": "#indel_realignment.cwl_2/bedtools_genomecov/input",
                            "source": "#indel_realignment.cwl_2/input_bam"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/bedtools_genomecov/option_bedgraph",
                            "source": "#indel_realignment.cwl_2/option_bedgraph"
                        }
                    ],
                    "out": [
                        {
                            "id": "#indel_realignment.cwl_2/bedtools_genomecov/bedtools_genomecove_bedgraph"
                        }
                    ],
                    "run": "#bedtools_genomecov_v2.28.0_cv2.cwl_2",
                    "label": "bedtools_genomecov",
                    "https://www.sevenbridges.com/x": 269.59375,
                    "https://www.sevenbridges.com/y": 1006.4609375
                },
                {
                    "id": "#indel_realignment.cwl_2/bedtools_merge",
                    "in": [
                        {
                            "id": "#indel_realignment.cwl_2/bedtools_merge/input",
                            "source": "#indel_realignment.cwl_2/bedtools_genomecov/bedtools_genomecove_bedgraph"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/bedtools_merge/distance_between_features",
                            "source": "#indel_realignment.cwl_2/distance_between_features"
                        }
                    ],
                    "out": [
                        {
                            "id": "#indel_realignment.cwl_2/bedtools_merge/bedtools_merge_bed"
                        }
                    ],
                    "run": "#bedtools_merge_v2.28.0_cv2.cwl_2",
                    "label": "bedtools_merge",
                    "https://www.sevenbridges.com/x": 635.5108642578125,
                    "https://www.sevenbridges.com/y": 1006.4609375
                },
                {
                    "id": "#indel_realignment.cwl_2/picard_fix_mate_information_4_1_8_1",
                    "in": [
                        {
                            "id": "#indel_realignment.cwl_2/picard_fix_mate_information_4_1_8_1/input",
                            "source": "#indel_realignment.cwl_2/abra2_2_22/abra_realigned_bam"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/picard_fix_mate_information_4_1_8_1/output_file_name",
                            "source": "#indel_realignment.cwl_2/output_file_name"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/picard_fix_mate_information_4_1_8_1/sort_order",
                            "source": "#indel_realignment.cwl_2/sort_order"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/picard_fix_mate_information_4_1_8_1/validation_stringency",
                            "source": "#indel_realignment.cwl_2/validation_stringency"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/picard_fix_mate_information_4_1_8_1/create_bam_index",
                            "source": "#indel_realignment.cwl_2/create_bam_index"
                        },
                        {
                            "id": "#indel_realignment.cwl_2/picard_fix_mate_information_4_1_8_1/temporary_directory",
                            "source": "#indel_realignment.cwl_2/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#indel_realignment.cwl_2/picard_fix_mate_information_4_1_8_1/picard_fix_mate_information_bam"
                        }
                    ],
                    "run": "#picard_fix_mate_information_4.1.8.1.cwl_2",
                    "label": "picard_fix_mate_information_4.1.8.1",
                    "https://www.sevenbridges.com/x": 1546.70458984375,
                    "https://www.sevenbridges.com/y": 978.328125
                }
            ],
            "requirements": [],
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:murphyc4@mskcc.org",
                    "https://schema.org/name": "Charlie Murphy"
                }
            ],
            "https://schema.org/citation": "",
            "https://schema.org/codeRepository": "https://github.com/msk-access/cwl_subworkflows/indel_realignment",
            "https://schema.org/contributor": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:murphyc4@mskcc.org",
                    "https://schema.org/name": "Charlie Murphy"
                }
            ],
            "https://schema.org/dateCreated": "2020-09-14",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0"
        },
        {
            "class": "Workflow",
            "id": "#uncollapsed_bam_generation.cwl",
            "doc": "This is the workflow is written using Common Workflow Language (CWL) version 1.0 (https://www.commonwl.org/v1.0/) and is used at Memorial Sloan Kettering Cancer Center for producing standard bam files from the NY state-approved MSK-ACCESS assay. It is meant to be run on a single sample paired-end read 1 and read 2 fastq's, from Illumina sequencing data, but may be generalizable to other sequencing platforms and sequencing panels as well, which produce paired-end data.",
            "label": "Uncollapsed BAM Generation",
            "inputs": [
                {
                    "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_umi-tag",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio FastqToBam: Tag in which to store molecular barcodes/UMIs.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2670.3125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_sort",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Fgbio FastqToBam: If true, queryname sort the BAM file, otherwise  preserve input order.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2777.125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/sequencing-center",
                    "type": "string",
                    "doc": "The sequencing center from which the data originated",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 427.25
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/sample",
                    "type": "string",
                    "doc": "The name of the sequenced sample.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 534.0625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/run-date",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Date the run was produced, to insert into the read group header  (Iso8601Date)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 640.875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/read-structures",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "doc": "Fgbio FastqToBam: Read structures, one for each of the FASTQs.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 854.5
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/read-group-id",
                    "type": "string",
                    "doc": "Read group ID to use in the file header.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 961.3125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_predicted-insert-size",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fgbio FastqToBam: Predicted median insert size, to insert into the read  group header",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2883.9375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/platform-unit",
                    "type": "string",
                    "doc": "Platform unit (e.g. \"..\")",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1388.5625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/platform-model",
                    "type": "string",
                    "doc": "Platform model to insert into the group header (ex. miseq, hiseq2500,  hiseqX)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1495.375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/platform",
                    "type": "string",
                    "doc": "Sequencing Platform.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1602.1875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fgbio FastqToBam: The output SAM or BAM file to be written.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2990.75
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/library",
                    "type": "string",
                    "doc": "The name/ID of the sequenced library.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2136.25
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/description",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Description of the read group.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4272.5
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/comment",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Comments to include in the output file\u2019s header.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4379.3125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/validation_stringency",
                    "type": "string",
                    "doc": "Validation stringency for all SAM files read by this program. Setting  stringency to SILENT can improve performance when processing a BAM file  in which variable-length data (read, qualities, tags) do not otherwise  need to be decoded. The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: STRICT or LENIENT or SILENT",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 0
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/unpaired_fastq_file",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Name of the Unpaired Fastq File",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 106.8125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_include_non_primary_alignments",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "\tIf true, include non-primary alignments in the output. Support of \n\tnon-primary alignments in SamToFastq is not comprehensive, so there \n\tmay be exceptions if this is set to true and there are paired reads \n\twith non-primary alignments.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2243.0625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_include_non_pf_reads",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Include non-PF reads from the SAM file into the output FASTQ files.  PF means 'passes filtering'. Reads whose 'not passing quality controls'  flag is set are non-PF reads. See GATK Dictionary for more info.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2349.875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/R1_output_fastq",
                    "type": "string",
                    "doc": "Name of the R1 output Fastq File",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1281.75
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/R2_output_fastq",
                    "type": "string",
                    "doc": "Name of the R2 Fastq File",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1174.9375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/reference_sequence",
                    "type": "File",
                    "doc": "Reference sequence file.  Please include \".fai\", \"^.dict\", \".amb\" , \".sa\", \".bwt\", \".pac\", \".ann\" as  secondary files if they are not present in the same location as the  \".fasta\" file",
                    "secondaryFiles": [
                        ".amb",
                        ".fai",
                        ".sa",
                        "^.dict",
                        ".ann",
                        ".bwt",
                        ".pac"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 747.6875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_unpaired2_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fastp: For PE input, if read2 passed QC but read1 not, it will be written to  unpaired2. If --unpaired2 is same as --unpaired1 (default mode), both  unpaired reads will be written to this same file.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3204.375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_unpaired1_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fastp: for PE input, if read1 passed QC but read2 not, it will be  written to unpaired1. Default is to discard it.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3311.1875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_read2_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fastp: Read2 output File Name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3418
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_read2_adapter_sequence",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fastp: The adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as  <adapter_sequence> (string)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3524.8125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_read1_output_file_name",
                    "type": "string",
                    "doc": "Fastp: Read1 output File Name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3631.625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_read1_adapter_sequence",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fastp: the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3738.4375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_minimum_read_length",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Fastp: reads shorter than length_required will be discarded, default is 15.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3845.25
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_json_output_file_name",
                    "type": "string",
                    "doc": "Fastp: the json format report file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3952.0625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_html_output_file_name",
                    "type": "string",
                    "doc": "Fastp: the html format report file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4058.875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_failed_reads_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Fastp: specify the file to store reads that cannot pass the filters.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4165.6875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/bwa_mem_Y",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "BWA MEM: use soft clipping for supplementary alignments",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4592.9375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/bwa_mem_T",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "BWA MEM: minimum score to output [30]",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4699.75
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/sort_order",
                    "type": "string",
                    "doc": "GATK: The order in which the reads should be output.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 320.4375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/bwa_mem_P",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "BWA MEM: skip pairing; mate rescue performed unless -S also in use",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4806.5625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/picard_addRG_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output BAM file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1922.625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/bwa_mem_output",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output SAM file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4913.375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/bwa_mem_M",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "BWA MEM: mark shorter split hits as secondary",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5020.1875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/bwa_mem_K",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "process INT input bases in each batch regardless of nThreads (for reproducibility)",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 5127
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/create_bam_index",
                    "type": "boolean",
                    "doc": "GATK: Generate BAM index file when possible",
                    "https://www.sevenbridges.com/x": 1456.3748779296875,
                    "https://www.sevenbridges.com/y": 2723.71875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Output BAM file name",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2563.5
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/optical_duplicate_pixel_distance",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Picard MarkDuplicates: The maximum offset between two duplicate clusters  in order to consider them optical duplicates. The default is appropriate  for unpatterned versions of the Illumina platform. For the patterned  flowcell models, 2500 is more appropriate. For other platforms and models,  users should experiment to find what works best.",
                    "https://www.sevenbridges.com/x": 2024.1610107421875,
                    "https://www.sevenbridges.com/y": 2056.65625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/duplicate_scoring_strategy",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Picard MarkDuplicates: The scoring strategy for choosing the  non-duplicate among candidates.",
                    "https://www.sevenbridges.com/x": 2024.1610107421875,
                    "https://www.sevenbridges.com/y": 3070.34375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/read_name_regex",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Picard MarkDuplicates: Regular expression that can be used to parse read  names in the incoming SAM file. Read names are parsed to extract three  variables: tile/region, x coordinate and y coordinate. These values are  used to estimate the rate of optical duplication in order to give a more  accurate estimated library size. Set this option to null to disable  optical duplicate detection, e.g. for RNA-seq or other data where  duplicate sets are extremely large and estimating library complexity is not an aim. Note that without optical duplicate counts, library size  estimation will be inaccurate. The regular expression should contain  three capture groups for the three variables, in order. It must match the  entire read name. Note that if the default regex is specified, a regex  match is not actually done, but instead the read name is split on colon  character. For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. For 7 element names (CASAVA 1.8), the 5th,  6th, and 7th elements are assumed to be tile, x and y values.",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1068.125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_mark_duplicates_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Picard MarkDuplicates: The output file to write marked records to",
                    "https://www.sevenbridges.com/x": 2024.1610107421875,
                    "https://www.sevenbridges.com/y": 2461.09375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_mark_duplicates_duplication_metrics_file_name",
                    "type": "string",
                    "doc": "Picard MarkDuplicates: File to write duplication metrics to",
                    "https://www.sevenbridges.com/x": 2024.1610107421875,
                    "https://www.sevenbridges.com/y": 2567.90625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_mark_duplicates_assume_sort_order",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Picard MarkDuplicates: If not null, assume that the input file has this  order even if the header says otherwise.",
                    "https://www.sevenbridges.com/x": 2024.1610107421875,
                    "https://www.sevenbridges.com/y": 2674.71875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_window_size",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "ABRA2: Processing window size and overlap (size,overlap)  (default: 400,200)",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 2603.6875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_soft_clip_contig",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "ABRA2: Soft clip contig args [maxcontigs,min_base_qual,frac  high_qual_bases,min_soft_clip_len] (default:16,13,80,15)",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 2710.5
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_scoring_gap_alignments",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 2817.3125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_output_bams",
                    "type": [
                        "string",
                        {
                            "type": "array",
                            "items": "string"
                        }
                    ],
                    "doc": "Required list of output sam or bam file",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 2924.125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_maximum_average_depth",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "ABRA2: Regions with average depth exceeding this value will be  downsampled (default: 1000)",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 3351.375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_bam_index",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "ABRA2: Generate BAM Index",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 3671.8125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_contig_anchor",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 3458.1875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_consensus_sequence",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "ABRA2: Contig anchor [M_bases_at_contig_edge,max_mismatches_near_edge]  (default:10,2)",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 3565
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/bedtools_merge_distance_between_features",
                    "type": [
                        "null",
                        "int"
                    ],
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 1989.25
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_maximum_mixmatch_rate",
                    "type": [
                        "null",
                        "float"
                    ],
                    "doc": "max allowed mismatch rate when mapping\nreads back to contigs (default: 0.05)",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 3244.5625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/bedtools_genomecov_option_bedgraph",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "bedtools genomecov: option flag parameter to choose output file format.  -bg refers to bedgraph format",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 2096.0625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/picard_fixmateinformation_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Picard FixMateInformation: The output BAM file to write to",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1709
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_no_sort",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "ABRA2: Do not attempt to sort final output",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 3030.9375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/abra2_no_edge_complex_indel",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "ABRA2: Prevent output of complex indels at read start or read end",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 3137.75
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/merge_sam_files_sort_order",
                    "type": "string",
                    "doc": "GATK MergeSamFiles: Sort order of output file",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2029.4375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_output_file_name",
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "GATK MergeSamFiles: SAM or BAM file to write merged result to",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 2456.6875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/bwa_number_of_threads",
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "BWA MEM: Number of threads",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 4486.125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_input",
                    "type": {
                        "type": "array",
                        "items": {
                            "items": "File",
                            "type": "array"
                        }
                    },
                    "doc": "Fgbio FastqToBam: Fastq files corresponding to each sequencing read ( e.g. R1, I1, etc.).",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 3097.5625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/picard_addRG_sort_order",
                    "type": "string",
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 1815.8125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/disable_trim_poly_g",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 1456.3748779296875,
                    "https://www.sevenbridges.com/y": 2510.09375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/disable_quality_filtering",
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "https://www.sevenbridges.com/x": 1456.3748779296875,
                    "https://www.sevenbridges.com/y": 2616.90625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/temporary_directory",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 0,
                    "https://www.sevenbridges.com/y": 213.625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fgbio_async_io",
                    "type": [
                        "null",
                        "string"
                    ],
                    "https://www.sevenbridges.com/x": 383.6767272949219,
                    "https://www.sevenbridges.com/y": 2878.35009765625
                }
            ],
            "outputs": [
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_unpaired_fastq",
                    "outputSource": [
                        "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_unpaired_fastq"
                    ],
                    "type": [
                        "null",
                        "File"
                    ],
                    "https://www.sevenbridges.com/x": 2024.1610107421875,
                    "https://www.sevenbridges.com/y": 2163.46875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_unpaired2_output",
                    "outputSource": [
                        "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_unpaired2_output"
                    ],
                    "type": [
                        "null",
                        "File"
                    ],
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 1562
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_unpaired1_output",
                    "outputSource": [
                        "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_unpaired1_output"
                    ],
                    "type": [
                        "null",
                        "File"
                    ],
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 1668.8125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_json_output",
                    "outputSource": [
                        "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_json_output"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 1775.625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_html_output",
                    "outputSource": [
                        "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_html_output"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 1882.4375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_metrics",
                    "outputSource": [
                        "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/picard_mark_duplicates_metrics"
                    ],
                    "type": "File",
                    "https://www.sevenbridges.com/x": 3310.43603515625,
                    "https://www.sevenbridges.com/y": 2377.09375
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/indel_realignment_bam",
                    "outputSource": [
                        "#uncollapsed_bam_generation.cwl/indel_realignment/indel_realignment_bam"
                    ],
                    "type": "File",
                    "doc": "This bam file will be used for collapsing",
                    "secondaryFiles": [
                        "^.bai"
                    ],
                    "https://www.sevenbridges.com/x": 3918.8408203125,
                    "https://www.sevenbridges.com/y": 2563.5
                }
            ],
            "steps": [
                {
                    "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0",
                    "in": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/input",
                            "source": [
                                "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_input"
                            ]
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/output_file_name",
                            "source": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/read-structures",
                            "source": [
                                "#uncollapsed_bam_generation.cwl/read-structures"
                            ]
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/sort",
                            "source": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_sort"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/umi-tag",
                            "source": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_umi-tag"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/read-group-id",
                            "source": "#uncollapsed_bam_generation.cwl/read-group-id"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/sample",
                            "source": "#uncollapsed_bam_generation.cwl/sample"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/library",
                            "source": "#uncollapsed_bam_generation.cwl/library"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/platform",
                            "source": "#uncollapsed_bam_generation.cwl/platform"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/platform-unit",
                            "source": "#uncollapsed_bam_generation.cwl/platform-unit"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/platform-model",
                            "source": "#uncollapsed_bam_generation.cwl/platform-model"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/sequencing-center",
                            "source": "#uncollapsed_bam_generation.cwl/sequencing-center"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/predicted-insert-size",
                            "source": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_predicted-insert-size"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/description",
                            "source": "#uncollapsed_bam_generation.cwl/description"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/comment",
                            "source": "#uncollapsed_bam_generation.cwl/comment"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/run-date",
                            "source": "#uncollapsed_bam_generation.cwl/run-date"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/temporary_directory",
                            "source": "#uncollapsed_bam_generation.cwl/temporary_directory"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/async_io",
                            "source": "#uncollapsed_bam_generation.cwl/fgbio_async_io"
                        }
                    ],
                    "out": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/fgbio_fastq_to_bam_ubam"
                        }
                    ],
                    "run": "#fgbio_fastq_to_bam_1.2.0.cwl",
                    "label": "fgbio_fastq_to_bam_1.2.0",
                    "scatter": [
                        "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/input"
                    ],
                    "scatterMethod": "dotproduct",
                    "https://www.sevenbridges.com/x": 477.953125,
                    "https://www.sevenbridges.com/y": 2451.5
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0",
                    "in": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/fastq",
                            "source": "#uncollapsed_bam_generation.cwl/R1_output_fastq"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/input",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_4_1_8_0/gatk_merge_sam_files_bam"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/include_non_pf_reads",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_include_non_pf_reads"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/include_non_primary_alignments",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_include_non_primary_alignments"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/reference_sequence",
                            "source": "#uncollapsed_bam_generation.cwl/reference_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/second_end_fastq",
                            "source": "#uncollapsed_bam_generation.cwl/R2_output_fastq"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/unpaired_fastq",
                            "source": "#uncollapsed_bam_generation.cwl/unpaired_fastq_file"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/validation_stringency",
                            "source": "#uncollapsed_bam_generation.cwl/validation_stringency"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/temporary_directory",
                            "source": "#uncollapsed_bam_generation.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_fastq"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_unpaired_fastq"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_second_end_fastq"
                        }
                    ],
                    "run": "#gatk_sam_to_fastq_4.1.8.0.cwl_2",
                    "label": "GATK-SamToFastq",
                    "https://www.sevenbridges.com/x": 1456.3748779296875,
                    "https://www.sevenbridges.com/y": 2347.28125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1",
                    "in": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/read1_input",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_fastq"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/read1_output_path",
                            "source": "#uncollapsed_bam_generation.cwl/fastp_read1_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/read2_input",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_sam_to_fastq_4_1_8_0/gatk_sam_to_fastq_second_end_fastq"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/read2_output_path",
                            "source": "#uncollapsed_bam_generation.cwl/fastp_read2_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/unpaired1_path",
                            "source": "#uncollapsed_bam_generation.cwl/fastp_unpaired1_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/unpaired2_path",
                            "source": "#uncollapsed_bam_generation.cwl/fastp_unpaired2_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/failed_reads_path",
                            "source": "#uncollapsed_bam_generation.cwl/fastp_failed_reads_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/read1_adapter_sequence",
                            "source": "#uncollapsed_bam_generation.cwl/fastp_read1_adapter_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/read2_adapter_sequence",
                            "source": "#uncollapsed_bam_generation.cwl/fastp_read2_adapter_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/minimum_read_length",
                            "source": "#uncollapsed_bam_generation.cwl/fastp_minimum_read_length"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/json_output_path",
                            "source": "#uncollapsed_bam_generation.cwl/fastp_json_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/html_output_path",
                            "source": "#uncollapsed_bam_generation.cwl/fastp_html_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/disable_quality_filtering",
                            "source": "#uncollapsed_bam_generation.cwl/disable_quality_filtering"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/disable_trim_poly_g",
                            "source": "#uncollapsed_bam_generation.cwl/disable_trim_poly_g"
                        }
                    ],
                    "out": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_json_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_html_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_read1_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_read2_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_unpaired1_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_unpaired2_output"
                        }
                    ],
                    "run": "#fastp_0.20.1.cwl",
                    "label": "fastp_0.20.1",
                    "https://www.sevenbridges.com/x": 2024.1610107421875,
                    "https://www.sevenbridges.com/y": 2872.53125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/alignment",
                    "in": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/create_bam_index",
                            "source": "#uncollapsed_bam_generation.cwl/create_bam_index"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/output_file_name",
                            "source": "#uncollapsed_bam_generation.cwl/picard_addRG_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/read_group_description",
                            "source": "#uncollapsed_bam_generation.cwl/description"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/read_group_identifier",
                            "source": "#uncollapsed_bam_generation.cwl/read-group-id"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/read_group_library",
                            "source": "#uncollapsed_bam_generation.cwl/library"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/read_group_platform_unit",
                            "source": "#uncollapsed_bam_generation.cwl/platform-unit"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/read_group_run_date",
                            "source": "#uncollapsed_bam_generation.cwl/run-date"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/read_group_sample_name",
                            "source": "#uncollapsed_bam_generation.cwl/sample"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/read_group_sequencing_center",
                            "source": "#uncollapsed_bam_generation.cwl/sequencing-center"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/read_group_sequencing_platform",
                            "source": "#uncollapsed_bam_generation.cwl/platform"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/sort_order",
                            "source": "#uncollapsed_bam_generation.cwl/picard_addRG_sort_order"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/validation_stringency",
                            "source": "#uncollapsed_bam_generation.cwl/validation_stringency"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/reference",
                            "source": "#uncollapsed_bam_generation.cwl/reference_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/reads",
                            "source": [
                                "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_read1_output",
                                "#uncollapsed_bam_generation.cwl/fastp_0_20_1/fastp_read2_output"
                            ]
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/output",
                            "source": "#uncollapsed_bam_generation.cwl/bwa_mem_output"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/P",
                            "source": "#uncollapsed_bam_generation.cwl/bwa_mem_P"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/M",
                            "source": "#uncollapsed_bam_generation.cwl/bwa_mem_M"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/T",
                            "source": "#uncollapsed_bam_generation.cwl/bwa_mem_T"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/Y",
                            "source": "#uncollapsed_bam_generation.cwl/bwa_mem_Y"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/K",
                            "source": "#uncollapsed_bam_generation.cwl/bwa_mem_K"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/bwa_number_of_threads",
                            "source": "#uncollapsed_bam_generation.cwl/bwa_number_of_threads"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/temporary_directory",
                            "source": "#uncollapsed_bam_generation.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/alignment/picard_add_or_replace_read_groups_bam"
                        }
                    ],
                    "run": "#alignment.cwl_2",
                    "label": "alignment",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 2349.875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_4_1_8_0",
                    "in": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_4_1_8_0/unmapped_bam",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_4_1_8_0/gatk_merge_sam_files_bam"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_4_1_8_0/reference",
                            "source": "#uncollapsed_bam_generation.cwl/reference_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_4_1_8_0/output_file_name",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_4_1_8_0/aligned_bam",
                            "source": [
                                "#uncollapsed_bam_generation.cwl/alignment/picard_add_or_replace_read_groups_bam"
                            ],
                            "valueFrom": "${ return [ self ]; }"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_4_1_8_0/validation_stringency",
                            "source": "#uncollapsed_bam_generation.cwl/validation_stringency"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_4_1_8_0/create_index",
                            "source": "#uncollapsed_bam_generation.cwl/create_bam_index"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_4_1_8_0/temporary_directory",
                            "source": "#uncollapsed_bam_generation.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_4_1_8_0/gatk_merge_bam_alignment_bam"
                        }
                    ],
                    "run": "#gatk_merge_bam_alignment_4.1.8.0.cwl_2",
                    "label": "GATK-MergeBamAlignment",
                    "https://www.sevenbridges.com/x": 2024.1610107421875,
                    "https://www.sevenbridges.com/y": 2312.28125
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1",
                    "in": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/input",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_merge_bam_alignment_4_1_8_0/gatk_merge_bam_alignment_bam"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/output_file_name",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_mark_duplicates_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/duplication_metrics",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_mark_duplicates_duplication_metrics_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/assume_sort_order",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_mark_duplicates_assume_sort_order"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/validation_stringency",
                            "source": "#uncollapsed_bam_generation.cwl/validation_stringency"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/create_bam_index",
                            "source": "#uncollapsed_bam_generation.cwl/create_bam_index"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/read_name_regex",
                            "source": "#uncollapsed_bam_generation.cwl/read_name_regex"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/duplicate_scoring_strategy",
                            "source": "#uncollapsed_bam_generation.cwl/duplicate_scoring_strategy"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/optical_duplicate_pixel_distance",
                            "source": "#uncollapsed_bam_generation.cwl/optical_duplicate_pixel_distance"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/temporary_directory",
                            "source": "#uncollapsed_bam_generation.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/picard_mark_duplicates_bam"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/picard_mark_duplicates_metrics"
                        }
                    ],
                    "run": "#picard_mark_duplicates_4.1.8.1.cwl",
                    "label": "picard_mark_duplicates_4.1.8.1",
                    "https://www.sevenbridges.com/x": 2538.884765625,
                    "https://www.sevenbridges.com/y": 1392.1875
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/indel_realignment",
                    "in": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/window_size",
                            "source": "#uncollapsed_bam_generation.cwl/abra2_window_size"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/soft_clip_contig",
                            "source": "#uncollapsed_bam_generation.cwl/abra2_soft_clip_contig"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/scoring_gap_alignments",
                            "source": "#uncollapsed_bam_generation.cwl/abra2_scoring_gap_alignments"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/reference_fasta",
                            "source": "#uncollapsed_bam_generation.cwl/reference_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/no_sort",
                            "source": "#uncollapsed_bam_generation.cwl/abra2_no_sort"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/maximum_mixmatch_rate",
                            "source": "#uncollapsed_bam_generation.cwl/abra2_maximum_mixmatch_rate"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/maximum_average_depth",
                            "source": "#uncollapsed_bam_generation.cwl/abra2_maximum_average_depth"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/input_bam",
                            "source": "#uncollapsed_bam_generation.cwl/picard_mark_duplicates_4_1_8_1/picard_mark_duplicates_bam"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/contig_anchor",
                            "source": "#uncollapsed_bam_generation.cwl/abra2_contig_anchor"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/consensus_sequence",
                            "source": "#uncollapsed_bam_generation.cwl/abra2_consensus_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/bam_index",
                            "source": "#uncollapsed_bam_generation.cwl/abra2_bam_index"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/option_bedgraph",
                            "source": "#uncollapsed_bam_generation.cwl/bedtools_genomecov_option_bedgraph"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/no_edge_complex_indel",
                            "source": "#uncollapsed_bam_generation.cwl/abra2_no_edge_complex_indel"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/distance_between_features",
                            "source": "#uncollapsed_bam_generation.cwl/bedtools_merge_distance_between_features"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/output_bams",
                            "source": [
                                "#uncollapsed_bam_generation.cwl/abra2_output_bams"
                            ]
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/validation_stringency",
                            "source": "#uncollapsed_bam_generation.cwl/validation_stringency"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/sort_order",
                            "source": "#uncollapsed_bam_generation.cwl/sort_order"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/output_file_name",
                            "source": "#uncollapsed_bam_generation.cwl/picard_fixmateinformation_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/create_bam_index",
                            "source": "#uncollapsed_bam_generation.cwl/create_bam_index"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/temporary_directory",
                            "source": "#uncollapsed_bam_generation.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/indel_realignment/indel_realignment_bam"
                        }
                    ],
                    "run": "#indel_realignment.cwl_2",
                    "label": "indel_realignment",
                    "https://www.sevenbridges.com/x": 3310.43603515625,
                    "https://www.sevenbridges.com/y": 2616.90625
                },
                {
                    "id": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_4_1_8_0",
                    "in": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_4_1_8_0/input",
                            "source": [
                                "#uncollapsed_bam_generation.cwl/fgbio_fastq_to_bam_1_2_0/fgbio_fastq_to_bam_ubam"
                            ]
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_4_1_8_0/output_file_name",
                            "source": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_output_file_name"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_4_1_8_0/reference_sequence",
                            "source": "#uncollapsed_bam_generation.cwl/reference_sequence"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_4_1_8_0/sort_order",
                            "source": "#uncollapsed_bam_generation.cwl/merge_sam_files_sort_order"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_4_1_8_0/validation_stringency",
                            "source": "#uncollapsed_bam_generation.cwl/validation_stringency"
                        },
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_4_1_8_0/temporary_directory",
                            "source": "#uncollapsed_bam_generation.cwl/temporary_directory"
                        }
                    ],
                    "out": [
                        {
                            "id": "#uncollapsed_bam_generation.cwl/gatk_merge_sam_files_4_1_8_0/gatk_merge_sam_files_bam"
                        }
                    ],
                    "run": "#gatk_merge_sam_files_4.1.8.0.cwl",
                    "label": "GATK-MergeSamFiles",
                    "https://www.sevenbridges.com/x": 1039.828125,
                    "https://www.sevenbridges.com/y": 2528.5
                }
            ],
            "requirements": [
                {
                    "class": "SubworkflowFeatureRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "",
                    "https://schema.org/name": "Ronak Shah"
                }
            ],
            "https://schema.org/citation": "",
            "https://schema.org/codeRepository": "https://github.com/msk-access/uncollapsed_bam_generation",
            "https://schema.org/contributor": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/email": "mailto:shahr2@mskcc.org",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9042-6213",
                    "https://schema.org/name": "Ronak Shah"
                }
            ],
            "https://schema.org/dateCreated": "2020-09-23",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0"
        }
    ],
    "cwlVersion": "v1.0",
    "$schemas": [
        "http://schema.org/version/latest/schemaorg-current-http.rdf"
    ]
}