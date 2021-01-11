---
description: Input files and parameters required to run workflow
---

# Inputs Description

{% hint style="warning" %}
 Common workflow language execution engines accept two types of input that are [JSON](https://json.org) or [YAML](https://yaml.org), please make sure to use one of these while generating the input file. For more information refer to: [http://www.commonwl.org/user\_guide/yaml/](http://www.commonwl.org/user_guide/yaml/)
{% endhint %}

## Parameter Used by Tools
### Common Parameters Across Tools

| **Argument Name** | Summary | Default Value |
| :---: | :---: | :---: |
| **sequencing-center** | The sequencing center from which the data originated | MSKCC |
| **sample** | The name of the sequenced sample.\(**Required**\) |  |
| **run-date** | Date the run was produced, to insert into the read group header \(Iso8601Date\) |  |
| **read-group-id** | Read group ID to use in the file header \(**Required**\) |  |
| **platform-unit** | Read-Group Platform Unit \(eg. run barcode\) \(**Required**\) |  |
| **platform-model** | Platform model to insert into the group header \(ex. miseq, hiseq2500, hiseqX\) | novaseq |
| **platform** | Read-Group platform \(e.g. ILLUMINA, SOLID\). | ILLUMINA |
| **library** | The name/ID of the sequenced library. \(**Required**\) |  |
| **description** | Description of the read group. |  |
| **comment** | Comments to include in the output file’s header. |  |
| **validation\_stringency** | Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data \(read, qualities, tags\) do not otherwise need to be decoded. The --VALIDATION\_STRINGENCY argument is an enumerated type \(ValidationStringency\), which can have one of the following values: STRICT or LENIENT or SILENT | LENIENT |
| **sort\_order** | GATK: The order in which the reads should be output. |  |
| **create\_bam\_index** | GATK: Generate BAM index file when possible |  |
| **reference\_sequence** | Reference sequence file. Please include ".fai", "^.dict", ".amb" , ".sa", ".bwt", ".pac", ".ann" as secondary files if they are not present in the same location as the ".fasta" file |  |
| **temporary\_directory** | Temporary directory to be used for all steps |  |
| **fgbio\_async\_io** | Fgbio asynchronous execution |  |

### Uncollapsed BAM Generation

#### Fgbio [FastqToBam](https://github.com/msk-access/cwl-commandlinetools/tree/develop/fgbio_fastq_to_bam_1.2.0)

| **Argument Name** | **Summary** | **Default Value** |
| :---: | :---: | :---: |
| **fgbio\_fastq\_to\_bam\_umi-tag** | Tag in which to store molecular barcodes/UMIs. |  |
| **fgbio\_fastq\_to\_bam\_sort** | If true, query-name sort the BAM file, otherwise preserve input order. |  |
| **fgbio\_fastq\_to\_bam\_input** | Fastq files corresponding to each sequencing read \( e.g. R1, I1, etc.\). Please refer to the [template file](inputs-description.md#template-inputs-file) to get this correct.  |  |
| **read-structures** | Read structures, one for each of the FASTQs. Refer to the [tool ](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures)for more details |  |
| **fgbio\_fastq\_to\_bam\_predicted-insert-size** | Predicted median insert size, to insert into the read group header |  |
| **fgbio\_fastq\_to\_bam\_output\_file\_name** | The output SAM or BAM file to be written. |  |

#### Picard [MergeSamFiles](https://github.com/msk-access/cwl-commandlinetools/tree/develop/gatk_merge_sam_files_4.1.8.0)

| **Argument Name** | **Summary** | **Default Value** |
| :---: | :---: | :---: |
| **gatk\_merge\_sam\_files\_output\_file\_name** | SAM or BAM file to write the merged result to \(**Required**\) |  |
| **merge\_sam\_files\_sort\_order** | Sort order of output file | queryname |

#### Picard [SamToFastq](https://github.com/msk-access/cwl-commandlinetools/tree/develop/gatk_sam_to_fastq_4.1.8.0)

| **Argument Name** | **Summary** | **Default Value** |
| :---: | :---: | :---: |
| **unpaired\_fastq\_file** | unpaired fastq output file name |  |
| **UBG\_picard\_SamToFastq\_R1\_output\_fastq** | Read1 fastq.gz output file name for uncollapsed bam generation \(**Required**\) |  |
| **UBG\_picard\_SamToFastq\_R2\_output\_fastq** | Read2 fastq.gz output file name for uncollapsed bam generation \(**Required**\) |  |
| **BC\_gatk\_sam\_to\_fastq\_output\_name\_R1** | Read1 fastq.gz output file name for bam collapsing \(**Required**\) |  |
| **BC\_gatk\_sam\_to\_fastq\_output\_name\_R2** | Read2 fastq.gz output file name for bam collapsing \(**Required**\) |  |
| **gatk\_sam\_to\_fastq\_include\_non\_primary\_alignments** | If true, include non-primary alignments in the output. Support of non-primary alignments in SamToFastq is not comprehensive, so there may be exceptions if this is set to true and there are paired reads with non-primary alignments. |  |
| **gatk\_sam\_to\_fastq\_include\_non\_pf\_reads** | Include non-PF reads from the SAM file into the output FASTQ files. PF means 'passes filtering'. Reads whose 'not passing quality controls' flag is set are non-PF reads. See GATK Dictionary for more info. |  |

#### [Fastp](https://github.com/msk-access/cwl-commandlinetools/tree/develop/fastp_0.20.1)

| **Argument Name** | **Summary** | **Default Value** |
| :---: | :---: | :---: |
| **fastp\_unpaired1\_output\_file\_name** | For PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it. |  |
| **fastp\_unpaired2\_output\_file\_name** | For PE input, if read2 passed QC but read1 not, it will be written to unpaired2. If --unpaired2 is same as --unpaired1 \(default mode\), both unpaired reads will be written to this same file. |  |
| **fastp\_read1\_adapter\_sequence** | the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. | GATCGGAAGAGC |
| **fastp\_read2\_adapter\_sequence** | The adapter for read2 \(PE data only\). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as  \(string\) | AGATCGGAAGAGC |
| **fastp\_read1\_output\_file\_name** | Read1 output File Name \(**Required**\) |  |
| **fastp\_read2\_output\_file\_name** | Read2 output File Name \(**Required**\) |  |
| **fastp\_minimum\_read\_length** | reads shorter than length\_required will be discarded | 25 |
| **fastp\_json\_output\_file\_name** | the json format report file name \(**Required**\) |  |
| **fastp\_html\_output\_file\_name** | the html format report file name \(**Required**\) |  |
| **disable\_trim\_poly\_g** | Disable Poly-G trimming. | True |
| **disable\_quality\_filtering** | Disable base quality filtering. | True |

#### [BWA MEM](https://github.com/msk-access/cwl-commandlinetools/tree/develop/bwa_mem_0.7.17)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **bwa\_mem\_Y** | Force soft-clipping rather than default hard-clipping of supplementary alignments | True |
| **bwa\_mem\_T** | Don’t output alignment with score lower than INT. This option only affects output. | 30 |
| **bwa\_mem\_P** | In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair. |  |
| **UBG\_bwa\_mem\_output** | Output SAM file name for uncollapsed bam generation \(**Required**\) |  |
| **BC\_bwa\_mem\_output** | Output SAM file name for bam collapsing \(**Required**\) |  |
| **bwa\_mem\_M** | Mark shorter split hits as secondary |  |
| **bwa\_mem\_K** | to achieve deterministic alignment results \(Note: this is a hidden option\) | 1000000 |
| **bwa\_number\_of\_threads** | Number of threads |  |

#### Picard [AddOrReplaceReadGroups](https://github.com/msk-access/cwl-commandlinetools/tree/develop/picard_add_or_replace_read_groups_4.1.8.1)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **UBG\_picard\_addRG\_output\_file\_name** | Output BAM file name for uncollapsed bam generation \(**Required**\) |  |
| **BC\_picard\_addRG\_output\_file\_name** | Output BAM file name for bam collapsing \(**Required**\) |  |
| **picard\_addRG\_sort\_order** | Sort order for the BAM file | queryname |

#### GATK [MergeBamAlignment](https://github.com/msk-access/cwl-commandlinetools/tree/develop/gatk_merge_bam_alignment_4.1.8.0)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **UBG\_gatk\_merge\_bam\_alignment\_output\_file\_name** | Output BAM file name for uncollapsed bam generation \(**Required**\) |  |
| **BC\_gatk\_merge\_bam\_alignment\_output\_file\_name** | Output BAM file name for bam collapsing \(**Required**\) |  |

#### Picard [MarkDuplicates](https://github.com/msk-access/cwl-commandlinetools/tree/develop/picard_mark_duplicates_4.1.8.1)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **optical\_duplicate\_pixel\_distance** | The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is more appropriate. For other platforms and models, users should experiment to find what works best. | 2500 |
| **read\_name\_regex** | Regular expression that can be used to parse read names in the incoming SAM file. Read names are parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used to estimate the rate of optical duplication in order to give a more accurate estimated library size. Set this option to null to disable optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are extremely large and estimating library complexity is not an aim. Note that without optical duplicate counts, library size estimation will be inaccurate. The regular expression should contain three capture groups for the three variables, in order. It must match the entire read name. Note that if the default regex is specified, a regex match is not actually done, but instead the read name is split on colon character. For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. For 7 element names \(CASAVA 1.8\), the 5th, 6th, and 7th elements are assumed to be tile, x and y values. |  |
| **duplicate\_scoring\_strategy** | The scoring strategy for choosing the non-duplicate among candidates. |  |
| **gatk\_mark\_duplicates\_output\_file\_name** | The output file to write marked records to \(**Required**\) |  |
| **gatk\_mark\_duplicates\_duplication\_metrics\_file\_name** | File to write duplication metrics to \(**Required**\) |  |
| **gatk\_mark\_duplicates\_assume\_sort\_order** | If not null, assume that the input file has this order even if the header says otherwise. |  |

#### bedtools [genomecov](https://github.com/msk-access/cwl-commandlinetools/tree/develop/bedtools_genomecov_v2.28.0_cv2)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **bedtools\_genomecov\_option\_bedgraph** | option flag parameter to choose output file format. -bg refers to bedgraph format | True |

#### bedtools [merge](https://github.com/msk-access/cwl-commandlinetools/tree/develop/bedtools_merge_v2.28.0_cv2)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **bedtools\_merge\_distance\_between\_features** | Maximum distance between features allowed for features to be merged. | 10 |

#### [ABRA2](https://github.com/msk-access/cwl-commandlinetools/tree/develop/abra2_2.22)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **abra2\_window\_size** | Processing window size and overlap \(size,overlap\) | "400,200" |
| **abra2\_soft\_clip\_contig** | Soft clip contig args \[max_contigs,min\_base\_qual,frac_ high\_qual\_bases,min\_soft\_clip\_len\] | "16,13,80,15" |
| **abra2\_scoring\_gap\_alignments** | Scoring used for contig alignments\(match, mismatch\_penalty,gap\_open\_penalty,gap\_extend\_penalty\) | "8,32,48,1" |
| **abra2\_no\_sort** | Do not attempt to sort final output | True |
| **abra2\_no\_edge\_complex\_indel** | Prevent output of complex indels at read start or read end | True |
| **abra2\_maximum\_mixmatch\_rate** | Max allowed mismatch rate when mapping reads back to contigs | 0.1 |
| **abra2\_maximum\_average\_depth** | Regions with average depth exceeding this value will be down-sampled | 1000 |
| **abra2\_contig\_anchor** | Contig anchor \[M\_bases\_at\_contig\_edge,max\_mismatches\_near\_edge\] | "10,2" |
| **abra2\_consensus\_sequence** | Use positional consensus sequence when aligning high quality soft clipping |  |
| **BC\_abra2\_output\_bams** | The output BAM file to write to \(**Required**\) |  |
| **UBG\_abra2\_output\_bams** | The output BAM file to write to \(**Required**\) |  |
#### Picard [FixMateInformation](https://github.com/msk-access/cwl-commandlinetools/tree/develop/picard_fix_mate_information_4.1.8.1)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **UBG\_picard\_fixmateinformation\_output\_file\_name** | The output BAM file to write to for uncollapsed bam generation \(**Required**\) |  |
| **BC\_picard\_fixmate\_information\_output\_file\_name** | The output BAM file to write to for bam collapsing \(**Required**\) |  |

### Base Quality Score Recalibration
#### GATK [BaseRecalibrator](https://github.com/msk-access/cwl-commandlinetools/tree/develop/gatk_base_recalibrator_4.1.8.1)	

| **Argument Name** | **Summary** | **Default Value** |	
| :---: | :--- | :---: |	
| **gatk\_base\_recalibrator\_known\_sites** | One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis \(**Required**\) |  |	
| **gatk\_bqsr\_read\_filter** | Read filters to be applied before analysis |  |	
| **base\_recalibrator\_output\_file\_name** | The output recalibration table file to create \(**Required**\) |  |	

#### GATK [ApplyBQSR](https://github.com/msk-access/cwl-commandlinetools/tree/develop/gatk_apply_bqsr_4.1.8.1)	

| **Argument Name** | **Summary** | **Default Value** |	
| :---: | :--- | :---: |	
| **apply\_bqsr\_output\_file\_name** | The output BAM file \(**Required**\) |  |	
| **gatk\_bqsr\_disable\_read\_filte**r | Read filters to be disabled before analysis |  |

### Collapsed BAM Generation

#### Fgbio [GroupReadsByUmi](https://github.com/msk-access/cwl-commandlinetools/tree/develop/fgbio_group_reads_by_umi_1.2.0)

| **Argument Name** | Summary | Default Value |
| :---: | :---: | :---: |
| **fgbio\_group\_reads\_by\_umi\_input** | The input BAM file |  |
| **fgbio\_group\_reads\_by\_umi\_strategy** | The UMI assignment strategy. \(identity, edit, adjacency, paired\) | paired |
| **fgbio\_group\_reads\_by\_umi\_raw\_tag** | The tag containing the raw UMI. | RX |
| **fgbio\_group\_reads\_by\_umi\_output\_file\_name** | The output BAM file name \(**Required**\) |  |
| **fgbio\_group\_reads\_by\_umi\_min\_umi\_length** | The minimum UMI length. If not specified then all UMIs must have the same length, otherwise, discard reads with UMIs shorter than this length and allow for differing UMI lengths. |  |
| **fgbio\_group\_reads\_by\_umi\_include\_non\_pf\_reads** | Include non-PF reads. | False |
| **fgbio\_group\_reads\_by\_umi\_family\_size\_histogram** | Optional output of tag family size counts. | True |
| **fgbio\_group\_reads\_by\_umi\_edits** | The allowable number of edits between UMIs. | 1 |
| **fgbio\_group\_reads\_by\_umi\_assign\_tag** | The output tag for UMI grouping. | MI |

#### Fgbio [CollectDuplexSeqMetrics](https://github.com/msk-access/cwl-commandlinetools/tree/develop/fgbio_collect_duplex_seq_metrics_1.2.0)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **fgbio\_collect\_duplex\_seq\_metrics\_intervals** | Optional set of intervals over which to restrict analysis. |  |
| **fgbio\_collect\_duplex\_seq\_metrics\_output\_prefix** | Prefix of output files to write. |  |
| **fgbio\_collect\_duplex\_seq\_metrics\_min\_ba\_reads** | Minimum BA reads to call a tag family a ‘duplex’. |  |
| **fgbio\_collect\_duplex\_seq\_metrics\_min\_ab\_reads** | Minimum AB reads to call a tag family a ‘duplex’. |  |
| **fgbio\_collect\_duplex\_seq\_metrics\_mi\_tag** | The output tag for UMI grouping. | MI |
| **fgbio\_collect\_duplex\_seq\_metrics\_duplex\_umi\_counts** | If true, produce the .duplex\_umi\_counts.txt file with counts of duplex UMI observations. | True |
| **fgbio\_collect\_duplex\_seq\_metrics\_description** | Description of data set used to label plots. Defaults to sample/library. |  |

#### Fgbio [CallDuplexConsensusReads](https://github.com/msk-access/cwl-commandlinetools/tree/develop/fgbio_call_duplex_consensus_reads_1.2.0)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **fgbio\_call\_duplex\_consensus\_reads\_trim** | If true, quality trim input reads in addition to masking low Q bases. |  |
| **fgbio\_call\_duplex\_consensus\_reads\_sort\_order** | The sort order of the output, if :none: then the same as the input. |  |
| **fgbio\_call\_duplex\_consensus\_reads\_read\_name\_prefix** | The prefix all consensus read names |  |
| **fgbio\_call\_duplex\_consensus\_reads\_read\_group\_id** | The new read group ID for all the consensus reads. |  |
| **fgbio\_call\_duplex\_consensus\_reads\_output\_file\_name** | Output SAM or BAM file to write consensus reads. |  |
| **fgbio\_call\_duplex\_consensus\_reads\_min\_reads** | The minimum number of input reads to a consensus read. | 1 1 0 |
| **fgbio\_call\_duplex\_consensus\_reads\_min\_input\_base\_quality** | Ignore bases in raw reads that have Q below this value. |  |
| **fgbio\_call\_duplex\_consensus\_reads\_max\_reads\_per\_strand** | The maximum number of reads to use when building a single-strand consensus. If more than this many reads are present in a tag family, the family is randomly downsampled to exactly max-reads reads. |  |
| **fgbio\_call\_duplex\_consensus\_reads\_error\_rate\_pre\_umi** | The Phred-scaled error rate for an error prior to the UMIs being integrated. |  |
| **fgbio\_call\_duplex\_consensus\_reads\_error\_rate\_post\_umi** | The Phred-scaled error rate for an error post the UMIs have been integrated. |  |

#### Fgbio [FilterConsensusReads](https://github.com/msk-access/cwl-commandlinetools/tree/develop/fgbio_filter_consensus_reads_1.2.0)

<table>
  <thead>
    <tr>
      <th style="text-align:center">Argument Name</th>
      <th style="text-align:center">Summary</th>
      <th style="text-align:center">Default Value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex</b>
      </td>
      <td style="text-align:center">Reverse [complement] per base tags on reverse strand reads.- Simplex+Duplex</td>
      <td
      style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_reverse_per_base_tags_duplex</b>
      </td>
      <td style="text-align:center">Reverse [complement] per base tags on reverse strand reads. - Duplex</td>
      <td
      style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_require_single_strand_agreement_simplex_duplex</b>
      </td>
      <td style="text-align:center">Mask (make N) consensus bases where the AB and BA consensus reads disagree
        (for duplex-sequencing only).</td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center">f<b>gbio_filter_consensus_read_require_single_strand_agreement_duplex</b>
      </td>
      <td style="text-align:center">Mask (make N) consensus bases where the AB and BA consensus reads disagree
        (for duplex-sequencing only).</td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_max_base_error_rate_duplex</b>
      </td>
      <td style="text-align:center">The maximum error rate for a single consensus base. (Max 3 values) - Duplex</td>
      <td
      style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_max_base_error_rate_simplex_duplex</b>
      </td>
      <td style="text-align:center">The maximum error rate for a single consensus base. (Max 3 values) - Simplex
        + Duplex</td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_max_no_call_fraction_duplex</b>
      </td>
      <td style="text-align:center">Maximum fraction of no- calls in the read after filtering - Duplex</td>
      <td
      style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_max_read_error_rate_duplex</b>
      </td>
      <td style="text-align:center">
        <p></p>
        <p>The maximum raw-read error rate across the entire consensus read. (Max
          3 values) - Duplex</p>
      </td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_max_no_call_fraction_simplex_duplex</b>
      </td>
      <td style="text-align:center">Maximum fraction of no- calls in the read after filtering - Simplex +
        Duplex</td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_max_read_error_rate_simplex_duplex</b>
      </td>
      <td style="text-align:center">The maximum raw-read error rate across the entire consensus read. (Max
        3 values) - Simplex + Duplex</td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_min_base_quality_duplex</b>
      </td>
      <td style="text-align:center">Mask (make N) consensus bases with quality less than this threshold. -
        Dupelx</td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_min_base_quality_simplex_duplex</b>
      </td>
      <td style="text-align:center">
        <p></p>
        <p>Mask (make N) consensus bases with quality less than this threshold. -
          Simplex+Dupelx</p>
      </td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_min_mean_base_quality_duplex</b>
      </td>
      <td style="text-align:center">The minimum mean base quality across the consensus read - Duplex</td>
      <td
      style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_min_mean_base_quality_simplex_duplex</b>
      </td>
      <td style="text-align:center">The minimum mean base quality across the consensus read - Simplex + Duplex</td>
      <td
      style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_min_reads_duplex</b>
      </td>
      <td style="text-align:center">
        <p></p>
        <p>The minimum number of reads supporting a consensus base/read. (Max 3 values)
          - Duplex</p>
      </td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_min_reads_simplex_duplex</b>
      </td>
      <td style="text-align:center">
        <p></p>
        <p>The minimum number of reads supporting a consensus base/read. (Max 3 values)</p>
        <p>-Simplex+Duplex</p>
      </td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_output_file_name_simplex_duplex</b>
      </td>
      <td style="text-align:center">Output BAM file name Simplex + Duplex</td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics</b>
      </td>
      <td style="text-align:center">Output file name Duplex alignment metrics</td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics</b>
      </td>
      <td style="text-align:center">Output file name Simplex alignment metrics</td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_output_file_name_duplex</b>
      </td>
      <td style="text-align:center">Output BAM file name - Duplex</td>
      <td style="text-align:center"></td>
    </tr>
    <tr>
      <td style="text-align:center"><b>fgbio_filter_consensus_read_min_simplex_reads</b>
      </td>
      <td style="text-align:center">
        <p></p>
        <p>The minimum number of reads supporting a consensus base/read. (Max 3 values)
          -</p>
        <p>Simplex+Duplex</p>
      </td>
      <td style="text-align:center"></td>
    </tr>
  </tbody>
</table>

#### Fgbio [Postprocessing](https://github.com/msk-access/cwl-commandlinetools/tree/develop/fgbio_postprocessing_simplex_filter_0.1.8)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **fgbio\_postprocessing\_output\_file\_name\_simplex** | Output BAM file name Simplex |  |

#### Picard [CollectAlignmentSummaryMetrics](https://github.com/mskcc/cwl-commandlinetools/tree/develop/picard_collect_alignment_summary_metrics_2.8.1)

| Argument Name | Summary | Default Value |
| :---: | :---: | :---: |
| **gatk\_collect\_alignment\_summary\_metrics\_output\_file\_name** | Output file name for metrics on collapsed BAM \(Duplex+Simplex+Singletons\) |  |

## Template Inputs File

{% code title="inputs.yaml" %}
```bash
BC_abra2_output_bams: null
BC_bwa_mem_output: null
BC_gatk_merge_bam_alignment_output_file_name: null
BC_gatk_sam_to_fastq_output_name_R1: null
BC_gatk_sam_to_fastq_output_name_R2: null
BC_picard_addRG_output_file_name: null
BC_picard_fixmate_information_output_file_name: null
UBG_abra2_output_bams: null
UBG_bwa_mem_output: null
UBG_gatk_merge_bam_alignment_output_file_name: null
UBG_picard_SamToFastq_R1_output_fastq: null
UBG_picard_SamToFastq_R2_output_fastq: null
UBG_picard_addRG_output_file_name: null
UBG_picard_fixmateinformation_output_file_name: null
abra2_bam_index: null
abra2_consensus_sequence: null
abra2_contig_anchor: null
abra2_maximum_average_depth: null
abra2_maximum_mixmatch_rate: null
abra2_no_edge_complex_indel: null
abra2_scoring_gap_alignments: null
abra2_soft_clip_contig: null
abra2_window_size: null
apply_bqsr_output_file_name: null
base_recalibrator_output_file_name: null
bedtools_genomecov_option_bedgraph: null
bedtools_merge_distance_between_features: null
bwa_mem_K: null
bwa_mem_T: null
bwa_mem_Y: null
create_bam_index: null
fastp_html_output_file_name: null
fastp_json_output_file_name: null
fastp_minimum_read_length: null
fastp_read1_adapter_sequence: null
fastp_read1_output_file_name: null
fastp_read2_adapter_sequence: null
fastp_read2_output_file_name: null
fgbio_async_io: null
fgbio_call_duplex_consensus_reads_min_reads: null
fgbio_call_duplex_consensus_reads_output_file_name: null
fgbio_collect_duplex_seq_metrics_duplex_umi_counts: null
fgbio_collect_duplex_seq_metrics_intervals: null
fgbio_collect_duplex_seq_metrics_output_prefix: null
fgbio_fastq_to_bam_input: null
fgbio_filter_consensus_read_min_base_quality_duplex: null
fgbio_filter_consensus_read_min_base_quality_simplex_duplex: null
fgbio_filter_consensus_read_min_reads_duplex: null
fgbio_filter_consensus_read_min_reads_simplex_duplex: null
fgbio_filter_consensus_read_output_file_name_duplex: null
fgbio_filter_consensus_read_output_file_name_duplex_aln_metrics: null
fgbio_filter_consensus_read_output_file_name_simplex_aln_metrics: null
fgbio_filter_consensus_read_output_file_name_simplex_duplex: null
fgbio_filter_consensus_read_reverse_per_base_tags_simplex_duplex: null
fgbio_group_reads_by_umi_family_size_histogram: null
fgbio_group_reads_by_umi_output_file_name: null
fgbio_group_reads_by_umi_strategy: null
fgbio_postprocessing_output_file_name_simplex: null
gatk_base_recalibrator_add_output_sam_program_record: null
gatk_base_recalibrator_known_sites:
  - class: File
    metadata: {}
    path: >-
      /Users/shahr2/Documents/test_reference/test_fastq_to_bam/known_sites/dbsnp_137_14_16.b37.vcf
    secondaryFiles:
      - class: File
        path: >-
          /Users/shahr2/Documents/test_reference/test_nucleo/known_sites/dbsnp_137_14_16.b37.vcf.idx
  - class: File
    metadata: {}
    path: >-
      /Users/shahr2/Documents/test_reference/test_fastq_to_bam/known_sites/Mills_and_1000G_gold_standard-14_16.indels.b37.vcf
    secondaryFiles:
      - class: File
        path: >-
          /Users/shahr2/Documents/test_reference/test_fastq_to_bam/known_sites/Mills_and_1000G_gold_standard-14_16.indels.b37.vcf.idx
gatk_collect_alignment_summary_metrics_output_file_name: null
gatk_mark_duplicates_duplication_metrics_file_name: null
gatk_mark_duplicates_output_file_name: null
gatk_merge_sam_files_output_file_name: null
library: null
merge_sam_files_sort_order: null
optical_duplicate_pixel_distance: null
picard_addRG_sort_order: null
platform: null
platform-model: null
platform-unit: null
read-group-id: null
read-structures: null
reference_sequence:
  class: File
  metadata: {}
  path: /Users/shahr2/Documents/test_reference/fasta/chr14_chr16.fasta
  secondaryFiles:
    - class: File
      path: ../../test_reference/fasta/chr14_chr16.fasta.amb
    - class: File
      path: ../../test_reference/fasta/chr14_chr16.fasta.ann
run-date: null
sample: null
sequencing-center: null
sort_order: null
temporary_directory: null
validation_stringency: null
```

