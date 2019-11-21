# Inputs Description

### UMI Processing Parameters:

| **Argument Name** | **Summary** | **Default Value** |
| :--- | :--- | :--- |
| **umi\_length** | Length of 5' UMI to use for collapsing | 3 |
| **fastq1** | Read 1 of the paired-end sample |  |
| **fastq2** | Read 2 of the paired-end sample |  |



### Standard Bam Processing Parameters:

| **Argument Name** | **Summary** | **Default Value** |
| :--- | :--- | :--- |
| **adapter** | Adapter sequence to be trimmed | GATCGGAAGAGC |
| **adapter2** | Optional adapter sequence to be trimmed off read 2 of paired-end files. This option requires '--paired' to be specified as well | AGATCGGAAGAGC |
| **quality** | Trim low-quality ends from reads in addition to adapter removal. For RRBS samples, quality trimming will be performed first, and adapter trimming is carried in a second round. Other files are quality and adapter trimmed in a single pass. The algorithm is the same as the one used by BWA \(Subtract INT from all qualities; compute partial sums from all indices to the end of the sequence; cut sequence at the index at which the sum is minimal\). | 1 |
| **stringency** | Overlap with adapter sequence required to trim a sequence. Defaults to a very stringent setting of '1', i.e. even a single bp of overlapping sequence will be trimmed of the 3' end of any read. | 3 |
| **length** | Discard reads that became shorter than length INT because of either quality or adapter trimming. A value of '0' effectively disables this behaviour. | 25 |
| **P** | Skip pairing for BWA; mate rescue performed unless -S also in use | true |
| **M** | Mark shorter split hits as secondary \(for Picard/GATK compatibility\) | true |
| **reference** | Reference fasta file |  |
| **read\_group\_identifier** | Read Group ID  Default value: 1. This option can be set to 'null' to clear the default value  Required |  |
| **read\_group\_sequencing\_center** | Read Group sequencing center name  Default value: null. Required. | MSKCC |
| **read\_group\_library** | Read Group Library.  Required |  |
| **read\_group\_platform\_unit** | Read Group platform unit \(eg. run barcode\)  Required. |  |
| **read\_group\_sequencing\_platform** | Read Group platform \(e.g. illumina, solid\)  Required. | ILLUMINA |
| **maximum\_average\_depth** | Regions with average depth exceeding this value will be downsampled | 1000 |
| **soft\_clip\_contig** | Soft clip contig args \[max\_contigs,min\_base\_qual,frac\_high\_qual\_bases,min\_soft\_clip\_len \] | 100,30,80,15 |
| **maximum\_mixmatch\_rate** | Max allowed mismatch rate when mapping reads back to contigs | 0.1 |
| **scoring\_gap\_alignments** | Scoring used for contig alignments \(match,mismatch\_penalty,gap\_open\_penalty,gap\_extend\_penalty\) | 8,32,48,1 |
| **contig\_anchor** | Contig anchor \[M\_bases\_at\_contig\_edge,max\_mismatches\_near\_edge\] | 10,1 |
| **window\_size** | Processing window size and overlap \(size,overlap\) | 800,700 |
| **consensus\_sequence** | Use positional consensus sequence when aligning high quality soft clipping | true |
| **output\_file\_name** | Required list of output sam or bam file \(s\) separated by comma |  |
| **ignore\_bad\_assembly** | Use this option to avoid parsing errors for corrupted assemblies | true |
| **bam\_index** | Enable BAM index generation when outputting sorted alignments \(may require additional memory\) | true |
| **no\_sort** | Do not attempt to sort final output |  |
| **known\_sites\_1 & known\_sites\_2** | One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis |  |
| **read\_filter** | Read filters to be applied before analysis | GoodCigarReadFilter |



### Bam Collapsing Parameters:

<table>
  <thead>
    <tr>
      <th style="text-align:left"><b>Argument Name</b>
      </th>
      <th style="text-align:left"><b>Summary</b>
      </th>
      <th style="text-align:left"><b>Default Value</b>
      </th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="text-align:left"><b>bed_file</b>
      </td>
      <td style="text-align:left">The genotype from positions in this bed file will be used as the consensus
        base if <b>min_consensus_percent</b> threshold is not reached. Otherwise,
        the reference base from the supplied <b>reference_fasta</b> will be used</td>
      <td
      style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>wobble</b>
      </td>
      <td style="text-align:left">Allowable left and right shift amount for grouping UMI families</td>
      <td
      style="text-align:left">1</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>reference_fasta</b>
      </td>
      <td style="text-align:left">Currently, hg19</td>
      <td style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>mismatches</b>
      </td>
      <td style="text-align:left">Allowable mismatch count in UMI bases for grouping UMI families</td>
      <td
      style="text-align:left">0</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>min_consensus_percent</b>
      </td>
      <td style="text-align:left">Percentage of bases that must be in agreement at each position in the
        consensus read before masking that base as &quot;N&quot;</td>
      <td style="text-align:left">90</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>min_base_quality</b>
      </td>
      <td style="text-align:left">? used for masking as &quot;N&quot;?</td>
      <td style="text-align:left">20</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>key</b>
      </td>
      <td style="text-align:left">SHOULD BE A DEFAULT!</td>
      <td style="text-align:left">
        <p></p>
        <ul>
          <li>6,6n</li>
          <li>8,8n</li>
        </ul>
      </td>
    </tr>
    <tr>
      <td style="text-align:left"><b>output_name_collapsed_gzip_R1</b>
      </td>
      <td style="text-align:left">
        <p></p>
        <p>Filename for collapsed read 1 fastq</p>
      </td>
      <td style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>output_name_collapsed_gzip_R2</b>
      </td>
      <td style="text-align:left">Filename for collapsed read 2 fastq</td>
      <td style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>picard_output_file_name</b>
      </td>
      <td style="text-align:left">Filename to be given to final re-aligned bam (SHOULD have DEFAULT?)</td>
      <td
      style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>aln_output_file_name</b>
      </td>
      <td style="text-align:left">Filename to be given to intermediate sam file (SHOULD have DEFAULT?)</td>
      <td
      style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>alignment_metrics_unfiltered</b>
      </td>
      <td style="text-align:left">Filename for alignment metrics TXT file generated by Picard CollectAlignmentMetrics
        for Unfilered BAM File.</td>
      <td style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>alignment_metrics_simplex</b>
      </td>
      <td style="text-align:left">Filename for alignment metrics TXT file generated by Picard CollectAlignmentMetrics
        for SIMPLEX BAM File.</td>
      <td style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>alignment_metrics_duplex</b>
      </td>
      <td style="text-align:left">Filename for alignment metrics TXT file generated by Picard CollectAlignmentMetrics
        for DUPLEX BAM File.</td>
      <td style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>assume_sorted</b>
      </td>
      <td style="text-align:left">Assume that the given bam file is coordinate sorted for picard tools</td>
      <td
      style="text-align:left">true</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>bqsr_read_filter</b>
      </td>
      <td style="text-align:left">GATK READ_FILTER option to apply defferent set of ReadFilter</td>
      <td style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>collapsing_aln_output_file_name</b>
      </td>
      <td style="text-align:left">Name of the SAM format output file created by bwa mem for collapsing step.</td>
      <td
      style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>collapsing_picard_output_file_name</b>
      </td>
      <td style="text-align:left">Name of the BAM format output file created by Picard AddOrReplaceReadGroups
        for collapsing step.</td>
      <td style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>consensus_sequence</b>
      </td>
      <td style="text-align:left">Use positional consensus sequence when aligning high quality soft clipping</td>
      <td
      style="text-align:left">true</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>create_bam_index</b>
      </td>
      <td style="text-align:left"></td>
      <td style="text-align:left">true</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>number_of_threads</b>
      </td>
      <td style="text-align:left">Number of threads for parallel execution of ABRA</td>
      <td style="text-align:left">16</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>sort_order</b>
      </td>
      <td style="text-align:left">How the BAM file should be sorted</td>
      <td style="text-align:left">coordinate</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>trim_galore_number_of_threads</b>
      </td>
      <td style="text-align:left">Number of threads to run Trim Galore with Cutadapt</td>
      <td style="text-align:left">4</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>validation_stringency</b>
      </td>
      <td style="text-align:left">Picard Validation Stringency while running Picard Tools</td>
      <td style="text-align:left">LENIENT</td>
    </tr>
    <tr>
      <td style="text-align:left"><b>option_bedgraph</b>
      </td>
      <td style="text-align:left"></td>
      <td style="text-align:left">true</td>
    </tr>
  </tbody>
</table>

### Parameters used across multiple steps:

<table>
  <thead>
    <tr>
      <th style="text-align:left">Argument Name</th>
      <th style="text-align:left">Summary</th>
      <th style="text-align:left">Default Value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="text-align:left"><b>reference_fasta</b>
      </td>
      <td style="text-align:left">
        <p></p>
        <ul>
          <li>Used in the follwing steps:
            <ul>
              <li>Waltz PileupMetrics</li>
              <li>Marianas Collapsing</li>
              <li>Re-alignment of collapsed fastqs with BWA mem</li>
              <li>Indel Realignment with Abra</li>
              <li>Picard CollectAlignmentSummaryMetrics</li>
            </ul>
          </li>
        </ul>
      </td>
      <td style="text-align:left"></td>
    </tr>
    <tr>
      <td style="text-align:left"><b>min_map_quality</b>
      </td>
      <td style="text-align:left">Used for masking as &quot;N&quot;?</td>
      <td style="text-align:left">1</td>
    </tr>
  </tbody>
</table>



### Example Inputs File:

```text
{
    "standard_aln_output_file_name": "test_standard.sam",
    "standard_picard_addrg_output_filename": "test_standard.bam",
    "collapsing_aln_output_file_name": "test_unfiltered.sam",
    "collapsing_picard_output_file_name": "test_unfiltered.bam",
    "output_name_collapsed_gzip_R1": "test_collapsed_R1_fastq.gz",
    "output_name_collapsed_gzip_R2": "test_collapsed_R2_fastq.gz",
    "sort_first_pass_output_file_name": "test_collapsing_first_pass.txt",
    
    "read_group_identifier": "test",
    "read_group_library": 1,
    "read_group_platform_unit": "IDT11",
    "read_group_sample_name": "test",
    
    "fastq1": {
        "class": "File",
        "path": "/path/to/test_R1_001.fastq.gz",
    },
    "fastq2": {
        "class": "File",
        "path": "/path/to/test_R2_001.fastq.gz",
    },
    "bed_file": {
        "class": "File",
        "path": "/path/to/panel.bed",
    },
    "known_sites_1": {
        "class": "File",
        "path": "/path/to/dbsnp_137_14_16.b37.vcf",
    },
    "known_sites_2": {
        "class": "File",
        "path": "/path/to/Mills_and_1000G_gold_standard-14_16.indels.b37.vcf",
    },
    "reference": {
        "class": "File",
        "path": "/path/to/Homo_sapiens_assembly19.fasta",
    },
}
```

