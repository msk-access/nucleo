---
description: Versions of tools in order of process
---

# Tools Description

| Tool | Version |
| :--- | :--- |
| [ProcessLoopUMIFastq ](https://github.com/mskcc/cwl-commandlinetools/tree/master/marianas_process_loop_umi_1.8.1)**\(Marianas\)** | **1.8.1** |
| \*\*\*\*[Trimgalore](https://github.com/mskcc/cwl-commandlinetools/tree/master/trim_galore_0.6.2) | **0.6.2** |
| \*\*\*\*[BWA mem](https://github.com/mskcc/cwl-commandlinetools/tree/master/bwa_mem_0.7.12) | **0.7.5a** |
| \*\*\*\*[AddOrReplaceReadGroups](https://github.com/mskcc/cwl-commandlinetools/tree/master/picard_add_or_replace_read_groups_1.96) **\(Picard\)** | **2.8.1** |
| \*\*\*\*[MarkDuplicates](https://github.com/mskcc/cwl-commandlinetools/tree/master/picard_mark_duplicates_2.8.1) **\(Picard\)** | **1.96** |
| \*\*\*\*[GenomeCov](https://github.com/mskcc/cwl-commandlinetools/tree/master/bedtools_genomecov_v2.28.0_cv2) **\(Bedtools\)** | **2.28.0\_cv2** |
| \*\*\*\*[Merge](https://github.com/mskcc/cwl-commandlinetools/tree/master/bedtools_merge_v2.28.0_cv2) **\(Bedtools\)** | **2.28.0\_cv2** |
| \*\*\*\*[ABRA](https://github.com/mskcc/cwl-commandlinetools/tree/master/abra2_2.17) | **2.17** |
| \*\*\*\*[FixMateInformation](https://github.com/mskcc/cwl-commandlinetools/tree/master/picard_fix_mate_information_1.96) **\(Picard\)** | **1.96** |
| \*\*\*\*[BaseRecalibrator](https://github.com/mskcc/cwl-commandlinetools/tree/master/gatk_BaseRecalibrator_4.1.2.0) **\(GATK\)** | **4.1.2.0** |
| \*\*\*\*[ApplyBQSR](https://github.com/mskcc/cwl-commandlinetools/tree/master/gatk_ApplyBQSR_4.1.2.0) **\(GATK\)** | **4.1.2.0** |
| [DuplexUMIBamToCollapsedFastqFirstPass](https://github.com/msk-access/cwl-commandlinetools/tree/master/marianas_collapsing_first_pass_1.8.1) **\(Marianas\)** | **1.8.1** |
| \*\*\*\*[DuplexUMIBamToCollapsedFastqSecondPass](https://github.com/msk-access/cwl-commandlinetools/tree/master/marianas_collapsing_second_pass_1.8.1) **\(Marianas\)** | **1.8.1** |
| \*\*\*\*[PileupMetrics](https://github.com/mskcc/cwl-commandlinetools/tree/master/waltz_pileupmatrices_3.1.1) **\(Waltz\)** | **1.0** |
| \*\*\*\*[CollectAlignmentSummaryMetrics](https://github.com/mskcc/cwl-commandlinetools/tree/develop/picard_collect_alignment_summary_metrics_2.8.1) **\(Picard\)** | **2.8.1** |

\*\*\*\*

UMI Processing:

* [Marianas](https://github.com/msk-access/cwl-commandlinetools/tree/master/marianas_process_loop_umi_1.8.1)

Standard Bam Processing Workflow:

* [Trimgalore](https://github.com/msk-access/cwl-commandlinetools/tree/master/trim_galore_0.6.2) 
* [BWA mem](https://github.com/msk-access/cwl-commandlinetools/tree/master/bwa_mem_0.7.12)
* [Picard AddOrReplaceReadGroups](https://github.com/msk-access/cwl-commandlinetools/tree/master/picard_add_or_replace_read_groups_1.96)
* [Picard MarkDuplicates](https://github.com/msk-access/cwl-commandlinetools/tree/master/picard_mark_duplicates_2.8.1)
* Bedtools
  * [GenomeCov](https://github.com/msk-access/cwl-commandlinetools/tree/master/bedtools_genomecov_v2.28.0_cv2)
  * [Merge](https://github.com/msk-access/cwl-commandlinetools/tree/master/bedtools_merge_v2.28.0_cv2)
* [ABRA](https://github.com/msk-access/cwl-commandlinetools/tree/master/abra2_2.17)
* [Picard FixMateInformation](https://github.com/msk-access/cwl-commandlinetools/tree/master/picard_fix_mate_information_1.96)
* GATK
  * [BaseRecalibrator](https://github.com/msk-accesss/cwl-commandlinetools/tree/master/gatk_BaseRecalibrator_4.1.2.0)
  * [ApplyBQSR](https://github.com/msk-access/cwl-commandlinetools/tree/master/gatk_ApplyBQSR_4.1.2.0)

Collapsing Workflow:

* [Marianas](https://github.com/msk-access/cwl-commandlinetools/tree/master/marianas_collapsing_first_pass_1.8.1) - 1.8.1
* [Waltz PileupMetrics](https://github.com/msk-access/cwl-commandlinetools/tree/master/waltz_pileupmatrices_3.1.1) - 1.0
* [BWA mem](https://github.com/msk-access/cwl-commandlinetools/tree/master/bwa_mem_0.7.5a) - 0.7.5a
* [Picard AddOrReplaceReadGroups](https://github.com/msk-access/cwl-commandlinetools/tree/master/picard_add_or_replace_read_groups_1.96) - 1.96
* [Picard CollectAlignmentSummaryMetrics](https://github.com/msk-access/cwl-commandlinetools/tree/develop/picard_collect_alignment_summary_metrics_2.8.1) - 2.8.1

