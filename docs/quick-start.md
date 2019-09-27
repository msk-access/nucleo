---
description: Running the ACCESS standard bam workflow
---

# Quick Start

### Step 1 -  Clone the project from GitHub:

```text
$ git clone --recursive https://github.com/msk-access/fastq_to_bam.git
```

### Step 2 - Install Toil and cwltool

```text
$ pip install toil-ionox0'[cwl]'==0.0.7
```

Note: you may use any CWL executor available, provided it works with your batch system

### Step 3 - Generate an inputs file

Using either json or yaml format, following these specifications:

It's also possible to create a "template" inputs file to be filled in using this command:

```text
$ cwltool --make-template fastq_to_bam.cwl > inputs.yaml
```

### Step 4 - Run the workflow with a CWL executor:

```text
$ cwltool fastq_to_bam.cwl inputs.json
```

Or, to run the workflow on a HPC cluster using toil:

```text
$ toil-cwl-runner \
    --singularity \
    --jobStore my_jobStore \
    --batchSystem lsf \
    --workDir `pwd` \
    --outdir `pwd` \
    --logFile cwltoil.log \
    --writeLogs `pwd` \
    --logLevel DEBUG \
    --retryCount 2 \
    --disableCaching \
    --maxLogFileSize 20000000000 \
    --stats \
    fastq_to_bam.cwl \
    inputs.yaml
```

{% hint style="info" %}
**Note:** To see help for the inputs for any cwl workflow you can use:

$ cwltool fastq\_to\_bam.cwl --help
{% endhint %}

