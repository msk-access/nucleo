#!/usr/bin/env bash


# Test data is hosted on Google Drive at:
# https://drive.google.com/open?id=1aQ5Hkm6XV7fk1qm1lTM5k-UYfG66KjXi

fileid=1aQ5Hkm6XV7fk1qm1lTM5k-UYfG66KjXi

filename=test_fastq_to_bam.tar.gz
foldername=test_fastq_to_bam

# Skip if already have test data
[[ -f $filename ]] && exit 0
[[ -d $foldername ]] && exit 0

curl -c ./cookie -s -k -L "https://drive.google.com/uc?export=download&id=$fileid" > /dev/null

curl -k -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}

# Suppress linux warnings for MacOS tar.gz files
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    tar --warning=no-unknown-keyword -xzvf $filename
elif [[ "$OSTYPE" == "darwin"* ]]; then
    tar -xzvf $filename
fi

rm $filename