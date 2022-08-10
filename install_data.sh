#!/usr/bin/env bash


# Test data is hosted on Google Drive at:
# https://drive.google.com/file/d/1GtT8jsBGwRoQC-5wHh06r8RFkiFBuirp/view?usp=sharing

fileid=1GtT8jsBGwRoQC-5wHh06r8RFkiFBuirp

filename=test_nucleo.tar.gz

pip install gdown
gdown https://drive.google.com/uc?id=$fileid 

# Suppress linux warnings for MacOS tar.gz files
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    tar --warning=no-unknown-keyword -xzvf $filename 
elif [[ "$OSTYPE" == "darwin"* ]]; then
    tar -xzvf $filename 
fi

rm $filename
