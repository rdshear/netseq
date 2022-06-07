#! /bin/bash
cd ~/temp || exit
rm -R ~/temp/*

conda activate cpa

# to check syntax 
# mamba run -n womtool womtool validate -i /Users/robertshear/Projects/netseq/test/local/inputs.json /Users/robertshear/Projects/netseq/netseq.wdl
cromwell run -i ~/Projects/netseq/test/local/inputs.json -t wdl \
    ~/Projects/netseq/netseq.wdl 
