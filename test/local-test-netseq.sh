#! /bin/bash
cd ~/temp || exit
rm -R ~/temp/*

conda activate cpa

cromwell run -i ~/Projects/netseq/test/inputs.json -t wdl \
    ~/Projects/netseq/netseq.wdl 
#    -o /Users/robertshear/Projects/netseq/test/options.json \
#    ~/Projects/netseq/netseq.wdl 