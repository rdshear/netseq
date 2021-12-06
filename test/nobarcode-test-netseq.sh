#! /bin/bash
cd ~/temp || exit
rm -R ~/temp/*

conda activate cpa

cromwell run -i /Users/robertshear/Projects/netseq/test/inputs_SRR072814.json -t wdl \
   -o /Users/robertshear/Projects/netseq/test/options.json \
   ~/Projects/netseq/netseq.wdl 