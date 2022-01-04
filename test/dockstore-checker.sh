#! /bin/bash
cd ~/temp || exit
rm -R ~/temp/*

conda activate cpa

cromwell run -i ~/Projects/netseq/test/test-2/checker_inputs.json -t wdl \
   -o /Users/robertshear/Projects/netseq/test/options.json \
   ~/Projects/netseq/checker_netseq.wdl 
