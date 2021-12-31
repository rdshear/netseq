#! /bin/bash
cd ~/temp || exit
rm -R ~/temp/*

conda activate cpa

cromwell run -i ~/Projects/netseq/test/test-2/local_inputs_big.json -t wdl \
   -o /Users/robertshear/Projects/netseq/test/options.json \
   ~/Projects/netseq/netseq.wdl 