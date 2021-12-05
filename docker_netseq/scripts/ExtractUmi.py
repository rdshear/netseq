#!/usr/bin/env python3
import sys

def ExtractUmi(InputBam, OutputBam, umi_length, umi_tag = 'RX'):
    # input: bam file with 3' adapter code location at XT tag
    # output: RX tag added to bam file
    tag_header = '\t' + umi_tag + ':Z:'
    with open(InputBam, mode="r") as infile:
        with open(OutputBam, mode="w") as outfile:
            while True:
                x = infile.readline()
                if x == '':
                    break
                if x[0] != '@':
                    x = x[:-1] + tag_header + x.split('\t')[9][0:umi_length] + '\n'
                outfile.write(x)
    return

if __name__ == '__main__':
    if len(sys.argv) < 5:
        sys.exit("4 command line arguments expected")
    ExtractUmi(InputBam=sys.argv[1], OutputBam=sys.argv[2], 
        umi_length=int(sys.argv[3]), umi_tag=sys.argv[4])
