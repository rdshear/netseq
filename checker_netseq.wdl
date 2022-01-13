version 1.0

import "netseq.wdl" as netseq_wf

workflow checker_netseq {
        meta {
        description: "Checker for rdshear/netseq workflow"
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

    input {
        String srrId = "SRR12840066"
        Int maxReads = 1000
        Int limitedReads = 50
        File truth_md5
    }


    call getSamples {
        input:
            srrId = srrId,
            maxReads = maxReads,
            limitedReads = limitedReads
    }
    # call netseq_wf.netseq as subject {
    #     input:
    #         sampleName = sampleName,
    #         inputFastQ = inputFastQ
    # }

    # call md5_filecheck {
    #     input:
    #         sampleName = sampleName,
    #         output_bam = subject.output_bam,
    #         bedgraph_neg = subject.bedgraph_neg,
    #         bedgraph_pos = subject.bedgraph_pos,
    #         truth_md5 = truth_md5
    # }
    
}

task getSamples {
    input {
        String srrId
        Int maxReads
        Int limitedReads
    }
    command <<<
        set -e

        # make sure that miniconda is properly initialized whether interactive or not
        . /bin/entrypoint.sh
        fastq-dump -X 1000 ~{srrId} --stdout  > testsample.fastq
        head -n ~{limitedReads*4} testsample.fastq > smallsample.fastq
        md5sum testsample.fastq smallsample.fastq > test.md5
        gzip testsample.fastq
        gzip smallsample.fastq
    >>>
    output {
        File fastqSample = "testsample.fastq.gz"
        File fastqShortSample = "smallsample.fastq.gz"
        File testmd5 = "test.md5"
    }
    runtime {
        docker: 'rdshear/netseq'
    }
}




task md5_filecheck {
    input {
        String sampleName
        File output_bam
        File bedgraph_neg
        File bedgraph_pos
        File truth_md5
    }

    command <<<
        set -e

        # make sure that miniconda is properly initialized whether interactive or not
        . /bin/entrypoint.sh

        samtools view ~{output_bam} > ~{sampleName}.sam
        gunzip -c ~{bedgraph_neg} > ~{sampleName}.neg.bedgraph
        gunzip -c ~{bedgraph_pos} > ~{sampleName}.pos.bedgraph

        md5sum  -c ~{truth_md5}
    >>>

    runtime {
        docker: 'rdshear/netseq'
    }
}