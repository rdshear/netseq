version 1.0

import "netseq.wdl" as netseq_wf

workflow checker_netseq {
        meta {
        description: "Quick check for quay.io/rdshear/netseq workflow"
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

    input {
        String refFasta
        String sampleName
        File inputFastQ
        File truth_md5
    }


    call netseq_wf.netseq as subject {
        input:
            refFasta = refFasta,
            inputFastQ = inputFastQ
    }

    call md5_filecheck {
        input:
            sampleName = sampleName,
            output_bam = subject.output_bam,
            bedgraph_neg = subject.bedgraph_neg,
            bedgraph_pos = subject.bedgraph_pos,
            truth_md5 = truth_md5
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
        . /usr/local/bin/_entrypoint.sh

        samtools view ~{output_bam} > ~{sampleName}.sam
        gunzip -c ~{bedgraph_neg} > ~{sampleName}.neg.bedgraph
        gunzip -c ~{bedgraph_pos} > ~{sampleName}.pos.bedgraph

        md5sum  -c ~{truth_md5}
    >>>

    runtime {
        docker: 'quay.io/rdshear/netseq'
    }
}