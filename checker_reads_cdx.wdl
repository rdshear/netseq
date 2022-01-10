version 1.0

import "reads_cdx.wdl" as wf

workflow checker_reads_cdx {
        meta {
        description: "Unit tests for reads_cdx.wdl"
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

    input {
        String sraId
        File truth_md5
    }


    call wf.reads_cdx as fq_all {
        input:
            sraId = sraId,
            OutputFileName = "SRR_full.fastq.gz"
    }

    call wf.reads_cdx as fq_100 {
        input:
            sraId = sraId,
            OutputFileName = "SRR_100.fastq.gz",
            MaxReadCount = 100
    }

    call wf.reads_cdx as bam_all {
        input:
            sraId = sraId,
            resultAsUbam = true,
            OutputFileName = "SRR_full.bam",
    }

    call wf.reads_cdx as bam_100 {
        input:
            sraId = sraId,
            resultAsUbam = true,
            OutputFileName = "SRR_100.bam",
            MaxReadCount = 100
    }


    call md5_reads_cdxfilecheck {
        input:
            sraFiles = [fq_all.OutputFile, fq_100.OutputFile, bam_all.OutputFile, bam_100.OutputFile],
            truth_md5 = truth_md5
    }
    
}

task md5_reads_cdxfilecheck {
    input {
        Array[File] sraFiles
        File truth_md5
    }

    command <<<
        set -e

        # make sure that miniconda is properly initialized whether interactive or not
        . /bin/entrypoint.sh
        cp ~{sep=" " sraFiles} .
        md5sum -c ~{truth_md5}
    >>>

    runtime {
        docker: 'rdshear/netseq'
    }
}