version 1.0

import "sraRead.wdl" as wf

workflow checker_sraRead {
        meta {
        description: "Unit tests for sraRead.wdl"
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

    input {
        String sraId
        File truth_md5
    }


    call wf.sraRead as fq_all {
        input:
            sraId = sraId,
            OutputFileName = "SRR_full.fastq.gz"
    }

    call wf.sraRead as fq_100 {
        input:
            sraId = sraId,
            OutputFileName = "SRR_100.fastq.gz",
            MaxReadCount = 100
    }

    call wf.sraRead as bam_all {
        input:
            sraId = sraId,
            resultAsUbam = true,
            OutputFileName = "SRR_full.bam",
    }

    call wf.sraRead as bam_100 {
        input:
            sraId = sraId,
            resultAsUbam = true,
            OutputFileName = "SRR_100.bam",
            MaxReadCount = 100
    }


    call md5_sraReadfilecheck {
        input:
            sraFiles = [fq_all.OutputFile, fq_100.OutputFile, bam_all.OutputFile, bam_100.OutputFile],
            truth_md5 = truth_md5
    }
    
}

task md5_sraReadfilecheck {
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