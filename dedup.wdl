version 1.0

workflow dedup {
        meta {
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

parameter_meta {
        inputFastQ: "Illumina Read file, FASTQ format"
#        maxSpotCount: "If not zero, then maximum number of fastQ record to read"

        #Outputs

        output_fastq: "input file deduped"

        # Environment
        netseq_docker: "Name of docker image"
        threads: "Number of CPUs to request for task runtimes"
        memory: "Memory required, in format appropriate to platform. Default is 8G"
        preemptible: "For terra.bio (GCP), Number of retries in premptable mode before running in on-demand mode. If 1, then attempt to run in preemptbile mode once"
    }
    input {

        File inputFastQ
#        Int maxSpotCount = 0

        # environment
        String netseq_docker = 'rdshear/bbtools'
        Int preemptible = 1
        # TODO...calculate it!
        String memory = "8G"
        Int threads = 8
    }



    call Dedup {
        input:
            Infile = inputFastQ,
#            maxSpotCount = maxSpotCount,

            threads = threads,
            docker = netseq_docker,
            memory = memory,
            preemptible = preemptible
    }

    output {
        File output_fastq = Dedup.FastqDeduped
        # TODO: Add log fle
    }
}

task Dedup {
    input {
        File Infile
#        Int maxSpotCount

        Int threads = 8
        String docker
        Int preemptible
        String memory
    }


    command <<<
        set -e

        . /root/bbmap/dedupe.sh -eoom ac=f in=~{Infile} out="todo.fastq.gz"
    >>>

    output {
        File FastqDeduped = "todo.fastq.gz"
    }

    runtime {
        docker: docker
        memory: memory
        cpu: threads
        disks: "local-disk 25 SSD"
        preemptible: preemptible
    }
}
