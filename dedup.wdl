version 1.0

workflow dedup_wf {
        meta {
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }
parameter_meta {
        Infile: "Illumina Read file, FASTQ or sam/bam format"
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

        File Infile
#        Int maxSpotCount = 0
        Float total_reads = 0
        Float total_bases = 0
        String sampleName

        # environment
        String netseq_docker = 'rdshear/bbtools'
        Int preemptible = 1
        String memory = "12 GB"
        Int threads = 4
    }

    Int lclmem = ceil((total_reads * 500 + total_bases) / 1000000000.0) + 2

    call Dedup {
        input:
            Infile = Infile,
#            maxSpotCount = maxSpotCount,
            sampleName = sampleName,
            threads = threads,
            docker = netseq_docker,
            # per bbtools documentation, we need reads * 500 + total bases of memory
            memory = if total_bases + total_reads > 0 then "~{lclmem} GB" else memory,
            preemptible = preemptible
    }

    output {
        File output_fastq = Dedup.FastqDeduped
        File Log = Dedup.Log
    }
}

task Dedup {
    input {
        File Infile
#        Int maxSpotCount
        String sampleName
        Int threads = 4
        String docker
        Int preemptible
        String memory
    }

    String outfileName = "~{sampleName}.deduped.fastq.gz"
    String logfileName = "~{sampleName}.dedup.log"

    command <<<
        set -e
        echo 'calc memory=~{memory}'
        . /root/bbmap/dedupe.sh -eoom ac=f in=~{Infile} out=~{outfileName} 2> ~{logfileName}
    >>>

    output {
        File FastqDeduped = outfileName
        File Log = logfileName
    }

    runtime {
        docker: docker
        memory: memory
        cpu: threads
        disks: "local-disk 25 SSD"
        preemptible: preemptible
    }
}
