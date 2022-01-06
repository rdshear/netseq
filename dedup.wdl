version 1.0

workflow dedup_wf {
        meta {
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }
parameter_meta {
        Infile: "Illumina Read file, FASTQ or sam/bam format"
#        maxReadCount: "If not zero, then maximum number of fastQ record to read"

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
#        Int maxReadCount = 0
        Float total_reads = 0
        Float total_bases = 0
        String sampleName

        # environment
        String netseq_docker = 'rdshear/bbtools'
        Int preemptible = 1
        Int threads = 4
    }

    call Dedup {
        input:
            Infile = Infile,
#            maxReadCount = maxReadCount,
            total_reads = total_reads,
            total_bases = total_bases,
            sampleName = sampleName,
            threads = threads,
            docker = netseq_docker,
            # per bbtools documentation, we need reads * 500 + total bases of memory
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
#        Int maxReadCount
        Float total_reads = 0
        Float total_bases = 0
        String sampleName
        Int threads = 4
        String docker
        Int preemptible
    }

    # memory required for bbtools dedupe (per documentation)
    Float bb_memory = if total_reads + total_bases == 0 then 2000000000 else total_reads * 500.0 + total_bases

    String java_memory = ceil(bb_memory / 1000000.0) + "m"
    String docker_memory = ceil(bb_memory / 0.85 / 1000000000) + 2 + " GB"

    String outfileName = "~{sampleName}.deduped.fastq.gz"
    String logfileName = "~{sampleName}.dedup.log"

    command <<<
        set -e
        echo 'calc java memory=~{java_memory} docker memory=~{docker_memory}'
        . /root/bbmap/dedupe.sh -Xmx~{java_memory} -eoom ac=f in=~{Infile} out=~{outfileName} 2> ~{logfileName}
    >>>

    output {
        File FastqDeduped = outfileName
        File Log = logfileName
    }

    runtime {
        docker: docker
        memory: docker_memory
        cpu: threads
        disks: "local-disk 25 SSD"
        preemptible: preemptible
    }
}
