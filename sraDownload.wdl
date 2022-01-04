version 1.0

task sraDownload {
    input {
        String sraId

        Int threads = 2
        String docker = "rdshear/sra-tools"
        Int preemptible = 1
        String memory = "2G"
    }
    
    String uBamFilename = "${sraId}.bam"
    command <<<
        sam-dump --gzip ~{sraId} > ~{uBamFilename}
    >>>

    output {
        File ubam_file = uBamFilename
    }

    runtime {
        docker: docker
        memory: memory
        cpu: threads
        disks: "local-disk 25 SSD"
        preemptible: preemptible
    }
}

workflow sraDownloadWorkflow {
    input {
        String sraId
    }
    call sraDownload {
        input:
            sraId = sraId
    }
}
