version 1.0

task reads_cdx {
    input {
        String? sraId
        File? inputFastQ
        Int? maxReadCount
        String OutputFileName
        
        
        # technical
        Int threads = 4
        String docker = "rdshear/sra-tools"
        Int preemptible = 1
        String memory = "G"
    }


    String sraCommand = if resultAsUbam then 
                            (if defined(maxReadCount) then "sam-dump ~{sraId} | head -n ~{maxReadCount} | gzip -c"
                                else "sam-dump --gzip ~{sraId}")
                        else ("fastq-dump --stdout " + (if defined(maxReadCount) then "-X " + maxReadCount + " " else "")
                                    + sraId + " | gzip -c")

    command <<<
        set -e
        ~{sraCommand} > ~{OutputFileName}
    >>>

    output {
        File OutputFile = OutputFileName
    }

    runtime {
        docker: docker
        memory: memory
        cpu: threads
        disks: "local-disk 25 SSD"
        preemptible: preemptible
    }
}

workflow reads_cdxWorkflow {
    input {
        String sraId
        String OutputFileName
        Boolean resultAsUbam = false
        Int? maxReadCount
    }
    call reads_cdx {
        input:
            sraId = sraId,
            OutputFileName = OutputFileName,
            resultAsUbam = resultAsUbam,
            maxReadCount = maxReadCount
    }

    output {
        File OutputFile = reads_cdx.OutputFile
    }
}
