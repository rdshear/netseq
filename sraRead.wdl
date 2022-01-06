version 1.0

task sraRead {
    input {
        String sraId
        Boolean resultAsUbam = false
        Int? MaxReadCount
        String OutputFileName
        
        
        # technical
        Int threads = 2
        String docker = "rdshear/sra-tools"
        Int preemptible = 1
        String memory = "2G"
    }


    String sraCommand = if resultAsUbam then 
                            (if defined(MaxReadCount) then "sam-dump ~{sraId} | head -n ~{MaxReadCount} | gzip -c"
                                else "sam-dump --gzip ~{sraId}")
                        else ("fastq-dump --stdout " + (if defined(MaxReadCount) then "-X " + MaxReadCount + " " else "")
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

workflow sraReadWorkflow {
    input {
        String sraId
        String OutputFileName
        Boolean resultAsUbam = false
        Int? MaxReadCount
    }
    call sraRead {
        input:
            sraId = sraId,
            OutputFileName = OutputFileName,
            resultAsUbam = resultAsUbam,
            MaxReadCount = MaxReadCount
    }

    output {
        File OutputFile = sraRead.OutputFile
    }
}
