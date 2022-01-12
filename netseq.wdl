version 1.0

workflow netseq {
        meta {
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

parameter_meta {
        # STAR index input
        refFasta: "Url to genome geference gile, FASTA format"
        genomeName: "Geneome name, UCSC version name, default sacCer3"

        # main input paramters
        inputFastQ: "Illumina Read file, FASTQ format"
        sampleName: "Sample name. If not specified, taken as base name of fastq input file"
        maxReadCount: "If defined, then maximum number of fastQ record to read"

        # STAR alignment parameters
        adapterSequence: "Adapter sequence to trim from 3' end"
        umiWidth: "Number of bases in UMI. Defaults to 6. If zero, no UMI deduplication occurs"

        #Outputs

        output_bam: "aligned, deduped BAM faile"
        bedgraph_pos: "Occupancy counts on + strand, bedgraph format"
        bedgraph_neg: "Occupancy counts on - strand, bedgraph format"

        # Environment
        netseq_docker: "Name of docker image"
        threads: "Number of CPUs to request for task runtimes"
        memory: "Memory required, in format appropriate to platform. Default is 8G"
        preemptible: "For terra.bio (GCP), Number of retries in premptable mode before running in on-demand mode. If 1, then attempt to run in preemptbile mode once"
    }
    input {

        # Genome source for STAR
        String refFasta
        String genomeName = "sacCer3"

        # Unprocessed reads
        File? inputFastQ
        String? sraRunId

        Int maxReadCount = 0
        String sampleName = basename(basename(select_first([inputFastQ, sraRunId, 'default']), ".gz"), ".fastq")
        String adapterSequence = "ATCTCGTATGCCGTCTTCTGCTTG"
        Int umiWidth = 6

        # environment
        String netseq_docker = 'rdshear/netseq'
        Int preemptible = 1
        String memory = "8G"
        Int threads = 4
    }

    call AlignReads {
        input:
            Infile = inputFastQ,
            sraRunId = sraRunId,
            refFasta = refFasta,
            sampleName = sampleName,
            genomeName = genomeName,
            maxReadCount = maxReadCount,
            adapterSequence = adapterSequence,
            umiWidth = umiWidth,
            threads = threads,
            docker = netseq_docker,
            memory = memory,
            preemptible = preemptible
    }

    output {
        File output_bam = AlignReads.BamFile
        File bedgraph_pos = AlignReads.CoverageBedgraph_Pos
        File bedgraph_neg = AlignReads.CoverageBedgraph_Neg
        File fastp_json = AlignReads.fastp_report_json
        File fastp_report_html = AlignReads.fastp_report_html
    }
}

task AlignReads {
    input {
        File? Infile
        String? sraRunId
        Int maxReadCount
        String refFasta
        String genomeName
        String sampleName
        String adapterSequence
        Int umiWidth

        Int DupCalcAccuracy = 3 # TODO DESCRIBE
        Int threads = 8
        String docker
        Int preemptible
        String memory
    }

    String bamFileName = "~{sampleName}.bam"

    Boolean UseFile = defined(Infile)

    command <<<
        set -e

        # make sure that miniconda is properly initialized whether interactive or not
        . /bin/entrypoint.sh

        wget --quiet ~{refFasta} -O ~{genomeName}.fa
        bwa index ~{genomeName}.fa

        
        cmd=$(if [[ ~{UseFile} == true ]]
            then echo -n fastp -i ~{Infile} --reads_to_process ~{maxReadCount}
            else echo -n fastq-dump $(if [[ ~{maxReadCount} -gt 0 ]]; then echo -n -X ~{maxReadCount}; fi) --stdout ~{sraRunId} "|" fastp --stdin
            fi)

        cmd=" $cmd -o ~{sampleName}.dedup.fastq.gz -D --dup_calc_accuracy ~{DupCalcAccuracy} \
            --adapter_sequence ~{adapterSequence} \
            --umi --umi_len ~{umiWidth} --umi_loc per_read \
            --umi_prefix umi \
            --html ~{sampleName}.fastp.html \
            --json ~{sampleName}.fastp.json"

        echo "Computed fastq pull command: $cmd" 

        eval "$cmd" 
        bwa aln sacCer3.fa ~{sampleName}.dedup.fastq.gz | bwa samse sacCer3.fa - ~{sampleName}.dedup.fastq.gz | samtools sort -O BAM > ~{bamFileName}

        bedtools genomecov -5 -bg -strand - -ibam ~{bamFileName} | bgzip > ~{sampleName}.pos.bedgraph.gz
        bedtools genomecov -5 -bg -strand + -ibam ~{bamFileName} | bgzip > ~{sampleName}.neg.bedgraph.gz
    >>>

    output {
        File fastp_report_html = "~{sampleName}.fastp.html"
        File fastp_report_json = "~{sampleName}.fastp.json"
        File BamFile = bamFileName
        File CoverageBedgraph_Pos = '~{sampleName}.pos.bedgraph.gz'
        File CoverageBedgraph_Neg = '~{sampleName}.neg.bedgraph.gz'
    }

    runtime {
        docker: docker
        memory: memory
        cpu: threads
        disks: "local-disk 25 SSD"
        preemptible: preemptible
    }
}
