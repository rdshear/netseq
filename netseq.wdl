version 1.0

workflow netseq {
        meta {
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

parameter_meta {
        # STAR index input
        refFasta: "Url to genome reference gile, FASTA format"
        genomeName: "Geneome name, UCSC version name, default sacCer3"

        # main input paramters
        inputFastQ: "Illumina Read file, FASTQ format"
        sampleName: "Sample name. If not specified, taken as base name of fastq input file"
        maxReadCount: "If defined, then maximum number of fastQ record to read"

        # STAR alignment parameters
        adapterSequence: "Adapter sequence to trim from 3' end"
        umiWidth: "Number of bases in UMI. Defaults to 6. If zero, no UMI deduplication occurs"
        outSAMmultNmax: "The number of alignments returned for each read. If 1, then no multimmappers are returned."
        outFilterMultiMax: "If a read has multimappers in excess of this paramter, then the read is disreagarded. Defaults"

        #Outputs

        output_bam: "aligned, deduped BAM faile"
        bedgraph_pos: "Occupancy counts on + strand, bedgraph format"
        bedgraph_neg: "Occupancy counts on - strand, bedgraph format"
        mask_pos: "Coverage of multimappers on + strand, bedgraph format"
        mask_neg: "Coverage of multimappers on - strand, bedgraph format"

        # Environment
        netseq_docker: "Name of docker image"
        threads: "Number of CPUs to request for task runtimes"
        memory: "Memory required, in format appropriate to platform. Default is 8G"
        preemptible: "For terra.bio (GCP), Number of retries in premptable mode before running in on-demand mode. If 1, then attempt to run in preemptbile mode once"
    }
    input {

        # Genome source for STAR
        String refFasta = "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz"
        String genomeName = "sacCer3"
        Int outSAMmultNmax = 1     # Default to outputting primary alignment only.
        Int outFilterMultiMax = 1   # Default to dropping reads with more than 1 alignment, implies ignore multi-mappers

        # Unprocessed reads
        File? inputFastQ
        String? sraRunId

        Int maxReadCount = 0
        String sampleName = basename(basename(select_first([inputFastQ, sraRunId, 'default']), ".gz"), ".fastq")
        String adapterSequence = "ATCTCGTATGCCGTCTTCTGCTTG"
        Int umiWidth = 6

        # environment
        # TODO: version the image
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
            outSAMmultNmax = outSAMmultNmax,
            outFilterMultiMax = outFilterMultiMax,
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
        File mask_pos = AlignReads.MaskBedggraph_Pos
        File mask_neg = AlignReads.MaskBedggraph_Neg
        File alignment_log = AlignReads.star_log_final
        File fastp_report_html = AlignReads.fastp_report_html
        File fastp_report_json = AlignReads.fastp_report_json
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
        Int outSAMmultNmax
        Int outFilterMultiMax
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

        # if pulling https://site/path/infile.fa.gz, then unzip it before storing it.
        lclname=$(basename ~{refFasta})
        extension=${lclname##*.}
        unzipFasta=$(if [[ ${extension} == "gz" ]]
            then 
                echo "gzip --decompress - "
            else
                echo "cat"
        fi)
        wget --quiet ~{refFasta} -O - | ${unzipFasta} > ~{genomeName}.fa
        samtools faidx ./~{genomeName}.fa

        STAR \
        --runMode genomeGenerate \
        --runThreadN ~{threads} \
        --genomeDir star_work \
        --genomeFastaFiles ./~{genomeName}.fa \
        --genomeSAindexNbases 10
 
        # force the temp directory to the container's local disks
        tempStarDir=$(mktemp -d)
        # star wants to create the directory itself
        rmdir "$tempStarDir"
        
        cmd=$(if [[ ~{UseFile} == true ]]
            then echo -n fastp -i ~{Infile} --reads_to_process ~{maxReadCount}
            else echo -n fastq-dump $(if [[ ~{maxReadCount} -gt 0 ]]; then echo -n -X ~{maxReadCount}; fi) --stdout ~{sraRunId} "|" fastp --stdin
            fi)
#TODO condition "dedup" processing on the presence of the umi   
        cmd=" $cmd --stdout -D --dup_calc_accuracy ~{DupCalcAccuracy} \
            --adapter_sequence ~{adapterSequence} \
            --umi --umi_len ~{umiWidth} --umi_loc per_read \
            --umi_prefix umi \
            --html ~{sampleName}.fastp.html \
            --json ~{sampleName}.fastp.json \
            --report_title \"~{sampleName} fastp report\""

        echo "Computed fastq pull command: $cmd" 

        eval "$cmd" | STAR --runMode alignReads \
            --genomeDir star_work \
            --runThreadN ~{threads} \
            --readFilesIn /dev/stdin  \
            --outTmpDir "$tempStarDir" \
            --outFileNamePrefix ~{sampleName}. \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped None \
            --outSAMmultNmax ~{outSAMmultNmax} \
            --outFilterMultimapNmax ~{outFilterMultiMax} \
            --clip3pNbases 0 \
            --clip5pNbases 0 \
            --limitBAMsortRAM 3221225472 \
            --alignIntronMax 1
        
        mv ~{sampleName}.Aligned.sortedByCoord.out.bam ~{bamFileName}
 
        tmpBAM=$(mktemp)
        # create occupancy bedGraphs
        samtools view -h --expr "[NH]==1" ~{bamFileName} -b > $tmpBAM
        bedtools genomecov -5 -bg -strand - -ibam $tmpBAM | bgzip > ~{sampleName}.pos.bedgraph.gz
        bedtools genomecov -5 -bg -strand + -ibam $tmpBAM | bgzip > ~{sampleName}.neg.bedgraph.gz

        # create multimapper coverage (for masking) as bedGraphs
        samtools view -h --expr "[NH]>1" ~{bamFileName} -b > $tmpBAM
        bedtools genomecov -bg -strand - -ibam $tmpBAM | bgzip > ~{sampleName}.mask_pos.bedgraph.gz
        bedtools genomecov -bg -strand + -ibam $tmpBAM | bgzip > ~{sampleName}.mask_neg.bedgraph.gz

        rm $tmpBAM
    >>>

    output {
        File star_log_final = "~{sampleName}.Log.final.out"
        File fastp_report_html = "~{sampleName}.fastp.html"
        File fastp_report_json = "~{sampleName}.fastp.json"
        File BamFile = bamFileName
        File CoverageBedgraph_Pos = '~{sampleName}.pos.bedgraph.gz'
        File CoverageBedgraph_Neg = '~{sampleName}.neg.bedgraph.gz'
        File MaskBedggraph_Pos = '~{sampleName}.mask_pos.bedgraph.gz'
        File MaskBedggraph_Neg = '~{sampleName}.mask_neg.bedgraph.gz'
    }

    runtime {
        docker: docker
        memory: memory
        cpu: threads
        disks: "local-disk 25 SSD"
        preemptible: preemptible
    }
}
