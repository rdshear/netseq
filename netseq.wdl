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

        # STAR alignment parameters
        inputFastQ: "Illumina Read file, FASTQ format"
        maxReadCount: "If defined, then maximum number of fastQ record to read"
        sampleName: "Sample name. If not specified, taken as base name of fastq input file"
        adapterSequence: "Adapter sequence to trim from 3' end"
        umiWidth: "Number of bases in UMI. Defaults to 6. If zero, no UMI deduplication occurs"
        umi_tag: "SAM tag for the UMI. Defaults to RX, which is the standard SAM recommended tag"
        outSAMmultNmax: "The number of alignments returned for each read. If 1, then no multimmappers are returned."
        OutFilterMultiMax: "If a read has multimappers in excess of this paramter, then the read is disreagarded. Defaults"

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
        Int outSAMmultNmax = 1     # Default to outputting primary alignment only. (Multimap count still available)
        Int OutFilterMultiMax = 10 # Default to dropping reads with more than 10 alignments

        # Unprocessed reads
        File inputFastQ
        String? sraRunId
        Float total_reads = 0
        Float total_bases = 0

        Int? maxReadCount
        String sampleName = basename(basename(select_first([inputFastQ, sraRunId, 'default']), ".gz"), ".fastq")
        String adapterSequence = "ATCTCGTATGCCGTCTTCTGCTTG"
        Int umiWidth = 6
        String umi_tag = "RX"

        # environment
        String netseq_docker = 'rdshear/netseq'
        Int preemptible = 1
        String memory = "8G"
        Int threads = 8
    }

    call AlignReads {
        input:
            Infile = inputFastQ,
            refFasta = refFasta,
            sampleName = sampleName,
            genomeName = genomeName,
            outSAMmultNmax = outSAMmultNmax,
            OutFilterMultiMax = OutFilterMultiMax,
            adapterSequence = adapterSequence,
            umiWidth = umiWidth,
            umi_tag = umi_tag,
            threads = threads,
            docker = netseq_docker,
            memory = memory,
            preemptible = preemptible
    }

    output {
        File output_bam = AlignReads.BamFile
        File bedgraph_pos = AlignReads.CoverageBedgraph_Pos
        File bedgraph_neg = AlignReads.CoverageBedgraph_Neg
        File alignment_log = AlignReads.star_log_final
    }
}

task AlignReads {
    input {
        File Infile
        Int maxReadCount = 0
        String refFasta
        String genomeName
        String sampleName
        Int outSAMmultNmax
        Int OutFilterMultiMax
        String adapterSequence
        Int umiWidth
        String umi_tag

        Int threads = 8
        String docker
        Int preemptible
        String memory
    }

    String bamFileName = "~{sampleName}.bam"

    command <<<
        set -e

        # make sure that miniconda is properly initialized whether interactive or not
        . /bin/entrypoint.sh

        wget --quiet ~{refFasta} -O ~{genomeName}.fa

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
        
        fastp -i ~{Infile} -o ~{sampleName}.fastq.gz -D --dup_calc_accuracy 4 --umi --umi_len ~{umiWidth} --umi_loc per_read \
            --umi_prefix umi --html ~{sampleName}.html --adapter_sequence ~{adapterSequence}
        STAR --runMode alignReads \
            --genomeDir star_work \
            --runThreadN ~{threads} \
            --readFilesIn /dev/stdin  \
            --outTmpDir "$tempStarDir" \
            --outStd SAM \
            --outFileNamePrefix ~{sampleName}. \
            --outReadsUnmapped None \
            --outSAMmultNmax ~{outSAMmultNmax} \
            --outFilterMultimapNmax ~{OutFilterMultiMax} \
            --clip3pAdapterSeq ~{adapterSequence} \
            --clip3pNbases 0 \
            --clip5pNbases ~{umiWidth}  \
            --alignIntronMax 1 \
        | samtools sort >  ~{bamFileName}

        bedtools genomecov -5 -bg -strand - -ibam ~{bamFileName} | bgzip > ~{sampleName}.pos.bedgraph.gz
        bedtools genomecov -5 -bg -strand + -ibam ~{bamFileName} | bgzip > ~{sampleName}.neg.bedgraph.gz
    >>>

    output {
        File star_log_final = "~{sampleName}.Log.final.out"
        File fastp_report = "~{sampleName}.html"
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
