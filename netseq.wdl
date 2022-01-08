version 1.0

import "dedup.wdl" as dedup_wf
import "sraRead.wdl" as sraRead_wf

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
        File? inputFastQ
        String? sraRunId
        Float total_reads = 0
        Float total_bases = 0

        Int? maxReadCount
#        String sampleName = if !defined(inputFastQ) then sraRunId else basename(basename(select_first([inputFastQ]), ".gz"), ".fastq")
        String inputFileFullName = select_first([inputFastQ, sraRunId, 'default'])
        String sampleName = basename(basename(inputFileFullName, ".gz"), ".fastq")
        String adapterSequence = "ATCTCGTATGCCGTCTTCTGCTTG"
        Int umiWidth = 6
        String umi_tag = "RX"

        # environment
        String netseq_docker = 'rdshear/netseq'
        Int preemptible = 1
        String memory = "8G"
        Int threads = 8
    }

    if (!defined(inputFastQ)) {
        call sraRead_wf.sraRead as sra {
            input: 
                sraId = select_first([sraRunId]),
                maxReadCount = maxReadCount,
                OutputFileName = sampleName + ".fastq.gz"
        }
    }

    if (umiWidth > 0) {
        call dedup_wf.dedup as dedupeResult {
            input:
                Infile = select_first([inputFastQ, sra.OutputFile]),
                sampleName = sampleName,
                total_reads = total_reads,
                total_bases = total_bases
        }
    }

    call AlignReads {
        input:
            Infile = select_first([dedupeResult.FastqDeduped, inputFastQ]),
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
        File output_bam = AlignReads.BamFileDeduped
        File bedgraph_pos = AlignReads.CoverageBedgraph_Pos
        File bedgraph_neg = AlignReads.CoverageBedgraph_Neg
        File? dedup_log = dedupeResult.Log
        File alignment_log = AlignReads.star_log_final
    }
}

task AlignReads {
    input {
        File? Infile
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

    String bamDedupName = "~{sampleName}.bam"

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
 
        # force the temp directory to the container's disks
        tempStarDir=$(mktemp -d)
        # star wants to create the directory itself
        rmdir "$tempStarDir"

        samtools import -0 ~{Infile} \
        | if [[ ~{maxReadCount} -ne 0 ]]
            then 
                head -n ~{maxReadCount} 
            else
                cat
        fi \
        | if [[ ~{umiWidth} -ne 0 ]]
            then
                python3 /home/micromamba/scripts/ExtractUmi.py \
                    /dev/stdin /dev/stdout ~{umiWidth} ~{umi_tag} 
            else
                cat
        fi \
        | STAR --runMode alignReads \
            --genomeDir star_work \
            --runThreadN ~{threads} \
            --readFilesIn /dev/stdin \
            --readFilesCommand samtools view \
            --readFilesType SAM SE \
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
        | samtools sort >  ~{bamDedupName}

        bedtools genomecov -5 -bg -strand - -ibam ~{bamDedupName} | bgzip > ~{sampleName}.pos.bedgraph.gz
        bedtools genomecov -5 -bg -strand + -ibam ~{bamDedupName} | bgzip > ~{sampleName}.neg.bedgraph.gz
    >>>

    output {
        File star_log_final = "~{sampleName}.Log.final.out"
        File star_log = "~{sampleName}.Log.out"
        File star_log_std =  "~{sampleName}.Log.std.out"
        File BamFileDeduped = bamDedupName
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
