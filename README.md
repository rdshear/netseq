# NET-seq Pipeline

The rdshear/netseq workflow transforms the raw sequence reads from the NET-seq assay to occupancy counts.

## Overview

NET-seq (native elongating transcript sequencing) reports the strand-specific density of active RNA elongation density at single nucleotide resolution (Churchman2011-me). 
The protocol results in a library of strand-specific short cDNA fragments, the 5'-end of which corresponds to the location of active elongation at the moment the sample was harvested.

The input to this pipeline is single-ended sequencing reads

![FASTQ sequence schematic](./fastq_schematic.drawio.png)

The pipeline takes as input a single-end FASTQ file, optionally gz compressed. 
Alternatively, given an [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) run identifier, the pipeline will automatically download the file from SRA.

The primary output of the pipeline is a pair of bedgraph files carrying the count at each nucleotide position.

### Processing Steps

1. [sra-tools/fastp-dump](https://github.com/ncbi/sra-tools) (optional)
    * Read public SRA archive data
2. [fastp](https://github.com/OpenGene/fastp)
    * Report on characteristics of FASTQ file
    * Remove reads of poor quality
    * Trim the 3'-end adapter
    * Trim the UMI sequence at the 5'-end
    * Remove duplicates (amplicons and optical duplicates)
3. [STAR](https://github.com/alexdobin/STAR)
    * Read and pre-process the reference genome
    * Align to genome
    * Remove ambiguous alignments ("multimappers)
    * Generate coordinate-sorted BAM file
4. [bedtools](https://github.com/arq5x/bedtools2/)
    * Create the bedgraph files

### Example inputs

```
{

}
```
### Time and cost estimates

|Sample Name|Reads|Bases|Duration|Cost|
|-----------|-----|-----|--------|----|
|SRR12840066|4.03G|52.97M|26 min|$0.02|

### Limitations

* This workflow has only been 
    * tested on terra.bio
    * aligned against the sacCer3 genome

* Only tested on sacCer3.
* Only tested on terra.bio
* Steps 1-3 are "piped". Requires extra memory

## Technical notes
### Deduplication approach

Dedup before alignment reduces possibility of interference by aligner
### STAR aligner configuration

No splicing. 44 nt average read length after trimming. Chance of a nascent co-transcrptionally spliced read is very low at best.

[sra-tools gs problem!](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-access-costs/)
## Detailed Instructions
### Input

|Parameter|Type|Description|Default Value|
|---------|----|-----------|-------------|
|||**Required Parameters**||
|refFasta|String|Uri to genome geference gile, FASTA format|[^1]|
|genomeName|String|UCSC genome version name|sacCer3|
|inputFastQ|File|The reads to be processed in fastq format Extension must be ".fastq", ".fastq.gz" or ".fastq.gz.1"||
|sraRunId|String|An SRA run identifier, e.g. SRR12840066. If inputFastQ is absent, then sraRunId must be present. If both are present, inputFastQ takes presidence.||
|sampleName|String|An identifier for the sample in a format compatible with file names.|[^2]|
|||**Preprocessing Parameters**||
|maxReadCount|Int|If defined and greater than zero, then only the first _n_ reads will be processed. Useful for testing, but not for downsampling|0|
|adapterSequence|String|The 3'-end adapter sequence|ATCTCGTATGCCGTCTTCTGCTTG|
|umiWidth|Int|The number of nucleotides in the UMI. Zero indicates that no UMI is present, in which case no duplication will be performed|6|
|||**STAR Aligner Parameters**||
|outSAMmultNmax|Int|Default to outputting primary alignment only. The number of alignments returned for each read. If 1, then no multimmappers are returned|1|
|OutFilterMultiMax|Int|If a read has multimappers in excess of this paramter, then the read is disreagarded|10|
|||**Environmental Parameters**||
|netseq_docker|String|Override the name of the Docker image invoked by the workflow.|rdshear/netseq|
|preemptible|Int|Applicable to Google cloud (GCP). The default value, 1, indicates that each task will be attempted first on a [spot instance](https://cloud.google.com/spot-vms). If spot instance is preemted, then it will be run on a standard on-demand instance   see [Cromwell Documentation](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#preemptible) for details. *Caution* setting the preemptible parameter to zero is likely to dramatically incrase the cost running the workflow.|1|
memory|String|The amount of RAM to use. The value is [dependent on the backend platform](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#memory). For [terra.bio]("https://terra.bio), the value _n_G indicates that it least _n_GiB are required to run the workflow |8G|
|threads|Int|The number of cpus requested to run the workflow|4|

[^1] https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz

[^2] For example, if sampleName is "wt-1", then the plus strand bedgraph output will have the filename "wt-1.pos.bedgraph.gz".|Will be constructed from the inputFastQ parameter file base name if defined, _viz._ "wt-1" from inputFastQ "wt-1.fastq.gz". If inputFastQ is not defined, then the sraRunId, _viz._ "SRR12840066"|

## Output

|Parameter|Type|Format|Example|
|---------|----|-----------|-------|
|output_bam|File|BAM|wt-1.bam|
|bedgraph_pos|File|bedgraph (+) strand|wt-1.pos.bedgraph.gz|
|bedgraph_neg|File|bedgraph (-) strand|wt-1.neg.bedgraph.gz|
|alignment_log|File|text|wt-1.Log.final.out|
|fastp_report_html|File|html report|wt-1.fastp.html|
|fastp_report_json|File|json|wt-1.fastp.json|

### Running the workflow

Example workflow configurations...

####Test Data

#### Steps

