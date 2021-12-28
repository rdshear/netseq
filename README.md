# NET-seq Pipeline

Process NET-seq short reads with UMIs in FASTQ format
## Alpha release. In active development

For S. cerevisiae. Currently hardwired to sacCer3
### Input

- Illumina reads - fastq.gz format or SRA Id

- Sample name - string
- Genome (fasta format)
- threads - number of threads to request for task containers

### Output
 
- bedgraph.gz occupancy counts
- bam files
- logs from STAR and deduplication task

