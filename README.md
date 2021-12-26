# NET-seq Pipeline

Process NET-seq short reads with UMIs in FASTQ format
## Alpha release. In active development

For S. cerevisiae. Currently hardwired to sacCer3
q`  
### Input

- Illumina reads - fastq.gz format
- Sample name - string
- Genome (fasta format)
- threads - number of threads to request for task containers

### Output
 
- bedgraph.gz occupancy counts
- bam files
- logs from STAR and deduplication task

### Acknowledgements

TBS
https://github.com/dockstore/checker-WDL-templates