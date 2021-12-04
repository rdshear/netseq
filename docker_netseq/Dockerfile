# /bin/bash -c docker build . -t rdshear/netseq
#############################################################
# Net-seq allignment (fastq to bedgraph)
#############################################################

FROM continuumio/miniconda3:latest
LABEL description="Components for NET-seq pipeline"

WORKDIR /root

# Install the needed apps
COPY environment.yml ./environment.yml
RUN conda env update -f environment.yml

# execution time scripts
COPY scripts/ /home/netseq/scripts/

CMD ["/bin/bash"]
