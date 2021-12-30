# /bin/bash -c docker build . -t rdshear/netseq
#############################################################
# Net-seq allignment (fastq to bedgraph)
#############################################################

FROM mambaorg/micromamba:0.18.2
LABEL description="Container for NET-seq pipeline"

# seperate environment for umi_tools due to prerequeisite conflicts
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        bedtools htslib=1.13 openssl=1.1 python=3.9 \
        samtools=1.13 star=2.7 wget && \
    micromamba create -y -n umi_tools -c bioconda -c conda-forge umi_tools && \
    micromamba clean --all --yes

WORKDIR /home/micromamba
COPY --chown=micromamba:micromamba scripts/ scripts/

CMD ["/bin/bash"]
