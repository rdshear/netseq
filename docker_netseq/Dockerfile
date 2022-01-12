# /bin/bash -c docker build . -t rdshear/netseq
#############################################################
# Net-seq allignment (fastq to bedgraph)
#############################################################

FROM mambaorg/micromamba:0.18.2
LABEL description="Container for NET-seq pipeline"

# TODO Remove star 
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        bedtools bwa fastp htslib openssl samtools sra-tools star wget && \
    micromamba clean --all --yes

WORKDIR /home/micromamba

RUN vdb-config --restore-defaults && \
    echo '/LIBS/GUID = "$(uuidgen)"' >> $HOME/.ncbi/user-settings.mkfg && \
    echo '/libs/cloud/report_instance_identity = "true"' >> $HOME/.ncbi/user-settings.mkfg


CMD ["/bin/bash"]
