# /bin/bash
#############################################################
# sra-tools
# Need bash and sra-toolkit initialization
#############################################################

FROM debian:bullseye-20211220-slim
LABEL description="sra-toolkit from NCIB"

ARG VERSION=current
RUN apt-get update && apt-get --quiet install --yes wget uuid-runtime && apt-get clean
RUN wget -O sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-ubuntu64.tar.gz
RUN tar -xzf sratoolkit.tar.gz -C /usr/local/ --strip-components 1
RUN rm sratoolkit.tar.gz
RUN vdb-config --restore-defaults && \
    echo '/LIBS/GUID = "$(uuidgen)"' >> /root/.ncbi/user-settings.mkfg && \
    echo '/libs/cloud/report_instance_identity = "true"' >> /root/.ncbi/user-settings.mkfg

CMD ["/bin/bash"]
