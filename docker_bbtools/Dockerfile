# /bin/bash
#############################################################
# bbtools -- from US DOE JGI
# see https://jgi.doe.gov/data-and-tools/bbtools/
#############################################################

FROM debian:bullseye-20211220-slim
LABEL description="BBTools from JGI DOE"


# TODO Make non-root user
# TODO path or bbtools command

WORKDIR /root


RUN apt-get update && \
    apt-get install --yes wget default-jre

RUN wget -O BBMap.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_38.94.tar.gz
RUN tar -xzf BBMap.tar.gz
RUN rm BBMap.tar.gz


CMD ["/bin/bash"]
