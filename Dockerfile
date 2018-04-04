FROM ubuntu:16.04

RUN apt-get update -y
RUN apt-get install -y python zip unzip wget curl bzip2


## Install Bowtie2
ENV ZIP=bowtie2-2.2.9-linux-x86_64.zip
ENV URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/
ENV FOLDER=bowtie2-2.2.9
ENV DST=/home/biodocker/bin

RUN wget $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER


## Install Cutadapt
RUN apt-get install -y \
    python-dev \
    build-essential \
    virtualenv \
    python-pip

RUN mkdir /opt/cutadapt

RUN pip install --install-option="--prefix=/opt/cutadapt" --upgrade cutadapt && \
    cp /opt/cutadapt/bin/cutadapt /usr/bin/cutadapt

RUN export PYTHONPATH=/opt/cutadapt/lib/python2.7/site-packages/

## Install bedtools

ARG PACKAGE_VERSION=2.27.1
ARG BUILD_PACKAGES="git openssl python build-essential zlib1g-dev"
ARG DEBIAN_FRONTEND=noninteractive

# Update the repository sources list
RUN apt-get update && \
    apt-get install --yes \
              $BUILD_PACKAGES && \
    cd /tmp && \
    git clone https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && \
    git checkout v$PACKAGE_VERSION && \
    make && \
    mv bin/* /usr/local/bin && \
    cd / && \
    rm -rf /tmp/* && \
    apt remove --purge --yes \
              $BUILD_PACKAGES && \
    apt autoremove --purge --yes && \
    apt clean && \
    rm -rf /var/lib/apt/lists/*


# Install bedGraphToBigWig
ENV BEDGRAPHTOBIGWIG_VERSION=287
ENV BEDGRAPHTOBIGWIG_URL=http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v${BEDGRAPHTOBIGWIG_VERSION}/bedGraphToBigWig
ENV DEST_DIR=/usr/local/bin/
RUN curl -SLo ${DEST_DIR}/bedGraphToBigWig ${BEDGRAPHTOBIGWIG_URL} && \
    chmod +x ${DEST_DIR}/bedGraphToBigWig


# ENV PATH $PATH:/opt/miniconda2/bin
CMD ["/bin/bash"]