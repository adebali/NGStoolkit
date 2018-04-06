FROM biocontainers/biocontainers:latest

RUN conda install samtools=1.3.1
RUN conda install bedtools=2.27.1
RUN conda install -c bioconda cutadapt

# RUN apt-get update -y
# RUN apt-get install -y python zip unzip wget curl bzip2 cpanminus
# RUN mkdir -p /home/biodocker/bin

# Bowtie2
ENV ZIP=bowtie2-2.3.4.1-linux-x86_64.zip
ENV URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.3.4.1/
ENV FOLDER=bowtie2-2.3.4.1-linux-x86_64
ENV DST=/home/biodocker/bin

RUN wget $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER

## Install Cutadapt
# RUN apt-get install -y \
# #     python-dev \
# #     build-essential \
# #     virtualenv \
# #     python-pip
# # USER root
# RUN mkdir -p /home/biodocker/cutadapt
# RUN pip install --user --install-option="--prefix=/home/biodocker/cutadapt" --upgrade cutadapt && \
#     cp /home/biodocker/cutadapt/bin/cutadapt /home/biodocker/bin
# RUN export PYTHONPATH=/home/biodocker/cutadapt/lib/python2.7/site-packages/
# USER biodocker

# Install bedGraphToBigWig
# ENV BEDGRAPHTOBIGWIG_VERSION=287
ENV BEDGRAPHTOBIGWIG_URL=http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
ENV DEST_DIR=/home/biodocker/bin/
RUN curl -SLo ${DEST_DIR}/bedGraphToBigWig ${BEDGRAPHTOBIGWIG_URL} && \
    chmod +x ${DEST_DIR}/bedGraphToBigWig

# SRA tools

USER root
WORKDIR /src
ENV VERSION 2.8.2-1
RUN wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-ubuntu64.tar.gz" && \
    tar zxfv sratoolkit.${VERSION}-ubuntu64.tar.gz && \
    cp -r sratoolkit.${VERSION}-ubuntu64/bin/* /home/biodocker/bin
USER biodocker

ADD ./ /toolkit
ENV PATH="/toolkit/bin:${PATH}"
ENV PYTHONPATH="/toolkit/bin/utils:${PYTHONPATH}"

WORKDIR /data
CMD ["/bin/bash"]