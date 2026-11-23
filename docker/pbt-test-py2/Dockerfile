FROM ubuntu:14.04

MAINTAINER Ryan Dale <dalerr@niddk.nih.gov>

RUN apt-get update && apt-get install -y \
    build-essential \
    bzip2 \
    ca-certificates \
    git \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender1 \
    mysql-client \
    wget \
    zlib1g-dev

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda-3.10.1-Linux-x86_64.sh && \
    echo "363f56f5608d1552325549e7371fcf460c5ed45484eb300058e3b99c997808b5  Miniconda-3.10.1-Linux-x86_64.sh" | sha256sum --check --strict - && \
    /bin/bash /Miniconda-3.10.1-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda-3.10.1-Linux-x86_64.sh && \
    /opt/conda/bin/conda install --yes conda==3.14.1
ENV PATH /opt/conda/bin:$PATH

RUN conda install -c daler \
    pip \
    cython \
    matplotlib \
    nose \
    numpydoc \
    pip \
    pandas \
    pyyaml \
    sphinx \
    pysam
RUN conda install -c daler \
    tabix \
    bedtools=2.25.0
ENV DISPLAY=:0
ENV LANG C.UTF-8
WORKDIR /opt/pybedtools
