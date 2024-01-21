FROM ghcr.io/dhslab/docker-python3:240110

LABEL David Spencer <dspencer@wustl.edu>

LABEL Base Spencerlab image with python3 and some standard bioinformatics modules

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /tmp/

RUN git clone -b dev https://github.com/dhspence/clinseqtools.git && \
    cd clinseqtools && \
    pip install . && \
    cd .. && \
    rm -rf clinseqtools \
