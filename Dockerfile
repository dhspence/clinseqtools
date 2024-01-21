FROM ghcr.io/dhslab/docker-python3:240110

LABEL David Spencer <dspencer@wustl.edu>

LABEL Base Spencerlab image with python3 and some standard bioinformatics modules

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install git+https://github.com/dhspence/clinseqtools.git@dev
