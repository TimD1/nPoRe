FROM ubuntu:20.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/bin:$PATH DEBIAN_FRONTEND=noninteractive

# install packages
RUN apt-get update --fix-missing && \
    yes | apt-get upgrade && \
    apt-get install -y \
        samtools \
        minimap2 \
        git \
        make \
        g++ \
        python3.8-venv \
        python3-setuptools \
        python3-pip

# clone repo
WORKDIR /opt/bin
RUN git clone https://github.com/TimD1/npore
WORKDIR /opt/bin/npore

# setup virtual environment
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install -r requirements.txt
RUN make npore
WORKDIR /opt/bin/npore/src
