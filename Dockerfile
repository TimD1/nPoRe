FROM ubuntu:20.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/bin:$PATH DEBIAN_FRONTEND=noninteractive

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

WORKDIR /opt/bin

RUN git clone https://github.com/TimD1/npore && cd npore && make

RUN /bin/bash -c "source /opt/bin/npore/venv3/bin/activate" && \
    echo ". /opt/bin/npore/venv3/bin/activate" > ~/.bashrc

WORKDIR /opt/bin/npore/src
