# This file was created by Cassidy Thompson of the BioFrontiers IT department as part of a collaboration between BioFrontiers IT and the Dowell Lab. It is now currently maintained by Margaret Gruca of the Dowell Lab.

MAINTAINER "Margaret Gruca"

FROM fedora:29

RUN dnf -y update && dnf install -y file gcc gcc-c++ git make wget which libtool python3-pip redhat-rpm-config python3-devel zlib-devel libstdc++ openmpi-devel libstdc++-static

RUN cd /tmp && \
    wget http://www.mpich.org/static/downloads/3.3/mpich-3.3.tar.gz && \
    gzip -dc mpich-3.3.tar.gz | tar xf - && \
    cd mpich-3.3 && \
    ./configure --disable-fortran --prefix=/usr/mpich-3.3 && \
    make && \
    make install

ENV PATH            /usr/mpich-3.3/bin:${PATH}
ENV LD_LIBRARY_PATH /usr/mpich-3.3/lib:${LD_LIBRARY_PATH}

RUN cd /usr && git clone https://github.com/Dowell-Lab/FStitch
RUN              cd /usr/FStitch/src && make clean && make
RUN pip3 install FStitch-Bidir --user
ENV PATH            /usr/FStitch/src:${PATH}
ENV PATH            /root/.local/bin:${PATH}

RUN cd /usr && git clone https://github.com/Dowell-Lab/Tfit
RUN              cd /usr/Tfit/src && make clean && make 
ENV PATH            /usr/Tfit/src:${PATH}

CMD /bin/bash

