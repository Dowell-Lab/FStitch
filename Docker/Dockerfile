FROM fedora:29

# This file was created by Cassidy Thompson of the BioFrontiers IT department as part of a collaboration between BioFrontiers IT and the Dowell Lab. It is now currently maintained by Margaret Gruca of the Dowell Lab.

MAINTAINER "Margaret Gruca"

RUN dnf -y update && dnf install -y file gcc gcc-c++ git make wget which

RUN cd /tmp && \
    wget http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz && \
    gzip -dc mpich-3.2.tar.gz | tar xf - && \
    cd mpich-3.2 && \
    ./configure --disable-fortran --prefix=/opt/mpich-3.2 && \
    make && \
    make install

ENV PATH            /opt/mpich-3.2/bin:${PATH}
ENV LD_LIBRARY_PATH /opt/mpich-3.2/lib:${LD_LIBRARY_PATH}

RUN cd /opt && git clone https://github.com/Dowell-Lab/FStitch
RUN              cd /opt/FStitch/src && make clean && make
RUN pip3 install FStitch-Bidir --user
ENV PATH            /opt/FStitch/src:${PATH}

RUN cd /opt && git clone https://github.com/Dowell-Lab/Tfit
RUN              cd /opt/Tfit/src && make clean && make
ENV PATH            /opt/Tfit/src:${PATH}

CMD /bin/bash

