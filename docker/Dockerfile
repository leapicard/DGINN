FROM ubuntu:latest as base
MAINTAINER Lea Picard lea.picard@gmail.com

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get  clean && apt-get update && apt-get install -y apt-utils

RUN TZ="Europe/Paris" apt-get -y install tzdata

RUN apt-get -y install build-essential \
	git \
	emboss \
	doxygen \
	pkg-config \
	automake \
	openssl \
	libssl-dev \
	openmpi-bin \
	libopenmpi-dev \
        libeigen3-dev
        
## Install most recent version of cMake
## RUN git clone https://github.com/Kitware/CMake && cd CMake && ./bootstrap && make && make install

RUN apt-get -y install build-essential cmake


#################
## Install PRANK
FROM base as prank
WORKDIR /opt

RUN git clone https://github.com/ariloytynoja/prank-msa.git && cd prank-msa/src && make


#################
## Install PhyML
FROM base as phyml
WORKDIR /opt

##RUN git clone  https://github.com/stephaneguindon/phyml && cd phyml && sh ./autogen.sh && ./configure --enable-phyml && make
RUN apt-get -y install phyml


#################
## Install Treerecs
FROM base as treerecs
WORKDIR /opt

RUN git  clone  --recurse-submodules https://gitlab.inria.fr/Phylophile/Treerecs && cd Treerecs && cmake -DCMAKE_BUILD_TYPE=MinSizeRel . && make && make install


#################
## Install BIO++
FROM base as biopp
WORKDIR /opt

RUN mkdir bpp && cd  bpp && mkdir sources && cd sources && for d in bpp-core; \
	do git clone https://github.com/BioPP/$d; \
	cd $d; \
	cmake  -DCMAKE_INSTALL_PREFIX=/opt/bpp -DCMAKE_LIBRARY_PATH=/opt/bpp/lib -DCMAKE_INCLUDE_PATH=/opt/bpp/include -DBUILD_TESTING=FALSE ./; \
	make -j 4; \
	make install; \
	cd ..; \
	done

RUN  cd bpp && cd sources && for d in bpp-seq; \
	do git clone https://github.com/BioPP/$d; \
	cd $d; \
	cmake  -DCMAKE_INSTALL_PREFIX=/opt/bpp -DCMAKE_LIBRARY_PATH=/opt/bpp/lib -DCMAKE_INCLUDE_PATH=/opt/bpp/include -DBUILD_TESTING=FALSE ./; \
	make -j 4; \
	make install; \
	cd ..; \
	done

RUN cd bpp && cd sources && for  d in bpp-phyl; \
	do git clone https://github.com/BioPP/$d; \
	cd $d; \
	cmake -DCMAKE_INSTALL_PREFIX=/opt/bpp -DCMAKE_LIBRARY_PATH=/opt/bpp/lib -DCMAKE_INCLUDE_PATH=/opt/bpp/include -DBUILD_TESTING=FALSE ./; \
	make -j 4; \
	make install; \
	cd ..; \
	done

RUN cd bpp && cd sources && for d in bpp-popgen; \
	do git clone https://github.com/BioPP/$d; \
	cd $d; \
	cmake -DCMAKE_INSTALL_PREFIX=/opt/bpp -DCMAKE_LIBRARY_PATH=/opt/bpp/lib -DCMAKE_INCLUDE_PATH=/opt/bpp/include -DBUILD_TESTING=FALSE ./; \
	make; \
	make install; \
	cd ..; \
	done

RUN cd bpp && git clone https://github.com/BioPP/bppsuite && cd bppsuite && cmake -DCMAKE_INSTALL_PREFIX=/opt/bpp ./ && make -j 4 && make install


#################
## Install HYPHY
FROM base as hyphy
WORKDIR /opt

RUN git clone https://github.com/veg/hyphy.git && cd hyphy && cmake . && make MPI && make MP && make install


#################
## Install paml
FROM base as paml
WORKDIR /opt
RUN git clone  https://github.com/abacus-gene/paml && cd paml && cd src && make


#################
## Install mafft
FROM base as mafft
#No WORKDIR /opt  because not used for apt-get
RUN apt-get clean  && apt-get update && apt-get install -y mafft


####################
## Install all
FROM base

WORKDIR /opt/

## Install python3, necessary packages and DGINN
RUN apt-get install -y python3-pip

RUN pip3 install --upgrade numpy \
	scipy \
	pandas \
	requests \
	setuptools \
	pyqt5==5.14 \
	six \
	lxml

RUN pip3 install biopython

RUN pip3 install --upgrade ete3

RUN git clone https://github.com/lgueguen/DGINN && cd DGINN && chmod +x DGINN.py
ENV PATH /opt/DGINN/:$PATH
ENV PYTHONPATH /opt/DGINN/:$PYTHONPATH


## All softwares
ENV PATH /opt/bin:$PATH

COPY --from=prank /opt/prank-msa/src/prank /opt/bin/prank

COPY --from=mafft /usr/bin/mafft /opt/mafft/bin/mafft
COPY --from=mafft /usr/bin/mafft-homologs /opt/mafft/bin/mafft-homologs
COPY --from=mafft /usr/bin/mafft-profile /opt/mafft/bin/mafft-profil
COPY --from=mafft /usr/lib/mafft /usr/lib/mafft   
# 
ENV PATH /opt/mafft/bin/:/opt/mafft/lib/:$PATH
#ENV LD_LIBRARY_PATH /opt/local/lib:$LD_LIBRARY_PATH
# RUN echo $PATH
#


COPY --from=phyml /usr/lib/phyml/bin/phyml /opt/bin/phyml


COPY --from=treerecs  /opt/Treerecs /opt/treerecs
ENV PATH /opt/treerecs/bin:$PATH

COPY --from=biopp  /opt/bpp/lib /opt/bpp/lib
COPY --from=biopp  /opt/bpp/bppsuite/bppSuite /opt/bpp/bppsuite/bppSuite

ENV PATH /opt/bpp/bppsuite/bppSuite/:$PATH
ENV LD_LIBRARY_PATH /opt/bpp/lib:$LD_LIBRARY_PATH

RUN apt-get install -y hyphy-mpi hyphy-common hyphy-pt
RUN echo $PATH

COPY --from=paml /opt/paml/src/codeml /usr/local/bin/

RUN apt-get install -y python3-biopython

#WORKDIR ~/
ENTRYPOINT ["DGINN.py"]
CMD ["-h"]
