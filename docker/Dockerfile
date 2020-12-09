FROM ubuntu:latest as base
MAINTAINER Lea Picard lea.picard@gmail.com

RUN apt-get clean && apt-get update

RUN DEBIAN_FRONTEND="noninteractive" TZ="Europe/Paris" apt-get -y install tzdata

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
        
## Install all the necessary python dependencies to run DGINN

WORKDIR /opt

## Install most recent version of cMake
## RUN git clone https://github.com/Kitware/CMake && cd CMake && ./bootstrap && make && make install

RUN apt-get -y install build-essential cmake

## Install PRANK
FROM base as prank
RUN git clone https://github.com/ariloytynoja/prank-msa.git && cd prank-msa/src && make

## Install PhyML
FROM base as phyml
RUN git clone  https://github.com/stephaneguindon/phyml && cd phyml && sh ./autogen.sh && ./configure --enable-phyml && make

## Install Treerecs
FROM base as treerecs
RUN git clone https://gitlab.inria.fr/Phylophile/Treerecs && cd Treerecs && cmake -DCMAKE_BUILD_TYPE=MinSizeRel . && make && make install

## Install BIO++
FROM base as biopp
RUN mkdir bpp && cd bpp && mkdir sources && cd sources && for d in bpp-core; \
	do git clone -b newlik https://github.com/BioPP/$d; \
	cd $d; \
	cmake -DCMAKE_INSTALL_PREFIX=/opt/bpp -DCMAKE_LIBRARY_PATH=/opt/bpp/lib -DCMAKE_INCLUDE_PATH=/opt/bpp/include -DBUILD_TESTING=FALSE ./; \
	make -j 4; \
	make install; \
	cd ..; \
	done

RUN  cd bpp && cd sources && for d in bpp-seq; \
	do git clone -b newlik https://github.com/BioPP/$d; \
	cd $d; \
	cmake  -DCMAKE_INSTALL_PREFIX=/opt/bpp -DCMAKE_LIBRARY_PATH=/opt/bpp/lib -DCMAKE_INCLUDE_PATH=/opt/bpp/include -DBUILD_TESTING=FALSE ./; \
	make -j 4; \
	make install; \
	cd ..; \
	done

RUN cd bpp &&  cd sources && for d in bpp-phyl; \
	do git clone -b newlik https://github.com/BioPP/$d; \
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

RUN cd bpp && git clone -b newlik  https://github.com/BioPP/bppsuite && cd bppsuite && cmake -DCMAKE_INSTALL_PREFIX=/opt/bpp ./ && make -j 4 && make install

## Install HYPHY
FROM base as hyphy
#RUN git clone https://github.com/veg/hyphy.git && cd hyphy && cmake . && make MPI && make MP && make install
WORKDIR /opt/hyphy
RUN git init && git remote add origin https://github.com/veg/hyphy/ && git fetch --depth 1 origin e1420457b4d0c9de8f206a78614c44525dc7be68 && git checkout FETCH_HEAD && cmake . && make HYPHYMPI && make HYPHYMP && make install

## Install paml
FROM base as paml
WORKDIR /opt
RUN git clone https://github.com/binlu1981/paml && cd paml && cd src && make

## Install mafft
FROM base as mafft
RUN apt-get clean  && apt-get update && apt-get install -y mafft

FROM base

## Install python3, necessary packages and DGINN
RUN apt-get install -y python3-pip

RUN pip3 install --upgrade numpy \
	scipy \
	pandas \
	requests \
	biopython \
	setuptools \
	pyqt5==5.14 \
	six \
	lxml

RUN pip3 install --upgrade ete3

RUN git clone https://github.com/lgueguen/DGINN && cd DGINN && chmod +x DGINN.py
ENV PATH /opt/DGINN/:$PATH
ENV PYTHONPATH /opt/DGINN/:$PYTHONPATH


## All softwares
COPY --from=prank /opt/prank-msa/src /opt/prank-msa/src
ENV PATH /opt/prank-msa/src:$PATH

COPY --from=phyml /opt/phyml/src /opt/phyml/src
ENV PATH /opt/phyml/src:$PATH

COPY --from=treerecs  /opt/Treerecs /opt/Treerecs
ENV PATH /opt/Treerecs:$PATH


COPY --from=biopp  /opt/bpp/lib /opt/bpp/lib
COPY --from=biopp  /opt/bpp/bppsuite/bppSuite /opt/bpp/bppsuite/bppSuite

ENV PATH /opt/bpp/bppsuite/bppSuite/:$PATH
ENV LD_LIBRARY_PATH /opt/bpp/lib:$LD_LIBRARY_PATH

COPY --from=hyphy /opt/hyphy /opt/hyphy
ENV PATH /opt/hyphy/:$PATH

COPY --from=paml /opt/paml/src/ /opt/paml/src/
ENV PATH /opt/paml/src/:$PATH

COPY --from=mafft /usr/bin/mafft /usr/bin/mafft



#WORKDIR ~/
ENTRYPOINT ["DGINN.py"]
CMD ["-h"]
