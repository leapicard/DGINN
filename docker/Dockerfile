FROM ubuntu:latest as base
MAINTAINER Lea Picard lea.picard@gmail.com

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get clean && apt-get update && apt-get install -y apt-utils && apt-get -y install default-jdk default-jre &&  apt-get -y install build-essential git emboss doxygen  pkg-config  automake  openssl  libssl-dev  openmpi-bin  libopenmpi-dev  libeigen3-dev wget

RUN TZ="Europe/Paris" apt-get -y install tzdata


## Install python3, necessary packages and DGINN

RUN apt-get install -y python3-numpy python3-scipy \
	python3-pandas \
	python3-requests \
	python3-setuptools \
	python3-pyqt5 \
	python3-six \
	python3-lxml python3-biopython python3-ete3




## Install most recent version of cMake

RUN apt-get -y install build-essential cmake

### SNAKEMAKE

RUN apt-get install -y snakemake

### Clean apt cache

RUN apt-get clean


############################################
##ALIGNERS

#################
## Install PRANK
FROM base as prank
WORKDIR /opt

RUN git clone https://github.com/ariloytynoja/prank-msa.git && cd prank-msa/src && make


#################
## Install mafft
FROM base as mafft
RUN apt-get install -y mafft


#################
## Install MACSE
FROM base as macse
WORKDIR /opt

RUN wget https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar


##################################################
### Trees

#################
## Install PhyML
FROM base as phyml
RUN apt-get -y install phyml

#################
## Install Iqtree
FROM base as iqtree
RUN apt-get  -y install iqtree


#################
## Install Treerecs
FROM base as treerecs
WORKDIR /opt

RUN git  clone  --recurse-submodules https://gitlab.inria.fr/Phylophile/Treerecs && cd Treerecs && cmake -DCMAKE_BUILD_TYPE=MinSizeRel . && make && make install

#################################################
### Selection

#################
## Install BIO++
FROM base as biopp
WORKDIR /opt

RUN mkdir bpp && cd  bpp && mkdir sources && cd sources
RUN     git clone https://github.com/BioPP/bpp-core; \
	cd bpp-core; \
        mkdir build; \
        cd build; \
        cmake -B . -S .. -DCMAKE_INSTALL_PREFIX=/opt/bpp/;\
        make -j 4; \
        make install;\ 
        cd ../..

RUN     git clone https://github.com/BioPP/bpp-seq; \
	cd bpp-seq; \
        mkdir build; \
        cd build; \
        cmake -B . -S .. -DCMAKE_INSTALL_PREFIX=/opt/bpp/;\
        make -j 4; \
        make install;\ 
        cd ../..
        
RUN     git clone https://github.com/BioPP/bpp-phyl; \
	cd bpp-phyl; \
        mkdir build; \
        cd build; \
        cmake -B . -S .. -DCMAKE_INSTALL_PREFIX=/opt/bpp/;\
        make -j 4; \
        make install;\ 
        cd ../..
        
RUN     git clone https://github.com/BioPP/bpp-popgen; \
	cd bpp-popgen; \
        mkdir build; \
        cd build; \
        cmake -B . -S .. -DCMAKE_INSTALL_PREFIX=/opt/bpp/;\
        make -j 4; \
        make install;\ 
        cd ../..

RUN cd bpp  && git clone https://github.com/BioPP/bppsuite && cd bppsuite && cmake -DCMAKE_INSTALL_PREFIX=/opt/bpp ./ && make -j 4 && make install



#################
## Install paml
FROM base as paml

RUN apt-get install -y paml


#################
### Install HYPHY
FROM base as hyphy

RUN apt-get install -y hyphy-mpi hyphy-common hyphy-pt



############################################################
############################################################
## Install all

#########################"
## DGINN

FROM base
WORKDIR /opt/
RUN git clone https://github.com/lgueguen/DGINN && cd DGINN


### Environnment

ENV PATH /opt/DGINN/:$PATH
ENV PYTHONPATH /opt/DGINN/:/opt/DGINN/scripts/:/opt/DGINN/lib/:$PYTHONPATH


## All softwares
ENV PATH /opt/bin:/usr/bin/$PATH

#################################
## ALIGNERS

COPY --from=prank /opt/prank-msa/src/prank /opt/bin/prank

COPY --from=mafft /usr/bin/mafft /opt/bin/mafft
COPY --from=mafft /usr/bin/mafft-homologs /opt/bin/mafft-homologs
COPY --from=mafft /usr/bin/mafft-profile /opt/bin/mafft-profil

COPY --from=mafft /usr/lib/mafft /usr/lib/mafft   
 
# ENV PATH /opt/mafft/bin/:/opt/mafft/lib/:$PATH
#

COPY --from=macse /opt/macse_v2.07.jar /opt/bin/macse.jar

###############################
## Trees

COPY --from=iqtree /usr/bin/iqtree2 /usr/bin/iqtree2

COPY --from=phyml /usr/bin/phyml /usr/bin/phyml
COPY --from=phyml /usr/lib/phyml/bin /usr/lib/phyml/bin

COPY --from=treerecs  /opt/Treerecs /opt/treerecs
ENV PATH /opt/treerecs/bin:$PATH


######################################
### Selection

COPY --from=biopp  /opt/bpp/lib /opt/bpp/lib
COPY --from=biopp  /opt/bpp/bppsuite/bppSuite /opt/bpp/bppsuite/bppSuite

COPY --from=paml /usr/bin/codeml  /usr/bin/codeml
COPY --from=paml /etc/paml  /etc/paml
COPY --from=paml /usr/lib/paml/bin /usr/lib/paml/bin

COPY --from=hyphy /usr/lib/hyphy/bin /usr/lib/hyphy/bin
COPY --from=hyphy /usr/bin/hyphy /usr/bin/hyphy
COPY --from=hyphy /usr/share/hyphy /usr/share/hyphy

ENV PATH /opt/bpp/bppsuite/bppSuite/:$PATH
ENV LD_LIBRARY_PATH /opt/bpp/lib:$LD_LIBRARY_PATH



### Set Up rights for snakemake

RUN mkdir /opt/DGINN/.snakemake
RUN chmod 777 /opt/DGINN/.snakemake
ENV XDG_CACHE_HOME /opt/DGINN/.snakemake

###

WORKDIR /opt/home

ENTRYPOINT ["/usr/bin/snakemake","-c1","--snakefile","/opt/DGINN/Snakefile"]


CMD ["--help"]
