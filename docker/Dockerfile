FROM ubuntu:latest
MAINTAINER Lea Picard lea.picard@gmail.com

RUN apt-get clean && apt-get update && apt-get install -y build-essential \
	git \
	emboss \
	doxygen \
	pkg-config \
	automake \
	openssl \
	libssl-dev \
	openmpi-bin \
	libopenmpi-dev

## Install all the necessary python dependencies to run DGINN

WORKDIR /opt

## Install most recent version of cMake
RUN git clone https://github.com/Kitware/CMake && cd CMake && ./bootstrap && make && make install

## Install PRANK
RUN git clone https://github.com/ariloytynoja/prank-msa.git && cd prank-msa/src && make
ENV PATH /opt/prank-msa/src:$PATH

## Install PhyML
RUN git clone https://github.com/stephaneguindon/phyml && cd phyml && sh ./autogen.sh && ./configure --enable-phyml && make
ENV PATH /opt/phyml/src:$PATH

## Install Treerecs
RUN git clone https://gitlab.inria.fr/Phylophile/Treerecs && cd Treerecs && cmake -DCMAKE_BUILD_TYPE=MinSizeRel . && make && make install
ENV PATH /opt/Treerecs:$PATH

## Install BIO++
RUN mkdir bpp && cd bpp && mkdir sources && cd sources && for d in bpp-core bpp-seq bpp-phyl bpp-popgen; \
	do git clone https://github.com/BioPP/$d; \
	cd $d; \
	cmake -DCMAKE_INSTALL_PREFIX=~/bpp -DCMAKE_LIBRARY_PATH=~/bpp/lib -DCMAKE_INCLUDE_PATH=~/bpp/include -DBUILD_TESTING=FALSE ./; \
	make; \
	make install; \
	cd ..; \
	done

RUN cd bpp && git clone http://biopp.univ-montp2.fr/git/bppsuite && cd bppsuite && cmake -DCMAKE_INSTALL_PREFIX=~/bpp ./ -DINFO=FALSE && make && make install
ENV PATH /opt/bpp/bppsuite/bppSuite/:$PATH
ENV CPATH /opt/bpp/include:$CPATH
ENV LIBRARY_PATH /opt/bpp/lib:$LIBRARY_PATH
ENV LD_LIBRARY_PATH /opt/bpp:/opt/bpp/include:/opt/bpp/lib:$LD_LIBRARY_PATH

## Install python3, necessary packages and DGINN
RUN apt-get install -y python3-pip

RUN pip3 install --target /usr/local/lib/python3.6/dist-packages numpy \
	scipy \
	pandas \
	requests \
	biopython \
	setuptools \
	pyqt5 \
	six \
	lxml

RUN pip3 install --target /usr/local/lib/python3.6/dist-packages --upgrade ete3
ENV PYTHONPATH /usr/local/lib/python3.6/dist-packages:$PYTHONPATH


## Install HYPHY
#RUN git clone https://github.com/veg/hyphy.git && cd hyphy && cmake . && make MPI && make MP && make install
WORKDIR /opt/hyphy
RUN git init && git remote add origin https://github.com/veg/hyphy/ && git fetch --depth 1 origin e1420457b4d0c9de8f206a78614c44525dc7be68 && git checkout FETCH_HEAD && cmake . && make HYPHYMPI && make HYPHYMP && make install
ENV PATH /opt/hyphy:$PATH

RUN ls /usr/local/lib/python3.6/dist-packages
WORKDIR /opt
RUN git clone https://github.com/leapicard/DGINN && cd DGINN && chmod +x DGINN.py
#RUN cd /opt/DGINN/lib/ && sed 's/mpirun -np 2/mpirun -np 2 --allow-run-as-root/g' AnalysisFunc.py > AnalysisFunc2.py && rm AnalysisFunc.py && mv AnalysisFunc2.py AnalysisFunc.py
ENV PATH /opt/DGINN/:$PATH
ENV PYTHONPATH /opt/DGINN/:$PYTHONPATH

#WORKDIR ~
ENV PATH /usr/local/lib/python3.6/dist-packages:$PATH
#RUN ete3
ENTRYPOINT ["DGINN.py"]
CMD ["-h"]
