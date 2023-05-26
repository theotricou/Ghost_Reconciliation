FROM ubuntu:22.04

RUN apt-get update && apt-get upgrade -y
ARG DEBIAN_FRONTEND=noninteractive


# pyton and ALE dependencies
RUN apt install cmake gcc g++ libboost-dev libboost-serialization-dev\
  libboost-mpi-dev libgsl-dev wget git dirmngr gnupg apt-transport-https\
  ca-certificates software-properties-common python3 python3-dev\
  python3-pip python3-numpy locales -y

RUN pip install ete3 scipy


# Zombi Horde
RUN mkdir /usr/tools

WORKDIR /usr/tools
RUN git clone https://github.com/theotricou/Zombi_Horde.git

ENV HOME /root
# ALE
RUN mkdir bpp
WORKDIR bpp

RUN git clone https://github.com/BioPP/bpp-core
RUN git clone https://github.com/BioPP/bpp-seq
RUN git clone https://github.com/BioPP/bpp-phyl

WORKDIR bpp-core
RUN git checkout tags/v2.4.1 -b version2.4.1
RUN sed -i "45s/^.*/#include <limits>/" src/Bpp/Graph/GlobalGraph.cpp
WORKDIR ../bpp-seq
RUN git checkout tags/v2.4.1 -b version2.4.1
WORKDIR ../bpp-phyl
RUN git checkout tags/v2.4.1 -b version2.4.1
WORKDIR ..

RUN mkdir bpp-core-build
RUN mkdir bpp-seq-build
RUN mkdir bpp-phyl-build

WORKDIR bpp-core-build/
RUN cmake ../bpp-core -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE
RUN make
RUN make install

WORKDIR ../bpp-seq-build/
RUN cmake ../bpp-seq -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE
RUN make
RUN make install

WORKDIR ../bpp-phyl-build/
RUN cmake ../bpp-phyl -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE
RUN make
RUN make install


WORKDIR ../..

RUN git clone https://github.com/ssolo/ALE.git
WORKDIR ALE
RUN mkdir build
WORKDIR build

RUN cmake ..
RUN make

ENV PATH="${PATH}:/usr/tools/ALE/build/bin"


# COPY . .
RUN hash -r

WORKDIR $HOME

# Set the locale
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

RUN apt install -y --no-install-recommends software-properties-common dirmngr
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt install -y r-base r-base-core r-recommended r-base-dev
RUN apt-get -y clean autoclean autoremove

RUN R --slave -e 'install.packages("phylter", repos="https://cloud.r-project.org/", dependencies = TRUE)'
RUN R --slave -e 'install.packages("ape", repos="https://cloud.r-project.org/", dependencies = TRUE)'
RUN R --slave -e 'install.packages("phangorn", repos="https://cloud.r-project.org/", dependencies = TRUE)'
RUN hash -r

WORKDIR $HOME
RUN chmod 755 /usr/tools/Zombi_Horde/Zombi.py

ENV PATH="$PATH:/usr/tools/ALE/build/bin"
ENV PATH="$PATH:/usr/tools/Zombi_Horde/"
RUN echo "export PATH=$PATH:/usr/tools/ALE/build/bin" >> ~/.bashrc
RUN echo "export PATH=$PATH:/usr/tools/Zombi_Horde/" >> ~/.bashrc

RUN echo $PATH

# sudo docker build -f ~/GitHub/Ghost_Reconciliation/scripts/Dockerfile --network=host -t treecoutheo/ghost_reconciliation:0.0.1-1 ./
# docker run -v $PWD:$PWD -v /usr/Zombi treecoutheo/ghost_reconciliation:0.0.1-1 ls -l /usr/Zombi


# docker run treecoutheo/ghost_reconciliation:0.0.1-1 python3 -c "import numpy; print(numpy.__version__)"










FROM ubuntu:22.04
RUN apt-get update && \
    apt-get clean && \
    apt-get install -qy cmake libboost-all-dev g++-11 git make python3 wget \
      ca-certificates openssl build-essential libeigen3-dev

WORKDIR /opt
# Install libboost
RUN mkdir -p {bpp/bpp-core-build,bpp/bpp-seq-build,bpp/bpp-phyl-build}
# Install Bio++ v2.4.1
WORKDIR /opt/bpp
# pull repositories
RUN git clone https://github.com/BioPP/bpp-core
RUN git clone https://github.com/BioPP/bpp-seq
RUN git clone https://github.com/BioPP/bpp-phyl
# freeze to version 2.4.1
WORKDIR /opt/bpp/bpp-core
RUN git checkout tags/v2.4.1 -b version2.4.1
RUN sed -i '45i#include <limits>' src/Bpp/Graph/GlobalGraph.cpp
WORKDIR /opt/bpp/bpp-seq
RUN git checkout tags/v2.4.1 -b version2.4.1
WORKDIR /opt/bpp/bpp-phyl
RUN git checkout tags/v2.4.1 -b version2.4.1
WORKDIR /opt/bpp/bpp-core-build
RUN cmake ../bpp-core -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_TESTING=FALSE
RUN make
RUN make install
WORKDIR /opt/bpp/bpp-seq-build
RUN cmake ../bpp-seq -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_TESTING=FALSE
RUN make
RUN make install
WORKDIR /opt/bpp/bpp-phyl-build
RUN cmake ../bpp-phyl -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_TESTING=FALSE
RUN make
RUN make install
# Install ALE
WORKDIR /opt
RUN mkdir ALE-build
RUN ln -s /usr/include/eigen3/Eigen /usr/include/Eigen
RUN git clone https://github.com/ssolo/ALE.git
WORKDIR /opt/ALE
RUN git checkout bppv241
WORKDIR /opt/ALE-build
RUN cmake ../ALE -DCMAKE_LIBRARY_PATH=/usr/local -DCMAKE_INCLUDE_PATH=/usr/local
RUN make
RUN for binary in $PWD/bin/*; do ln -s $binary /usr/local/bin/; done
ENV LD_LIBRARY_PATH /usr/local/lib


RUN apt-get update && \
    apt-get clean && \
    apt-get install -qy python3 python3-dev

RUN pip install ete3 scipy
