#!/bin/bash -e

# sudo yum update && \
# sudo yum install -y build-essential \
# libseccomp-dev pkg-config squashfs-tools cryptsetup

# sudo rm -r /usr/local/go

# export VERSION=1.13.15 OS=linux ARCH=amd64  # change this as you need

# wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz && \
# sudo tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz

# echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
# echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
# source ~/.bashrc

# curl -sfL https://install.goreleaser.com/github.com/golangci/golangci-lint.sh |
# sh -s -- -b $(go env GOPATH)/bin v1.21.0

# mkdir -p ${GOPATH}/src/github.com/sylabs && \
# cd ${GOPATH}/src/github.com/sylabs && \
# git clone https://github.com/sylabs/singularity.git && \
# cd singularity

# git checkout v3.7.3

# cd ${GOPATH}/src/github.com/sylabs/singularity && \
# ./mconfig && \
# cd ./builddir && \
# make && \
# sudo make install

sudo yum update -y && \
    sudo yum groupinstall -y 'Development Tools' && \
    sudo yum install -y \
    openssl-devel \
    libuuid-devel \
    libseccomp-devel \
    wget \
    squashfs-tools

export GOVERSION=1.13.15 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$GOVERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$GOVERSION.$OS-$ARCH.tar.gz && \
    rm go$GOVERSION.$OS-$ARCH.tar.gz

echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

#curl -sfL https://install.goreleaser.com/github.com/golangci/golangci-lint.sh |
#sh -s -- -b $(go env GOPATH)/bin v1.21.0

export VERSION=3.7.3 && \
    mkdir -p $GOPATH/src/github.com/sylabs && \
    cd $GOPATH/src/github.com/sylabs && \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
    tar -xzf singularity-${VERSION}.tar.gz && \
    cd ./singularity && \
    ./mconfig 

./mconfig && \
    make -C ./builddir && \
    sudo make install -C ./builddir