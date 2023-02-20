#!/bin/bash

set -ex

export CARPETXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

cd Cactus
./simfactory/bin/sim setup-silent --optionlist "$CARPETXSPACE/scripts/debian.cfg"

# for Formaline
git config --global user.email "carpetx@einsteintoolkit.org"
git config --global user.name "Github Actions"

NCPUS=$(nproc)
time ./simfactory/bin/sim build -j$NCPUS --thornlist repos/cactusamrex/azure-pipelines/carpetx.th 2>&1 | tee build.log
