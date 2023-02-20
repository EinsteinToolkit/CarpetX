#!/bin/bash

set -ex

pwd
echo $PWD
ls -l
ls -l $PWD/CarpetX

export CARPETXSPACE="$PWD/CarpetX"
export WORKSPACE="$PWD/workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

# Check out Cactus
wget https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
chmod a+x GetComponents
./GetComponents --no-parallel --shallow https://bitbucket.org/einsteintoolkit/manifest/raw/master/einsteintoolkit.th

# Unshallow Carpetx
pushd "$CARPETXSPACE"
git fetch --unshallow
echo "At git ref: $(git rev-parse HEAD)"
popd

cd Cactus

# Download old CarpetX thorns for the external libraries there
# TODO: Move these into ExternalLibraries
pushd repos
git clone https://bitbucket.org/eschnett/cactusamrex
popd
pushd arrangements/ExternalLibraries
ln -s ../../repos/cactusamrex/ADIOS2 .
ln -s ../../repos/cactusamrex/AMReX .
ln -s ../../repos/cactusamrex/Boost .
ln -s ../../repos/cactusamrex/NSIMD .
ln -s ../../repos/cactusamrex/Silo .
ln -s ../../repos/cactusamrex/openPMD_api .
ln -s ../../repos/cactusamrex/yaml_cpp .
popd

# Create a link to the CarpetX repository
ln -s "$CARPETXSPACE" repos
# Create links for the CarpetX thorns
mkdir -p arrangements/CarpetX
pushd arrangements/CarpetX
ln -s ../../repos/CarpetX/* .
popd
