#!/bin/bash

set -ex

export CARPETXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

# Check out Cactus
wget https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
chmod a+x GetComponents
./GetComponents --no-parallel --shallow "$CARPETXSPACE/scripts/carpetx.th"

cd Cactus

# Create a link to the CarpetX repository
ln -s "$CARPETXSPACE" repos
# Create links for the CarpetX thorns
mkdir -p arrangements/CarpetX
pushd arrangements/CarpetX
ln -s ../../repos/CarpetX/* .
popd
