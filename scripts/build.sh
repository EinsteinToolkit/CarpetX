#!/bin/bash

set -ex

export CARPETXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

cd Cactus

# Set up SimFactory
cp "$CARPETXSPACE/scripts/actions-$REAL_PRECISION.cfg" ./simfactory/mdb/optionlists
cp "$CARPETXSPACE/scripts/actions-$REAL_PRECISION.ini" ./simfactory/mdb/machines
cp "$CARPETXSPACE/scripts/actions-$REAL_PRECISION.run" ./simfactory/mdb/runscripts
cp "$CARPETXSPACE/scripts/actions-$REAL_PRECISION.sub" ./simfactory/mdb/submitscripts
cp "$CARPETXSPACE/scripts/carpetx.th" .
cp "$CARPETXSPACE/scripts/defs.local.ini" ./simfactory/etc

# For Formaline
git config --global user.email "carpetx@einsteintoolkit.org"
git config --global user.name "Github Actions"

# Build
# The build log needs to be stored for later.
time ./simfactory/bin/sim --machine="actions-$REAL_PRECISION" build -j $(nproc) sim 2>&1 | tee "build-$REAL_PRECISION".log
