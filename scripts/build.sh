#!/bin/bash

set -ex

export CARPETXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

cd Cactus

# Set up SimFactory
cp "$CARPETXSPACE/scripts/actions.cfg" ./simfactory/mdb/optionlists
cp "$CARPETXSPACE/scripts/actions.ini" ./simfactory/mdb/machines
cp "$CARPETXSPACE/scripts/actions.run" ./simfactory/mdb/runscripts
cp "$CARPETXSPACE/scripts/actions.sub" ./simfactory/mdb/submitscripts
cp "$CARPETXSPACE/scripts/actions.th" .
cp "$CARPETXSPACE/scripts/defs.local.ini" ./simfactory/etc

# # For Formaline
# git config --global user.email "carpetx@einsteintoolkit.org"
# git config --global user.name "Github Actions"

# Build
time ./simfactory/bin/sim --machine=actions build -j $(nproc) sim
