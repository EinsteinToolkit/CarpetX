#!/bin/bash

set -ex

export CARPETXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

cd Cactus

# Set up SimFactory
cp "$SCRIPTSPACE/scripts/actions.cfg" ./simfactory/mdb/optionlists
cp "$SCRIPTSPACE/scripts/actions.ini" ./simfactory/mdb/machines
cp "$SCRIPTSPACE/scripts/actions.run" ./simfactory/mdb/runscripts
cp "$SCRIPTSPACE/scripts/actions.sub" ./simfactory/mdb/submitscripts
cp "$SCRIPTSPACE/scripts/actions.th" .
cp "$SCRIPTSPACE/scripts/defs.local.ini" ./simfactory/etc

# # For Formaline
# git config --global user.email "carpetx@einsteintoolkit.org"
# git config --global user.name "Github Actions"

# Build
time ./simfactory/bin/sim --machine=actions build -j $(nproc) sim
