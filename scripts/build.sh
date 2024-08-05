#!/bin/bash

# set -ex

export CARPETXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

cd Cactus

# Set up SimFactory
cp "$CARPETXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.cfg" ./simfactory/mdb/optionlists
cp "$CARPETXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.ini" ./simfactory/mdb/machines
cp "$CARPETXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.run" ./simfactory/mdb/runscripts
cp "$CARPETXSPACE/scripts/actions-$ACCELERATOR-$REAL_PRECISION.sub" ./simfactory/mdb/submitscripts
cp "$CARPETXSPACE/scripts/carpetx.th" .
cp "$CARPETXSPACE/scripts/defs.local.ini" ./simfactory/etc

# For Formaline
git config --global user.email "carpetx@einsteintoolkit.org"
git config --global user.name "Github Actions"

case "$MODE" in
    debug) mode='--debug';;
    optimize) mode='--optimise';;
    *) exit 1;;
esac

# Build
# The build log needs to be stored for later.
time ./simfactory/bin/sim --machine="actions-$ACCELERATOR-$REAL_PRECISION" build "$mode" --jobs $(nproc) sim 2>&1 |
    tee build.log

# TODO
du -sh $WORKSPACE
du -sh $WORKSPACE/configs/sim
du -sh $WORKSPACE/configs/sim/build
du -sh $WORKSPACE/configs/sim/lib
du -sh $WORKSPACE/configs/sim/scratch
du -sh $WORKSPACE/configs/sim/*
du -sh $WORKSPACE/configs/sim/build/*
du -sh $WORKSPACE/configs/sim/lib/*
du -sh $WORKSPACE/configs/sim/scratch/*

# Check whether the executable exists and is executable
test -x exe/cactus_sim

exit $?                         # TODO
