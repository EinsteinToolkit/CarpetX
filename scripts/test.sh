#!/bin/bash

set -ex

export CARPETXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

cd Cactus

# Somehow /usr/local/lib is not in the search path
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
# OpenMPI does not like to run as root (even in a container)
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

time ./simfactory/bin/sim create-run TestJob01_temp_1 --cores 1 --num-threads 2 --testsuite --select-tests=CarpetX
ONEPROC_DIR="$(./simfactory/bin/sim --machine=actions get-output-dir TestJob01_temp_1)/TEST/sim"

time ./simfactory/bin/sim create-run TestJob01_temp_2 --cores 2 --num-threads 1 --testsuite --select-tests=CarpetX
TWOPROC_DIR="$(./simfactory/bin/sim --machine=actions get-output-dir TestJob01_temp_2)/TEST/sim"

#TODO # parse results, generate plots
#TODO cd $PAGESSPACE
#TODO python3 $SCRIPTSPACE/bin/store.py $WORKSPACE/Cactus/repos/cactusamrex $ONEPROC_DIR $TWOPROC_DIR
#TODO 
#TODO python3 $SCRIPTSPACE/bin/logpage.py $WORKSPACE/Cactus/repos/cactusamrex
#TODO 
#TODO # store HTML results
#TODO git add docs
#TODO git add records
#TODO git add test_nums.csv
#TODO if git commit -m "updated html file"; then
#TODO     git config --local user.email "maintainers@einsteintoolkit.org"
#TODO     git config --local user.name "github runner"
#TODO     git push
#TODO fi
