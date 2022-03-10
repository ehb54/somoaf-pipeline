#!/bin/bash

# start mongod
md=/tmp/mongo
mkdir $md
singularity run ~/somoafpipe.sif bash -c "mongod --dbpath $md --logpath $md/log --fork"
