#!/bin/bash

# start mongod
md=`mktemp -d $SCRATCH/af/work/af_tmp/mongo.XXXXXXXXXXXXXXX`
sleep 10
singularity run ~/somoafpipe.sif bash -c "mongod --dbpath $SCRATCH/af/work/af_tmp/mongo.xX8LyW5Jg4f6iCk --logpath $SCRATCH/af/work/af_tmp/mongo.xX8LyW5Jg4f6iCk/log --fork"
