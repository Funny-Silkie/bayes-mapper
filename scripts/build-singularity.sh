#!/bin/bash

cd $(dirname $0)
cd ../
singularity build --fakeroot package/bayes-mapper.sif package/singularity.def
