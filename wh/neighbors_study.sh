#! /bin/env bash

set -o nounset
set -o errexit

#directories
source jobid.sh

export jobid=$jobid8
rake kNNSplitted
python NeighborsAnalyzer.py

export jobid=$jobid7
rake kNNSplitted
python NeighborsAnalyzer.py

