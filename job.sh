#!/bin/bash

WORKDIR="/home/de3u14/lib/projects/MWTC/MWTC-ParameterScan"

module load gcc/6.1.0
source /home/de3u14/lib/build/miniconda/build/envs/py27/bin/activate py27

cd ${WORKDIR}
python -V
python mwtc-xsec-scan.py job
