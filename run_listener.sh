#!/bin/bash

export ROOT_DIR=$PWD
#conda activate /data/des80.a/data/imcmahon/micromamba/envs/mi38
conda activate mi38
python listener --mode $1
