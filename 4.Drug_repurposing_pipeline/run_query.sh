#!/bin/bash

# Set the working directory 
#$ -cwd

# Input / output streams
#$ -e out.err
#$ -o out.out

dir=""

python3 query_controller.py \
  --up ${dir}/up.grp \
  --down ${dir}/down.grp \
  --outdir ${dir}/ \
  --njobs 10 \
  --all \
  --included ${dir}/included.grp \
  --login ${dir}/login.cfg

python3 quantile.py \
  --pickle ${dir}/results.pickle \
  --outdir ${dir}/ \
  --min_n_sig 10
