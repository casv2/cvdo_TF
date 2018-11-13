#!/bin/sh
#$ -N tf
#$ -pe mem_smp 16 
#$ -l vmem=512G
#$ -l h_rt=24:00:00
#$ -R y
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -V

export OMP_NUM_THREADS=1

python analyse.py Si PIP_Si_4BBAenv_sw3_reg

