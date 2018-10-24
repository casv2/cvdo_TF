#!/bin/bash
#$ -pe smp 16
#$ -l h_rt=10:00:00
#$ -q 'any|tomsk|orinoco'
#$ -S /bin/bash
#$ -N LocalTest
#$ -j yes
#$ -cwd

#source ~/venv/bin/activate
export OMP_NUM_THREADS=1

python analyse.py Si PIP_Si_4BBAenv_sw3
