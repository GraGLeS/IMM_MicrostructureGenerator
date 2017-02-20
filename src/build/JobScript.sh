#!/usr/bin/env zsh
#BSUB -J MICROGEN
#BSUB -o MICROGEN.%J
#BSUB -e MICROGEN.e%J
#BSUB -M 200000
#BSUB -W 2:00
#BSUB -u miessen@imm.rwth-aachen.de
#BSUB -N
#BSUB -n 128
#BSUB -a "bcs openmp"
#BSUB -P rwth0157
no_numa_balancing MicrostructureGenerator CR_55.xml

