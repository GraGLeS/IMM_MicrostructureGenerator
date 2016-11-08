#!/usr/bin/env zsh
#BSUB -J MICROGEN
#BSUB -o MICROGEN.%J
#BSUB -e MICROGEN.e%J
#BSUB -M 200000
#BSUB -W 2:00
#BSUB -u miessen@imm.rwth-aachen.de
#BSUB -N
#BSUB -n 32
#BSUB -a "openmp" 
#BSUB -P jara0076
no_numa_balancing MicrostructureGenerator parameters.xml

