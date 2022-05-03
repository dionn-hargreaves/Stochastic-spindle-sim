#!/bin/bash --login
#$ -cwd               # Application will run from current working directory
#$ -N julia_stochasticSpindle    # Name given to batch job (optional)
#$ -m bea
#$ -M dionn.hargreaves@postgrad.manchester.ac.uk


module load apps/binapps/julia/1.6.2


julia Run_Main.jl
