#!/bin/bash
#SBATCH -J YT0_C
#SBATCH -o YT0_C.out
#SBATCH -e YT0_C.err
#SBATCH -t 36:00:00
#SBATCH --account=project_2010143
#SBATCH -N 1
#SBATCH --mem=68524
#SBATCH --partition=small
#SBATCH --cpus-per-task=36
export GAUSS_SCRDIR=$PWD/$SLURM_JOB_ID.GAUSS_SCRDIR
mkdir -p $GAUSS_SCRDIR

module purge
export g16root=/appl/soft/chem/gaussian/G16RevC.02
source $g16root/g16/bsd/g16.profile
export OMP_NUM_THREADS=1
srun g16 < YT0_C.com >& YT0_C.log

trap "rm -rf $GAUSS_SCRDIR" EXIT
