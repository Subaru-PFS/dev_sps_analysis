#!/bin/bash
#SBATCH -J Detrend                                                   # run's name
#SBATCH -N 1                                                         # request 1 node
#SBATCH -c 1                                                         # request 1 cpu per task
#SBATCH --mem=10GB                                                   # request 10GB
#SBATCH -t 00:15:00                                                  # request 12 hours walltime
#SBATCH -o /net/SRVSTK20C/drp/cluster/logs/Out_%A_%a.txt             # output file name
#SBATCH -e /net/SRVSTK20C/drp/cluster/logs/Err_%A_%a.txt             # error file name

#cd $PBS_O_WORKDIR # change to submission jobs directory
#date > ${SLURM_ARRAY_TASK_ID}.txt

#cd $PBS_O_WORKDIR # change to submission jobs directory
umask 2

echoDo ()
{
    echo "$@"
    eval "$@"
}


echo $repo
echo $calib
echo $rerun
echo $dark
echo $bias
echo $flat
echo $defect


export LD_LIBRARY_PATH=""
export OMP_NUM_THREADS=1

source ~/.pfs_profile

detrend.py $repo --calib $calib --rerun $rerun --id visit=${SLURM_ARRAY_TASK_ID}  --no-versions -c isr.doBias=$bias -c isr.doDark=$dark -c isr.doFlat=$flat -c isr.doDefect=$defect -c isr.doSaturationInterpolation=True repair.cosmicray.nCrPixelMax=250000 -c repair.cosmicray.keepCRs=True

sleep 2
