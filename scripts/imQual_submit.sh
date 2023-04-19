#!/bin/bash
#SBATCH -J pfs_imQuality 
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 00:10:00
#SBATCH --mem=4GB
#SBATCH -o /net/SRVSTK20C/drp/cluster/logs/Out_%A_%a.txt

#date > ${SLURM_ARRAY_TASK_ID}.txt

#cd $PBS_O_WORKDIR # change to submission jobs directory


echoDo ()
{
    echo "$@"
    eval "$@"
}

echo $cam
echo $rerun
echo $peak
echo $roi_size

# init environment

export LD_LIBRARY_PATH=""

source ~/.pfs_profile

cd /drp/devel/fmadec/lam_sps_devel/notebooks/devel/fmadec/SM/getImqual
#cd /net/SRVSTK20C/SHARED/devel/fmadec/lam_sps_analysis/notebooks/optical
if [ -v $seek_size];
then
	seek_size=40
fi
if [ -v $doFit];
then
        doFit=False
fi
if [ -v $doLSF];
then
        doLSF=False
fi

echoDo python Cluster_GetImqual2csv.py -c $cam -r $rerun -p $peak --roi_size $roi_size --seek_size $seek_size --doFit $doFit --doLSF $doLSF -v ${SLURM_ARRAY_TASK_ID}

sleep 2
