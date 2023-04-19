#!/bin/bash
#SBATCH -J pfs_imQuality 
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 00:10:00
#SBATCH --mem=4GB
#SBATCH -o /net/SRVSTK20C/drp/cluster/logs/Out_%A_%a.txt
#    #SBATCH -e /net/SRVSTK20C/drp/cluster/logs/Err_%A_%a.txt

#date > ${SLURM_ARRAY_TASK_ID}.txt

#cd $PBS_O_WORKDIR # change to submission jobs directory


echoDo ()
{
    echo "$@"
    eval "$@"
}

echo $visits
echo $cam
echo $rerun
echo $peak
echo $roi_size
echo $doBck
echo $experimentId

# init environment

export LD_LIBRARY_PATH=""

source ~/.pfs_profile

cd /drp/devel/fmadec/lam_sps_devel/notebooks/devel/fmadec/nir
#cd /net/SRVSTK20C/SHARED/devel/fmadec/lam_sps_analysis/notebooks/optical
if [ -v $seek_size];
then
	echoDo python Cluster_MClipGetImqual2csv.py -c $cam -r $rerun -p $peak --roi_size $roi_size --visits $visits --experimentId $experimentId
else
	echoDo python Cluster_MClipGetImqual2csv.py -c $cam -r $rerun -p $peak --roi_size $roi_size --seek_size $seek_size --visits $visits --experimentId $experimentId
fi

sleep 2
