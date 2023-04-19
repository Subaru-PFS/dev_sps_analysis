#!/bin/bash
#SBATCH -J pfs_imQuality 
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 00:10:00
#SBATCH --mem=4GB
#SBATCH -o /net/SRVSTK20C/drp/cluster/logs/Out_%A_%a.txt
#SBATCH -e /net/SRVSTK20C/drp/cluster/logs/Err_%A_%a.txt

#date > ${SLURM_ARRAY_TASK_ID}.txt

#cd $PBS_O_WORKDIR # change to submission jobs directory


echoDo ()
{
    echo "$@"
    eval "$@"
}

echo $cam
echo $rerun
echo $basePath
echo $roi_size
# init environment

export LD_LIBRARY_PATH=""

source ~/.pfs_profile

cd /drp/devel/fmadec/lam_sps_devel/notebooks/devel/fmadec/SM/getImqual

if [ -v $filtered_waves];
then
	echoDo python Cluster_GetBestFocusPlanefromCSV.py -c $cam -r $rerun -b $basePath --roi_size $roi_size -e ${SLURM_ARRAY_TASK_ID}
else
	echoDo python Cluster_GetBestFocusPlanefromCSV.py -c $cam -r $rerun -b $basePath --roi_size $roi_size --filtered_waves $filtered_waves -e ${SLURM_ARRAY_TASK_ID}
fi

sleep 2
