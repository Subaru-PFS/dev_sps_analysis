#!/bin/bash
#SBATCH -J pfs_detrend_array 
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 00:10:00
#SBATCH --mem=2GB
#SBATCH -o /net/SRVSTK20C/drp/cluster/logs/out/Out_%A_%a.txt
#SBATCH -e /net/SRVSTK20C/drp/cluster/logs/err/Err_%A_%a.txt

#date > ${SLURM_ARRAY_TASK_ID}.txt

#cd $PBS_O_WORKDIR # change to submission jobs directory

echoDo ()
{
    echo "$@"
    eval "$@"
}


echo $repo
echo $rerun

# init environment
if [ -d "/software/bin" ] ; then
    PATH="/software/bin:$PATH"
    fi

export LD_LIBRARY_PATH=""
export OMP_NUM_THREADS=1

source scl_source enable devtoolset-8 rh-git218


source /software/drp/loadLSST.bash
setup pfs_pipe2d


DRP_detrend_LAM.sh -t $repo -r $rerun -v ${SLURM_ARRAY_TASK_ID}

sleep 2




