#!/bin/bash
#SBATCH -J pfs_imQuality 
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 00:10:00
#SBATCH --mem=4GB
#SBATCH -o /net/SRVSTK20C/drp/fmadec/logs/Out_%A_%a.txt
#SBATCH -e /net/SRVSTK20C/drp/fmadec/logs/Err_%A_%a.txt

#date > ${SLURM_ARRAY_TASK_ID}.txt

#cd $PBS_O_WORKDIR # change to submission jobs directory


echoDo ()
{
    echo "$@"
    eval "$@"
}

echo $cam
echo $rerun
echo $doBck

# init environment
if [ -d "/software/bin" ] ; then
    PATH="/software/bin:$PATH"
    fi

export LD_LIBRARY_PATH=""

source scl_source enable devtoolset-8 rh-git218

#export ENGINE_PFS_PATH="///net/SRVSTK20C/data/ait/experimentLog_subaru.db"

#export ENGINE_PFS_PATH="///net/SRVSTK20C/data/ait/experimentLog_subaru_fromMay.db"


source /software/drp/loadLSST.bash
setup pfs_pipe2d

export PYTHONPATH=/drp/devel/lib/:$PYTHONPATH

EUPS_PATH=$EUPS_PATH:/software/mhs/products
setup -v spt_operational_database

cd /drp/devel/ait-notebook/jan-2021/clustercluster/

echoDo python Cluster_GetIMqual2csv.py -c $cam -r $rerun -b $doBck -v ${SLURM_ARRAY_TASK_ID}

sleep 2
