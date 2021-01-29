#!/bin/bash

# change 2020/11/30
# add echoDo
# add pfsDesign
#
# change 2020/11/20
# update RAW_DATA_DIR
#
# change 2021/01/20
# modif to support cluster repo
# add an argument cluster to set folder according


echoDo ()
{
            echo "$@"
            printf "\n"
            eval "$@"
}



usage() {
    echo "perform ingest" 1>&2
    echo "" 1>&2
    echo "Usage: $0 -d <DATE> -t <TARGET> -c"
    echo "  -c : set cluster folders" 1>&2
    echo "" 1>&2
    exit 1
}

CLUSTER='false'
DRY_RUN='false'

while getopts ":d:t:c:n:" opt; do
    case "${opt}" in
        d)
            DATE=${OPTARG}
            ;;
        t)
            TARGET=${OPTARG}
            ;;
        c)
            CLUSTER='true'
            ;;
        n)
            DRY_RUN='true'
            ;;
        *)
            usage
            ;;
    esac
done


RAW_DATA_DIR="/data/raw"
PFS_CONFIG_DIR="/data/drp/pfsDesign"

#RERUN="pfs"  # Rerun name to use
#TARGET="/drp/lam-test" # Directory name to give data repo


if [[ "${CLUSTER}" == true ]]; then
    RAW_DATA_DIR="/net/SRVSTK20C/data/raw"
    PFS_CONFIG_DIR="/net/SRVSTK20C/drp/sps-cluster/pfsDesign"
fi

batchArgs="--batch-type=none --doraise"

echo "$DATE"
#echo " ingestPfsImages.py $TARGET --mode=link \
# $RAW_DATA_DIR/$DATE/sps/*.fits \
# -c clobber=True register.ignore=True"

echoDo ingestPfsImages.py $TARGET --mode=link \
        $RAW_DATA_DIR/$DATE/sps/*.fits \
        --pfsConfigDir $PFS_CONFIG_DIR \
        -c clobber=True register.ignore=True
