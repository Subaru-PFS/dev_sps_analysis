#!/bin/bash

usage() {
    echo "perform detrend" 1>&2
    echo "" 1>&2
    echo "Usage: $0 -v <VISIT> -t <TARGET> -r <RERUN> -d <doDark> -f <doFlat> -c <keepCR>"
    echo "  -v <VISIT>" 1>&2
    echo "" 1>&2
    exit 1
}

ODARK="isr.doDark=False"
OFLAT="isr.doFlat=False"
OCR="repair.cosmicray.keepCRs=False"

while getopts ":v:t:r:dfc" opt; do
    case "${opt}" in
        v)
            VISIT=${OPTARG}
            ;;
        t)
            TARGET=${OPTARG}
            ;;
        r)
            RERUN=${OPTARG}
            ;;
        d)
            ODARK="isr.doDark=True"
            ;;
        f)
            OFLAT="isr.doFlat=True"
            ;;
        c)
            OCR="repair.cosmicray.keepCRs=True"
            ;;
        *)
            usage
            ;;
    esac
done


RAW_DATA_DIR="/data/pfs/pfs"
#RERUN="pfs"  # Rerun name to use
#TARGET="/drp/lam-test" # Directory name to give data repo

batchArgs2="-c repair.cosmicray.nCrPixelMax=100000 -j 4"
#batchArgs="--clobber-versions  -c isr.doFlat=False -c isr.doDark=False -c repair.cosmicray.keepCRs=False"
batchArgs="--clobber-versions  -c $OFLAT -c $ODARK -c $OCR"



echo " $VISIT"
echo "detrend.py $TARGET --calib $TARGET/CALIB --rerun $RERUN/detrend --id visit=$VISIT $batchArgs $batchArgs2"

detrend.py $TARGET --calib $TARGET/CALIB --rerun $RERUN/detrend --id visit=$VISIT $batchArgs $batchArgs2

