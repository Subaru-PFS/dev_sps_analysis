#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import time
import numpy as np
import sys
import argparse
import os


from pfs.lam.opdb import *
from pfs.lam.fileHandling import *
from pfs.lam.detFocusAnalysis import *
from pfs.lam.detAnalysis import *
from pfs.lam.linePeaksList import filterPeakList
from pfs.lam.analysisPlot import plotRoiPeak, plotPeaksBrightness

import lsst.daf.persistence as dafPersist

import matplotlib
matplotlib.use('Agg')



def main():
    
    parser = argparse.ArgumentParser(description="Cluster_GetImqual2csv.py argument parser")
    
    parser.add_argument("-v","--visit", type=int, help="data visitId", required=True)
    parser.add_argument("-p","--peak", type=str, help="Peaklist file", required=True)
    parser.add_argument("-c","--cam", type=str, help="camera to process like 'r3'", required=True)
    parser.add_argument("-r","--rerun", type=str, help="data rerun. could be 'ginga/detrend' or your rerun folder ", required=True)
    parser.add_argument("--experimentId", type=int, default=None,help="experimentId or visit_set_id")
    parser.add_argument("--outpath", type=str, default="/data/drp/analysis",help="path where to store the results")
    parser.add_argument("--drpPath", type=str, default="/data/drp",help="main drp folder")
    parser.add_argument("--repo", type=str, default="sps",help="drp repository")
    parser.add_argument("--roi_size", type=int, default=24,help="roi_size in px used to calculate the total flux")
    parser.add_argument("--seek_size", type=int, default=60,help="distance in px where to seek a peak")
    parser.add_argument("--doBck", action="store_true" ,default=True,help="local bck substraction. default=True")
    parser.add_argument("--roiPlot", action="store_true" ,default=True,help="save an roi plot. default=True")
    parser.add_argument("--plotPeaksFlux", action="store_true" ,default=True,help="save a peak flux plot. default=True")
    
    args = parser.parse_args()
    
    visit = args.visit
    peaklist = args.peak
    cam = args.cam
    rerun = args.rerun
    experimentId = args.experimentId
    outpath = args.outpath
    drpPath = args.drpPath
    repo = args.repo
    roi_size = args.roi_size
    seek_size = args.seek_size
    doBck = args.doBck
    roiPlot = args.roiPlot
    plotPeaksFlux = args.plotPeaksFlux
    
    
    com = True  # Center Of Mass
    head = 0
    tail = 0
    verbose = False
    doPrint = False
    arm = cam[0]
    specId = int(cam[1])    

  # get experimentId from visit to create CSV folder
    if experimentId is None:
#        try:
    #            experimentId = Logbook.visitExperimentId(visit=visit)
#                experimentId = get_Visit_Set_Id(visit)
#        except:
        try:
            print(f"Try to get expId for visit={visit}")
            experimentId = get_Visit_Set_Id_fromWeb(visit)
            print(f"{experimentId} found for visit {visit}")
        except:
            experimentId = np.nan
            raise(f"Unable to get experimentId from logbook for visit {visit}")

    csvPath = os.path.join(outpath,f"sm{specId}",f"Exp{experimentId}",rerun,f"roi{roi_size}",f"doBck{doBck}")
    print(csvPath)
    print(f'Start visit {visit} of {cam} with {peaklist}')
    print(f"{drpPath}/{repo}/rerun/{rerun}")

    # start the butler 

    repoRoot = f"{drpPath}/{repo}"
    print(repoRoot)
    print(os.path.join(repoRoot,  "rerun", rerun))
    print(os.path.join(repoRoot, "CALIB"))

    butler = dafPersist.Butler( os.path.join(repoRoot, "rerun", rerun), calibRoot=os.path.join(repoRoot, "CALIB"))
    #butler.getKeys('raw')

    # define dataID

    dataId = dict(arm=arm, spectrograph=specId)

    dataId.update(visit=int(visit))


    # get lamp used to filter the list of peaks
    lamps = butler.queryMetadata('raw', ['lamps'], dataId) 
    calExp = butler.get("calexp", visit=visit, arm=cam[0])
    [exptime] = butler.queryMetadata('raw', ['exptime'], dataId)
    # DCB line Wheel hole size
    lwh = calExp.getMetadata().toDict()['W_AITLWH']

    peaks = filterPeakList(peaklist, arm, lamps)

    waves = peaks.wavelength.unique()
    print(waves)
    # !!!!! Wave filter !!!!!!!
    # filter wave if needed
    #print("!!!!!!!!!!!!!!!!!!!!!!!")
    #print(f"wave {waves[-1]} filtered")
    #peaks = peaks[peaks.wavelength != waves[-1]]   

    if not os.path.exists(csvPath):
        os.makedirs(csvPath)
    if roiPlot:
        plotRoiPeak(calExp.image.array, peaks, roi_size=roi_size, savePlotFile=os.path.join(csvPath,f"{cam}_{visit}"),raw=True,doSave=True)

    df = ImageQualityToCsv(butler, dataId, peaks, csv_path=csvPath, com=com, doBck=doBck, EE=[3,5],seek_size=seek_size,doFit=False, doLSF=False,  doSep=True,mask_size=20, threshold= 50, subpix = 5 , maxPeakDist=80,maxPeakFlux=40000, minPeakFlux=2000,doPlot=False, doPrint=doPrint)
    if plotPeaksFlux:
        plotPeaksBrightness(df, doSave=True, savePlotFile=os.path.join(csvPath,f"{cam}_{visit}_fluxes_{'_'.join(lamps)}{exptime:.0f}s_lwh{lwh}"), plot_title=f"{cam}_{visit} - {'_'.join(lamps)} exptime {exptime}s lwh {lwh}")


if __name__ == "__main__":
    start = time.time()
    main()
    finish = time.time()
    elapsed = finish - start
    print(f"Time elapsed: {elapsed}")
