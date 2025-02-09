#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import time
import numpy as np
import sys
import argparse
import os

# remove FutureWarning: Gen2 Butler has been deprecated (PfsButler). It will be removed sometime after v23.0 but no earlier than the end of 2021.
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from pfs.lam.opdb import *
from pfs.lam.fileHandling import *
from pfs.lam.detFocusAnalysis import *
from pfs.lam.detAnalysis import *
from pfs.lam.linePeaksList import filterPeakList
from pfs.lam.analysisPlot import plotRoiPeak, plotPeaksBrightness

import lsst.daf.persistence as dafPersist

import matplotlib
matplotlib.use('Agg')



def main(visit, peaklist, cam, rerun, experimentId, outpath, drpPath, repo, roi_size, seek_size, doBck, roiPlot, plotPeaksFlux, doFit, doLSF, doPrint, detMap=None, fiberType="DCB"):
    

    
    com = True  # Center Of Mass
    head = 0
    tail = 0
    verbose = False
    doPrint = doPrint
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
            experimentId = get_Visit_Set_Id_fromWeb(visit, url="http://133.40.164.16/sps-logs/index.html")
            print(f"{experimentId} found for visit {visit}")
        except:
            experimentId = np.nan
            raise(f"Unable to get experimentId from logbook for visit {visit}")

    csvPath = os.path.join(outpath,f"sm{specId}",f"Exp{experimentId}",rerun,f"roi{roi_size}",f"doBck{doBck}")
    if doPrint:
        print(csvPath)
        print(f'Start visit {visit} of {cam} with {peaklist}\n')

    # start the butler 

    repoRoot = f"{drpPath}/{repo}"
    if doPrint:
        print(f"drp base: {repoRoot}")
        print(f"rerun {os.path.join(repoRoot,  'rerun', rerun)}")
#    print(os.path.join(repoRoot, "CALIB"))

    butler = dafPersist.Butler( os.path.join(repoRoot, "rerun", rerun)) #, calibRoot=os.path.join(repoRoot, "CALIB"))
    #butler.getKeys('raw')

    # define dataID

    dataId = dict(arm=arm, spectrograph=specId)

    dataId.update(visit=int(visit))


    # get lamp used to filter the list of peaks
    lamps = butler.queryMetadata('raw', ['lamps'], dataId) 
    # IIS lamp name has "_eng"
    if fiberType == "ENGINEERING":
        lamps = [l.replace('_eng', '') for l in lamps]
    
    #calExp = butler.get("calexp", visit=visit, arm=cam[0])
    calExp = butler.get("calexp", dataId)
    [exptime] = butler.queryMetadata('raw', ['exptime'], dataId)
    # DCB line Wheel hole size
    if fiberType == 'DCB':
        lwh = calExp.getMetadata().toDict()['W_AITLWH']
    else:
        lwh = None

    if peaklist is None:
        peaks = None
    else:
        peaks = filterPeakList(peaklist, arm, lamps) if peaklist is not None else None

        waves = peaks.wavelength.unique()
        if doPrint:
            print(waves)
    # !!!!! Wave filter !!!!!!!
    # filter wave if needed
    #print("!!!!!!!!!!!!!!!!!!!!!!!")
    #print(f"wave {waves[-1]} filtered")
    #peaks = peaks[peaks.wavelength != waves[-1]]   

    if not os.path.exists(csvPath):
        os.makedirs(csvPath)
    if roiPlot and (peaklist is not None):
        RoiPlotTitle = f"Peaklist roiPlot {cam.upper()} Exp{experimentId} - visit{visit} - roi_size={roi_size}\n"
        plotRoiPeak(calExp.image.array, peaks, roi_size=roi_size, savePlotFile=os.path.join(csvPath,f"{cam}_{visit}_rawPeak"),raw=True,doSave=True, title=RoiPlotTitle )

    df = ImageQualityToCsv(butler, dataId, peaks, csv_path=csvPath, com=com, doBck=doBck, EE=[3,5],seek_size=seek_size,doFit=doFit, doLSF=doLSF,  doSep=True,mask_size=20, threshold= 50, subpix = 5 , maxPeakDist=80,maxPeakFlux=40000, minPeakFlux=2000,doPlot=roiPlot, doPrint=doPrint, experimentId=experimentId, detMap=detMap, fiberType=fiberType)
    if plotPeaksFlux:
        plotPeaksBrightness(df, doSave=True, savePlotFile=os.path.join(csvPath,f"{cam}_{visit}_fluxes_{'_'.join(lamps)}{exptime:.0f}s_lwh{lwh}"), plot_title=f"{cam}_{visit} - {'_'.join(lamps)} exptime {exptime}s lwh {lwh}")


if __name__ == "__main__":
    start = time.time()
    
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
    parser.add_argument("--seek_size", type=int, default=None,help="distance in px where to seek a peak")
    parser.add_argument("--doBck", action="store_true" ,default=True,help="local bck substraction. default=True")
    parser.add_argument("--roiPlot", action="store_true" ,default=True,help="save an roi plot. default=True")
    parser.add_argument("--plotPeaksFlux", action="store_true" ,default=True,help="save a peak flux plot. default=True")
    parser.add_argument("--doFit", action="store_true" ,default=False,help="Do 2d gaussian fit. default=False")
    parser.add_argument("--doLSF", action="store_true" ,default=False,help="Calculate LSF. default=False")
    parser.add_argument("--detMap", type=str, help="detectorMap filepath", default=None)
    parser.add_argument("--fiberType", type=str, help="type of fiber to use", default="DCB")

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
    doFit = args.doFit
    doLSF = args.doLSF
    detMap = args.detMap
    fiberType = args.fiberType
    
    main(visit, peaklist, cam, rerun, experimentId, outpath, drpPath, repo, roi_size, seek_size, doBck, roiPlot, plotPeaksFlux, doFit, doLSF, detMap, fiberType)
    
    finish = time.time()
    elapsed = finish - start
    print(f"Time elapsed: {elapsed}")

