#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import time
import numpy as np
import sys
import getopt

# In[ ]:


import os


# In[ ]:


from pfs.lam.opdb import *
from pfs.lam.fileHandling import *
from pfs.lam.detFocusAnalysis import *
from pfs.lam.detAnalysis import *
from pfs.lam.linePeaksList import filterPeakList
from pfs.lam.analysisPlot import plotRoiPeak, plotPeaksBrightness

# In[ ]:


import lsst.daf.persistence as dafPersist


# In[ ]:

import matplotlib
matplotlib.use('Agg')



# In[ ]:


def main(argv):
    visit = ''
    peak = ''
    outpath = os.path.join('/data/drp/analysis')
    drpPath = "/data/drp"
    repo = "sps"
    repoRoot = f"{drpPath}/{repo}"
    experimentId = None

    try:
        opts, args = getopt.getopt(argv,"hv:p:c:o:r:e:",["visit=","peak=", "cam=", "outpath=", "rerun=","experimentId="])
    except getopt.GetoptError:
        print('argument error')
        print('GetIMqual2csv-standalone.py -v <visit> -p <peakfile> -c <cam> -r <rerun> -o <outpath>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('GetIMqual2csv-standalone.py -v <visit> -p <peakfile> -c <cam>')
            sys.exit()
        elif opt in ("-v", "--visit"):
            visit = int(arg)
        elif opt in ("-p", "--peak"):
            peaklist = arg
        elif opt in ("-c", "--cam"):
            cam = arg
        elif opt in ("-o", "--outpath"):
            outpath = arg
        elif opt in ("-r", "--rerun"):
            rerun = arg
        elif opt in ("-e", "--experimentId"):
            experimentId = arg
           

        
    roiPlot = True
    plotPeaksFlux = True
    #outpath = "output\\" if outpath is None else outpath

    # define defaut parameters
    roi_size = 30
    seek_size = 60

    com = True  # Center Of Mass
    doBck = True
    head = 0
    tail = 0
    criteria = 'EE5'
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

    butler = dafPersist.Butler( os.path.join(repoRoot, "rerun", rerun),                            calibRoot=os.path.join(repoRoot, "CALIB"))
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


# In[ ]:


if __name__ == "__main__":
    start = time.time()
    main(sys.argv[1:])
    finish = time.time()
    elapsed = finish - start
    print(f"Time elapsed: {elapsed}")

