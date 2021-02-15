import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ginga.util import iqcalc
from scipy.interpolate import interp1d
from astropy.io import fits


from pfs.lam.sacFileHandling import Logbook, constructFilelist
import pfs.lam.imageAnalysis as imeas


def slitFindPeak_save(data, radius=60, threshold=300, com=False, doPrint=False):
    calc = iqcalc.IQCalc(None)
    peaks = calc.find_bright_peaks(data, threshold=threshold, radius=radius)
    objlist = calc.evaluate_peaks(peaks, data, fwhm_radius=radius, cb_fn=None, ev_intr=None, fwhm_method='gaussian')
    if doPrint: 
        print(len(objlist))
       

    objlist = [elem for elem in objlist if elem['fwhm'] > 15]
    objlist = [elem for elem in objlist if (elem['fwhm_x'] > 15) and (elem['fwhm_y'] > 15)]
    objlist = [elem for elem in objlist if threshold < elem['brightness'] < 50000]
    if doPrint: 
        print(f"Object detected after filtering: {len(objlist)}")
    if not objlist:
        raise ValueError('peak has not been properly detected')

    maxi = np.argmax([imeas.getEE(image=data, cx=peak['oid_x'], cy=peak['oid_y'], ee_size=20, roi_size=300) for peak in objlist])
    
    return objlist[maxi]['oid_x'], objlist[maxi]['oid_y']

def slitFindPeak(data, radius=60, threshold=300, com=False, doPrint=False):
    calc = iqcalc.IQCalc(None)
    peaks = calc.find_bright_peaks(data, threshold=threshold, radius=radius)
    objlist = calc.evaluate_peaks(peaks, data, fwhm_radius=radius, cb_fn=None, ev_intr=None, fwhm_method='gaussian')
    if doPrint: 
        print(len(objlist))
       

    objlist = [elem for elem in objlist if elem['fwhm'] > 15]
    objlist = [elem for elem in objlist if (elem['fwhm_x'] > 15) and (elem['fwhm_y'] > 15)]
    objlist = [elem for elem in objlist if threshold < elem['brightness'] < 50000]
    if doPrint: 
        print(f"Object detected after filtering: {len(objlist)}")
    if not objlist:
        print('peak has not been properly detected')
        return np.nan, np.nan
    else:
        maxi = np.argmax([imeas.getEE(image=data, cx=peak['oid_x'], cy=peak['oid_y'], ee_size=20, roi_size=300)[0] for peak in objlist])
    
    return objlist[maxi]['oid_x'], objlist[maxi]['oid_y']

def getPeakData(data, roi_size=150, doPlot=False, com=False, doBck=False, fwhm_radius=60, fwhm_method='gaussian', doPrint=False, **kwargs):
    cx, cy = slitFindPeak(data, doPrint=doPrint)
    if doPrint:
        print(f"cx: {cx:.2f}  cy: {cy:.2f}")
    if np.isnan(cx):
        obj ={"px": np.nan,
              "py": np.nan,
              "oid_x": np.nan,
              "oid_y": np.nan,
              "ee20" : np.nan
             }
        return dict(obj)
    
    return imeas.getPeakData(data, cx, cy, EE=[20], roi_size=roi_size, doPlot=doPlot, com=com, fwhm_radius=fwhm_radius, fwhm_method=fwhm_method, **kwargs)



def getSlitPosFromMove_2(experimentId):
    low = float(Logbook.getParameter(experimentId=experimentId,param='lowBound'))
    up = float(Logbook.getParameter(experimentId=experimentId,param='upBound'))
    nbPos = int(Logbook.getParameter(experimentId=experimentId,param='nbPosition'))
    step = (up - low)/(nbPos - 1)
    
    return np.array([low + step * i for i in range(nbPos)])

def getSlitPosFromMove(experimentId):
    low, up, nbPos = Logbook.getParameter(experimentId=experimentId,param='position').split(',')
    low  = float(low)
    up = float(up)
    nbPos = int(nbPos)
    step = (up - low)/(nbPos - 1)
    
    return np.array([low + step * i for i in range(nbPos)])

def stackedImage2(filelist, ind, duplicate, bck=None):
    sublist = filelist[ind*duplicate:ind*duplicate+duplicate]
    first = sublist[0]
    img = fits.open(first)
    hdr =  img[0].header
    data = img[0].data
    res = np.zeros(data.shape, dtype=np.float32)
    
    for filepath in sublist:
        img = fits.open(filepath)
        res += np.copy(img[0].data)
        
    res = res/duplicate
    
    if bck is not None:
        return hdr, res - bck
        
    return hdr, res

def stackedImage(filelist, ind, duplicate, doBck=False):
    sublist = filelist[ind*duplicate:ind*duplicate+duplicate]
    first = sublist[0]
    img = fits.open(first)
    hdr =  img[0].header
    data = img[0].data
    res = np.zeros(data.shape, dtype=np.float64)
    
    for filepath in sublist:
        img = fits.open(filepath)
        res += np.copy(img[0].data.astype('float64'))
        
    res = res/duplicate
    if doBck:
        res-=np.median(res)
        
    return hdr, res

def getSlitTF(experimentId, com=False, doBck=False, doPlot=False, doPrint=False, head=0, tail=0):
    visitStart, visitEnd = Logbook.visitRange(experimentId=experimentId)
    filelist = constructFilelist(visitStart=visitStart, visitEnd=visitEnd)
    duplicate = int(Logbook.getParameter(experimentId=experimentId, param='duplicate'))
    filelist = filelist[head*duplicate:len(filelist)-tail*duplicate]
    guessedPos = getSlitPosFromMove(experimentId)
  
    res = []
    for i in range(len(filelist) // duplicate):
        if doPrint:
            print(f"{i} => {filelist[i]}")
        hdr, data = stackedImage(filelist=filelist, ind=i, duplicate=duplicate, doBck=doBck)
        peak = getPeakData(data, com=com, doBck=False, doPlot=doPlot, doPrint=doPrint)
        peak['experimentId'] = experimentId
        try:
            fca_x = hdr['FCA_X']
        except KeyError:
            fca_x = guessedPos[i]
            
        peak['fca_x'] = fca_x
        res.append(peak)
        if doPlot:
            plt.show()
        if doPrint:
            print("\n")
    return pd.DataFrame(res)

def getSlitTF2(experimentId, com=False, doBck=False, doPlot=False, bck_expId=None, doPrint=False, head=0, tail=0):
    visitStart, visitEnd = Logbook.visitRange(experimentId=experimentId)
    filelist = constructFilelist(visitStart=visitStart, visitEnd=visitEnd)
    duplicate = int(Logbook.getParameter(experimentId=experimentId, param='duplicate'))
    filelist = filelist[head*duplicate:len(filelist)-tail*duplicate]
    guessedPos = getSlitPosFromMove(experimentId)
    bck_data = None
    if bck_expId is not None:
        bck_visitStart, bck_visitEnd = Logbook.visitRange(experimentId=bck_expId)
        bck_filelist = constructFilelist(visitStart=bck_visitStart, visitEnd=bck_visitEnd)
        bck_duplicate = int(Logbook.getParameter(experimentId=bck_expId, param='duplicate'))
        bck_hdr, bck_data = stackedImage2(filelist=bck_filelist, ind=0, duplicate=bck_duplicate, bck=None)

    res = []
    for i in range(len(filelist) // duplicate):
        hdr, data = stackedImage2(filelist=filelist, ind=i, duplicate=duplicate, bck=bck_data)
        peak = getPeakData(data, com=com, doBck=doBck, doPlot=doPlot, doPrint=doPrint)
        peak['experimentId'] = experimentId
        try:
            fca_x = hdr['FCA_X']
        except KeyError:
            fca_x = guessedPos[i]
            
        peak['fca_x'] = fca_x
        res.append(peak)
    
    return pd.DataFrame(res)


def getFocus(series, criteria, index, corrector=False, doPrint=False):
    if 'EE' in criteria:
        if corrector:
            thfoc, parfit = imeas.fitgauss1D(series[index].values, series[criteria].values)  
        else:
            thfoc = imeas.fitparabola(series[index].values, series[criteria].values, deg=15)  
            
    elif criteria=='fwhm':
        thfoc = imeas.fitparabola(series[index].values, series[criteria].values, deg=3)  
        
    else:
        thfoc, parfit = imeas.fitgauss1D(series[index].values, series[criteria].values) 
    
    return thfoc.rename(index=str, columns={"x": index, "y": criteria})


def fitFocusData(cube, corrector=False, doPlot=False, index='fca_x'):
    thfoc_data = []
    
    for experimentId, series in cube.groupby('experimentId'):
        series = series.dropna()
        thfoc = getFocus(series, 'EE20', index, corrector=corrector)
        #for criteria in ['brightness', 'fwhm']:
        #    thfoc[criteria] = getFocus(series, criteria, index, corrector=corrector)[criteria]

        thfoc['px'] = np.interp(thfoc[index], series[index], series['px'])
        thfoc['py'] = np.interp(thfoc[index], series[index], series['py'])
        thfoc['experimentId'] = experimentId
        thfoc_data.append(thfoc)
        
    thfoc_data = pd.concat(thfoc_data)

    if doPlot:
        kwargs = dict(grid=True, figsize=(14,10), legend=True, subplots=True)
        #criterias = ['EE20', 'brightness', 'fwhm']
        criterias = ['EE20']
        for experimentId, fit in thfoc_data.groupby('experimentId'):
            raw = cube.query("experimentId==%d"%(experimentId))
            axes = fit.set_index(index)[criterias].plot(**kwargs)
            for i, criteria in enumerate(criterias):
                axes[i].plot(raw[index].values, raw[criteria].values, 'o')
                
    return thfoc_data


def getFocusModel(fitdata, index='fca_x'):
    data = []
    for experimentId, series in fitdata.groupby('experimentId'):
        series = series.dropna()
        for criteria in ['EE20']: #, 'brightness', 'fwhm']:
            ixmax = series[criteria].idxmax() if criteria !='fwhm' else series[criteria].idxmin()
            focus = series[index][ixmax]
            px = series.px[ixmax]
            py = series.py[ixmax]
            mat = [experimentId, criteria, px, py, focus]
            data.append(tuple(mat))
    
    return pd.DataFrame(data, columns=['experimentId', 'criteria', 'px', 'py', index])
