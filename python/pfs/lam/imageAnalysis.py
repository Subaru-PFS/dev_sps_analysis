import pandas as pd
from astropy.io import fits
from datetime import datetime as dt
import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import center_of_mass
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from ginga.util import iqcalc
from scipy import interpolate
from scipy.interpolate import interp1d
import matplotlib.patches as patches
from matplotlib.colors import LogNorm


from astropy.modeling import models, fitting
from astropy.stats import gaussian_sigma_to_fwhm


# # Fit

class TFocusDf(pd.DataFrame):
    def __init__(self, data, focus='max'):
        pd.DataFrame.__init__(self, data=data, columns=['x', 'y'])
        self.focusMethod = focus

    @property
    def focus(self):
        return self.x[self.indfocus]

    @property
    def indfocus(self):
        if self.focusMethod == 'max':
            return self.y.idxmax()
        else:
            return self.y.idxmin()

    @property
    def vline(self):
        return dict(x=self.focus, ymin=min(self.y), ymax=max(self.y))

class FitGauss1D(pd.DataFrame):
    def __init__(self, amplitude, mean, sigma, offset):
        fwhm = sigma*2*np.log(2 + np.sqrt(3))
        pd.DataFrame.__init__(self, data=[(amplitude, mean, sigma, offset,fwhm)], columns=['amplitude', 'mean', 'sigma', 'offset','fwhm'])


def fitgauss1D(x, y):
    offset = np.median(y)
    fy = y-offset
    amp = np.max(fy)
    mean = x[np.argmin(np.abs(fy-amp))]
    hmean = x[np.argmin(np.abs(fy-amp/2))]
    sig =  np.abs(mean-hmean) / (np.sqrt(2 * np.log(2)))

    popt1, pcov = curve_fit(oneD_Gaussian, x, y, p0=[amp,mean,sig,offset], maxfev=10000)

    newx = np.linspace(np.min(x), np.max(x), 10000)
    data = np.zeros((len(newx), 2))
    data[:,0] = newx     
    data[:,1] = oneD_Gaussian(newx, *popt1)

    return TFocusDf(data), FitGauss1D(*popt1)


def parabola(x, a, b, c):
    return a*x**2 + b*x + c

def fitparabola(x, y, deg=2, focus='min'):
    c = np.polyfit(x, y, deg)
    newx = np.linspace(np.min(x), np.max(x), 10000)
    data = np.zeros((len(newx), 2))
    data[:,0] = newx     
    data[:,1] = np.polyval(c, newx)
        
    return TFocusDf(data=data, focus=focus)


def interpdata(x, y, criteria):
    focus = 'min' if 'fwhm' in criteria else 'max'
    f = interpolate.interp1d(x, y, kind='cubic')
    newx = np.linspace(np.min(x), np.max(x), 1000)
    
    data = np.zeros((len(newx), 2))
    data[:,0] = newx     
    data[:,1] = f(newx)
        
    return TFocusDf(data=data, focus=focus)

def oneD_Gaussian(x, amp, mu, sig, offset):
    return offset + amp * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def openImage(image):
    if type(image) is str:
        hdulist = fits.open(image, "readonly")
        image = hdulist[1].data
        
    return image

def selectRoi(image, cx, cy, roi_size, doBck=False, nRows=5, doPlot=False):
    data = openImage(image)

    indx = int(cy)
    indy = int(cx)
    half = int(roi_size/2)
    
    nRows = nRows if doBck else 0
    broi = np.copy(data[indx - half - nRows: indx + half + 1 + nRows, indy - half : indy + half + 1 ])
    roi = broi[nRows:-nRows,:] if doBck else broi

    edges = np.concatenate([ broi[:nRows,:], broi[-nRows:,:]])

    if doBck:
        for i in range(roi.shape[1]):
            broi[:,i] -= np.median(edges[:,i])

    if doPlot:
        f, ax = plt.subplots(1, 1)
        y, x = center_of_mass(roi)  
        if doBck:
            offset = nRows
            im = broi
            ax.add_patch(patches.Rectangle((-0.5,-0.5), im.shape[1], nRows, fill=False, edgecolor='r'))
            ax.add_patch(patches.Rectangle((-0.5,im.shape[0]-0.5-nRows), im.shape[1], nRows, fill=False, edgecolor='r'))
        else:
            im = roi
            offset = 0

        ax.scatter(x, y+offset, marker='x', s=40, c='w')
        ax.add_patch(patches.Rectangle((roi.shape[1]/2-3, roi.shape[0]/2-3+offset), 5, 5, fill=False, edgecolor='w'))
        ax.add_patch(patches.Rectangle((roi.shape[1]/2-2, roi.shape[0]/2-2+offset), 3, 3, fill=False, edgecolor='m'))
        ax.imshow(im, origin='lower')
            
    return roi

def estimateCOM(image, cx, cy, roi_size=30, doBck=True, nRows=5, seek_size=None):
    """
    estimate center of mass 
    if seek_size is given, recentration will be tried using the seek_size
    if doBck = True it used the rRowns outside of the roi_size to evaluate the bck level and substract it to the roi sub_image
    Return centroid of the center of mass in px
    """
    if seek_size is not None:
        (cx,cy) = estimateCOM(image, cx, cy, roi_size=seek_size, doBck=doBck, nRows=nRows)       
    
    indx = int(cy)
    indy = int(cx)
    half = int(roi_size/2)
    
    outer_data = selectRoi(image, cx, cy, roi_size=roi_size, doBck=doBck, nRows=nRows) 
    # center_of_mass does not support nan so replace with 0
    
    y, x = center_of_mass(outer_data)

    return indy - half + x + 0.5, indx - half + y + 0.5

def getEE(image, cx, cy, ee_size=5, roi_size=100, doBck=False, nRows=5):
    roi = selectRoi(image, cx, cy, roi_size=roi_size, doBck=doBck, nRows=nRows)
    ee = selectRoi(roi, roi.shape[1]/2, roi.shape[0]/2, roi_size=ee_size, doBck=False)
               
    return np.sum(ee)/np.sum(roi), np.sum(roi)

def EELsf(image, px, py, ee_size, roi_size=16, doBck=False):
    roi = selectRoi(image, px, py, roi_size=roi_size, doBck=doBck)
    lsf = roi.sum(axis=1)
    cmass, = center_of_mass(lsf)
    grid = np.arange(len(lsf))
    ngrid = np.arange(len(lsf)-1) + cmass%1
    fint = interp1d(grid, lsf, kind='cubic')
    nlsf = fint(ngrid)
    imax = np.argmax(nlsf)
    half = int(ee_size/2)
    sroi = nlsf[imax-half:imax+half+1]

    return np.sum(sroi)/np.sum(nlsf)

def getPeakData(image, cx, cy, EE=None,roi_size=30, seek_size=None, \
                doPlot=False,com=True,doFit=False,\
                doBck=False,doLSF=False,fwhm_radius=10,fwhm_method='gaussian',**kwargs):
    
    image = openImage(image)
    
    comx, comy = estimateCOM(image, cx, cy, roi_size=roi_size, seek_size=seek_size)
    if doFit :
        calc = iqcalc.IQCalc(None)
        [obj] = calc.evaluate_peaks([(comx-0.5,comy-0.5)], image, fwhm_radius=fwhm_radius,fwhm_method=fwhm_method,\
                                      cb_fn=None, ev_intr=None, **kwargs)
        obj['oid_x'] = comx
        obj['oid_y'] = comy
        obj['objx'] = obj['objx'] + 0.5
        obj['objy'] = obj['objy'] + 0.5

        obj["px"] = obj['oid_x'] if com else obj['objx']
        obj["py"] = obj['oid_y'] if com else obj['objy']
    else:
        obj ={"px": comx,
              "py": comy,
              "oid_x": comx,
              "oid_y": comy}

    EE = [3,5] if EE is None else EE
    for ee in EE:
        obj["EE%d"%ee],obj["TotEE%d"%ee]  = getEE(image, obj["px"], obj["py"], ee_size=ee, roi_size=roi_size, doBck=doBck)
        if doLSF: obj["EL%d"%ee] = EELsf(image, obj["px"], obj["py"], ee_size=ee, roi_size=roi_size, doBck=doBck)
    if doLSF: obj["fwhm_lsf"] = getLSF(image, obj["px"], obj["py"], half_width=roi_size/2, half_spec = roi_size/2, doPlot=doPlot)
        
    if doPlot:
        selectRoi(image, obj["px"], obj["py"], roi_size=roi_size, doBck=doBck, nRows=5, doPlot=True)
        
    return dict(obj)

def getLSF(filename, lsf_x, lsf_y, half_width=5, half_spec = 5, doPlot=False):
    if type(filename) is list and len(filename) == 1 :
        filename = filename[0]

    if type(filename) is str:     
        hdulist = fits.open(filename, "readonly")
        data = hdulist[1].data
    elif type(filename) is np.ndarray:
        data = filename
    else :
        raise Exception("getLSF: file or image type not correct")
    
    # Fit the data using a Gaussian
    g_init = models.Gaussian1D(amplitude=1000., mean=10, stddev=1.)
    fit_g = fitting.LevMarLSQFitter()
    x = np.arange(2*half_spec)
    y = data[int(lsf_y-half_spec):int(lsf_y+half_spec),int(lsf_x-half_width):int(lsf_x+half_width)].sum(axis=1)
    g = fit_g(g_init, x, y)
    fwhm = g.stddev* gaussian_sigma_to_fwhm
    
    if doPlot:
        # Plot the data with the best-fit model
        plt.figure(figsize=(8,5))
        plt.plot(x, y, 'ko')
        newx = np.arange(np.min(x), np.max(x) + 0.001, 0.001)

        plt.plot(newx, g(newx), label='Gaussian')
        plt.xlabel('Position')
        plt.ylabel('Flux')
        plt.legend(loc=2)
    
    return fwhm


def getRois(image, cx, cy, inner_size=5, outer_size=100, doBck=False, nRows=5):
    indx = int(cy)
    indy = int(cx)
    
    if type(image) is str:
        hdulist = fits.open(image, "readonly")
        image = hdulist[1].data
        
    data = np.copy(image)
    
    halfout = int(outer_size/2)
    halfin = int(inner_size/2)
    
    outer_data = data[indx - halfout : indx + halfout + 1,
                      indy - halfout : indy + halfout + 1]
    
    inner_data = data[indx - halfin : indx + halfin + 1,
                      indy - halfin : indy + halfin + 1]
    
    
    if doBck:
        lowRows = image[indx - halfout-nRows:indx - halfout, indy - halfout : indy + halfout + 1]
        upRows = image[indx + halfout + 1:indx + halfout + 1 +nRows, indy - halfout : indy + halfout + 1]
      
        allRows = np.zeros((2*nRows, outer_data.shape[1]))
        allRows[:nRows,:] = lowRows
        allRows[-nRows:,:] = upRows
        for i in range(outer_size):
            outer_data[:,i] -= np.median(allRows[:,i])
            
    return outer_data, inner_data


def neighbor_outlier_filter(df, column, thres, absolute=False):
    """
    Filter data of the Dataframe column by neighbors comparison
    add a <column>_nbh_flag 
    First and last point are compare the previous and following point respectively 
    This flag can then be used to filter the data:
    df[df.<column>_nbh_flag] return filtered Dataframe so values that are greater than the threshold
    """
    dfc = df.copy()
    dfc.loc[:,f"{column}_nbh_diff"] = dfc[column] - (dfc[column].shift(-1) + dfc[column].shift(1))/2 
    dfc[f"{column}_nbh_diff"].fillna(0, inplace=True)
    dfc[f"{column}_nbh_diff"].iloc[0] = (dfc[column] - dfc[column].shift(-1)).iloc[0]
    dfc[f"{column}_nbh_diff"].iloc[-1] = (dfc[column] - dfc[column].shift(1)).iloc[-1]
    if absolute:
        dfc.loc[:,f"{column}_nbh_flag"] = (abs(dfc[f"{column}_nbh_diff"])<abs(thres))
    else:
        dfc.loc[:,f"{column}_nbh_flag"] = (dfc[f"{column}_nbh_diff"]> thres)
        
    return dfc
