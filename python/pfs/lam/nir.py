# -*- coding: utf-8 -*-
import numpy as np
import os
import matplotlib.colors as mpColors
import pandas as pd
import matplotlib.pyplot as plt
from pfs.lam.imageAnalysis import *
#from pfs.lam.detAnalysis import getFullImageQuality
from datetime import datetime



def getImageArrayFromButler(butler:object, visit:int, arm:str, specId:int) -> object:
    dataId = dict(visit=visit, arm=arm, spectrograph=specId)
    exp = butler.get("calexp", dataId)
    return exp.image.array


def getImageFromButler(butler:object, visit:int, arm:str, specId:int, masked=True) -> object:
    dataId = dict(visit=visit, arm=arm, spectrograph=specId)
    exp = butler.get("calexp", dataId)
    if masked is True:
        return exp.maskedImage
    else :
        return exp.image

def getMedianImage(butler:object ,arm:str, specId:int, visits:list=None, expId:int=None) -> object:
    """
    """
    if expId is None and visits is None:
        raise TypeError("Either visits or expId argument require, but not both")
    if expId is not None:
        visitStart, visitEnd = getVisitRange(expId)
    else:
        visitStart, visitEnd = visits
        
    image_concat = [getImageArrayFromButler(butler, visit, arm, specId) for visit in range(visitStart,visitEnd+1)]
    imgs = np.asarray(image_concat)

    # Take the median over the first dim
    med = np.median(imgs, axis=0)
    
    return med

# equivalent LSST:
# from lsst.afw.math import statisticsStack, MEDIAN
# imageList = [getImageFromButler(butler, visit, arm, specId) for visit in range(visitStart,visitEnd+1)]
# stack = statisticsStack(imageList, MEDIAN)


    
def isRoiBad(calexp:object, X:float, Y:float , roi_size:int=16, getMaskName=False)->bool:
#    badMaskList = ['BAD', 'BAD_FIBERTRACE', 'BAD_FLAT', 'BAD_FLUXCAL', 'BAD_SKY',  'CR', 'FIBERTRACE', 'NO_DATA', 'SAT', 'SUSPECT', 'UNMASKEDNAN']
#    badMaskList = ['BAD', 'CR', 'NO_DATA', 'SAT']
    badMaskList = ['CR', 'NO_DATA', 'SAT']

    allX = np.arange(int(X - roi_size/2), int(X + roi_size/2), 1)
    allY = np.arange(int(Y - roi_size/2), int(Y + roi_size/2), 1)
    allpM = []
    flag = False
    for x in allX:
        for y in allY:
            pixMasks = calexp.mask.getAsString(x,y)
            for m in pixMasks.split(","):
                if m in badMaskList:
                    flag = True
                    allpM.extend(pixMasks.split(","))
    if len(allpM) >1 :
        allpM = set(allpM)
    if getMaskName is True:
        return flag, allpM
    else:
        return flag
    
    
def getSlitOffsets(md:dict, pix2mm:list=[0.035, 0.035], verbose:bool=False) -> tuple:
    """
    Retrieve slit move values using header though Metadata, 
    and convert it to values in detector plane using pix2mm factors
    To be implemented: use pfs.instdata to get pix2mm factors
    md: medatadata from butler : exp.getMetadata()
    """
    yslitOff = md["W_ENFCAY"]
    xslitOff = md["W_ENFCAZ"]
    coeffX, coeffY = pix2mm
    ypix = np.round(-1 * yslitOff / coeffY,2)
    xpix = np.round(xslitOff / coeffX,2)
    if verbose is True:
        print(f"Slit move: z={xslitOff}mm  y={yslitOff}mm so \nDet move:  x={xpix}px y={ypix}px")  
        
    return xpix, ypix



#
# Symbolic names for mask/line colors.  N.b. ds9 supports any X11 color for masks
#
WHITE = "white"
BLACK = "black"
RED = "red"
GREEN = "green"
BLUE = "blue"
CYAN = "cyan"
MAGENTA = "magenta"
YELLOW = "yellow"
ORANGE = "orange"
IGNORE = "ignore"

_defaultMaskPlaneColor = dict(
        BAD=RED,
        CR=MAGENTA,
        EDGE=YELLOW,
        INTERPOLATED=GREEN,
        SATURATED=GREEN,
        DETECTED=BLUE,
        DETECTED_NEGATIVE=CYAN,
        SUSPECT=YELLOW,
        NO_DATA=ORANGE,
        # deprecated names
        INTRP=GREEN,
        SAT=GREEN,
    )
    
def getMaskPlaneColor(name):

    color = _defaultMaskPlaneColor.get(name)

    return color



def plotRoiPeakMaskedImage(exp, peak_list, roi_size=20, raw=False, scale=True, verbose=False, savePlotFile=None, doSave=False):
    
   
    image = exp.image.array
    
    dataArr = exp.mask.getArray()

    maskPlanes = exp.mask.getMaskPlaneDict()
    nMaskPlanes = max(maskPlanes.values()) + 1

    planes = {}                      # build inverse dictionary
    for key in maskPlanes:
        planes[maskPlanes[key]] = key

    planeList = range(nMaskPlanes)

    maskArr = np.zeros_like(dataArr, dtype=np.int32)

    colorNames = ['black']
    #colorGenerator = self.display.maskColorGenerator(omitBW=True)
    for p in planeList:
        color = getMaskPlaneColor(planes[p]) if p in planes else None
        if not color:
            color = 'black'

        colorNames.append(color)

    colors = mpColors.to_rgba_array(colorNames)    

    alphaChannel = 3            # the alpha channel; the A in RGBA
    colors[0][alphaChannel] = 0.0      # it's black anyway
    for i, p in enumerate(planeList):
        if colorNames[i + 1] == 'black':
            alpha = 0.0
        else:
            alpha = 1 - 0.4

        colors[i + 1][alphaChannel] = alpha

    cmap = mpColors.ListedColormap(colors)
    norm = mpColors.NoNorm()


       
     
    plist = pd.read_csv(peak_list) if type(peak_list) is str else peak_list
    # 2023/02/03
    # y a des duplicates et je ne sais pas pourquoi
    plist = plist.drop_duplicates(subset=["fiber", "wavelength"])
    
#    plist = plist.sort_values(["X","Y"], ascending=[True,False])
    # Use X,Y when it the list of peak, px,py otherwise
    ind_x = "px"
    ind_y = "py"
    if raw :
        ind_x = "X"
        ind_y = "Y"        
    half = int(roi_size/2)
    
    nbfiber = len(plist.fiber.unique())
    nbwave = len(plist.wavelength.unique())
    # have a list of fiber from 0 to nbfiber
    listwavestoplot = pd.DataFrame(plist.wavelength.unique(), columns=["wavelength"])

    listfiberstoplot = pd.DataFrame(plist.fiber.unique(), columns=["fiber"])
    
    print("#Fiber= %d and #wavelength= %d"%(nbfiber, nbwave))
    f, axarr = plt.subplots(nbwave, nbfiber,  sharex='col', sharey='row',figsize=(12,8))
#    print(axarr.shape)
    vmin=None
    vmax=None
    
    for (wave, fiber), group in plist.groupby(['wavelength','fiber']):
        k = listwavestoplot[listwavestoplot.wavelength == wave].index.tolist()[0]
        i = listfiberstoplot[listfiberstoplot.fiber == fiber].index.tolist()[0]
        if verbose:
            print(k,i)
            print(f"px {group[ind_x]}    py: {group[ind_y]}")
        #cut_data = image[int(indx-roi_size/2):int(indx+roi_size/2), int(indy-roi_size/2):int(indy+roi_size/2)]
        
        try:
            cut_data = selectRoi(image, group[ind_x], group[ind_y], roi_size=roi_size)
        except:
            print(group[ind_x], group[ind_y])
        maskArr_roi = selectRoi(maskArr, group[ind_x], group[ind_y], roi_size=roi_size)
        dataArr_roi = selectRoi(dataArr, group[ind_x], group[ind_y], roi_size=roi_size)
        if nbwave == 1 and nbfiber == 1:
            axarr.set_title(f"{str(fiber)}, {str(wave)}")
            axarr.imshow(cut_data,interpolation="none", origin="lower", vmin=vmin, vmax=vmax)
        else:
            #axarr[nbwave -1 -k, nbfiber - i -1].set_title(f"{fiber}, {wave:.2f}")
            #axarr[nbwave -1 -k, nbfiber - i -1].label_outer()
            axarr[nbwave -1 -k, nbfiber - i -1].imshow(cut_data,interpolation="none", origin="lower", vmin=vmin, vmax=vmax)
            #axarr[nbwave -1 -k, nbfiber - i -1].set_ylabel(f"{wave:.2f}")
            #axarr[nbwave -1 -k, nbfiber - i -1].set_xlabel(fiber)
            axarr[nbwave -1 -k, nbfiber - i -1].grid(False)
            for i, p in reversed(list(enumerate(planeList))):
                if colors[i + 1][alphaChannel] == 0:  # colors[0] is reserved
                    continue

                bitIsSet = (dataArr_roi & (1 << p)) != 0
                if bitIsSet.sum() == 0:
                    continue

                maskArr_roi[bitIsSet] = i + 1  # + 1 as we set colorNames[0] to black
                axarr[nbwave -1 -k, nbfiber - i -1].imshow(maskArr_roi, origin='lower', interpolation='nearest',
                               cmap=cmap, norm=norm)
                maskArr_roi[:] = 0

            
            

                
    f.subplots_adjust(hspace=0.5,wspace=0.5)

    for ax, wave in zip(axarr[:,0], listwavestoplot.sort_index(ascending=False).wavelength.values) :
            ax.set_ylabel(f"{wave:.2f}", rotation='horizontal', ha='right', fontsize=20)
#            ax.set_xlabel(ax.get_xlabel(), rotation='vertical', fontsize=20)
            ax.set_yticklabels('')
            ax.set_xticklabels('')
            ax.set_frame_on(False)
    for ax, fiber in zip(axarr[-1,:], listfiberstoplot.sort_index(ascending=False).fiber.values):
#            ax.set_ylabel(ax.get_ylabel(), rotation='horizontal', ha='right', fontsize=20)
            ax.set_xlabel(fiber, rotation='vertical', fontsize=20)
            ax.set_yticklabels('')
            ax.set_xticklabels('')
            ax.set_frame_on(False)

    plt.gcf().set_facecolor('w')
    if doSave:
        f.patch.set_alpha(0.5)
        plt.savefig(savePlotFile+f"_roi_all.png")
                  
#    plt.show()
    
    

def MaskedImageQualityToCsv(maskedImage, dataId, peaksList, csv_path=".",\
                      roi_size=16, EE=[3,5], seek_size=None,\
                      com=True, doBck=True, doFit=True, doLSF=False, doSep=False,fullSep=False,\
                      doPlot=False, doPrint=False, \
                      mask_size=50, threshold= 50, subpix = 5 , maxPeakDist=80,\
                      maxPeakFlux=40000, minPeakFlux=2000, experimentId=None):
    """
    Calculate quality (EE, and sep info) for each peak from peaklist for a given visit defined in dataId dict
    butler is required to access the data
    use getImagequality and getImageEncerclEnergy(sep) functions
    imageInfo must a dict: dict(arm="b", spectrograph=2, visit=1234, filename="visit_file_path")
    csv file is <csv_path>\Imquality_{cam}_Exp{experimentId}_{visit}_{date_time}.csv
    Returns a pandas Dataframe and write csv file 
    
    """
    
    exp = maskedImage
    visits = dataId["visits"]
    visit = visits[0]
    cam = f"{dataId['arm']}{dataId['spectrograph']}"
    if experimentId is None:
        try:
            experimentId = get_Visit_Set_Id(visit)
        except:
            try:
                experimentId = get_Visit_Set_Id_fromWeb(visit)
            except:
                experimentId = np.nan
                print(f"Unable to get experimentId from logbook for visit: {visit}")
                raise(f"Unable to get experimentId from logbook for visit: {visit}")

    imageInfo = dict(dataId)
    imageInfo.update(experimentId=experimentId)    
    """
    # Update peaklist centroid use known slit positions in the header
    md = exp.getMetadata()
    # values are good for SM2
    pixOffsets = getSlitOffsets(md, pix2mm=[ 0.033914, 0.036498 ], verbose=True)
    if type(peaksList) is str:
        peaksList = filterPeakList(peaksList, dataId["arm"],butler.queryMetadata('raw', ['lamps'], dataId))
    peaksList["X"] = peaksList["X"] + pixOffsets[0]
    peaksList["Y"] = peaksList["Y"] + pixOffsets[1]    
    """   
    
    data = getFullImageQuality(exp.image.array, peaksList, imageInfo=imageInfo,\
                      roi_size=roi_size, EE=EE, seek_size=seek_size,\
                      com=com, doBck=doBck, doFit=doFit, doLSF=doLSF, doSep=doSep,fullSep=fullSep,\
                      doPlot=doPlot, doPrint=doPrint, \
                      mask_size=mask_size, threshold= threshold, subpix = subpix , maxPeakDist=maxPeakDist,\
                      maxPeakFlux=maxPeakFlux, minPeakFlux=minPeakFlux, calexpMask=exp)
    
    now = datetime.now() # current date and time\n",
    date_time = now.strftime("%Y%m%dT%Hh%M")

    csvName = f"Imquality_{cam}_Exp{experimentId}_{visits[0]}-{visits[1]}_meanClip_{date_time}.csv"
    if not os.path.exists(csv_path):
        os.makedirs(csv_path,exist_ok =True)
    data.to_csv(os.path.join(csv_path, csvName))
    if doPrint:
        print(os.path.join(csv_path, csvName))

    return data
