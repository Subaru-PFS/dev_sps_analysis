# ================================================================================================================================= #
# ======================================== 	Description of slitFocusAnalysis.py 	=============================================== #
# ================================================================================================================================= #
# Author and Editor :                                                                                                               #
#   									                                                                                            #
#    - MADEC Fabrice                : February 2021                                                                                 #
#    - LHOUSSAINE BEN BRAHIM Romain : March 18, 2021                                              								    #
#                                                                                                                                   #
# Purpose :                                                                                                                         #
#   									                                                                                            #
#    - Gather all the functions used for the slit alignment : in particular for the slit through focus experiment                   #
#   									                                                                                            #
# Functions :  																					                                    #
#   									                                                                                            #
#   		- slitFindPeak							                                                                                #
#   		- getPeakData							                                                                                #
#   		- getSlitPosFromMove							                                                                        #
#   		- stackedImage							                                                                                #
#   		- getSlitTF							                                                                                    #
#           - getFocus																		                                        #                          
#           - fitFocusData																		                                    #
#           - getFocusModel																		                                    #
# ================================================================================================================================= #
# ================================================ Librairies ===================================================================== #
# ================================================================================================================================= #
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ginga.util import iqcalc
from scipy.interpolate import interp1d
from astropy.io import fits
from pfs.lam.sacFileHandling import Logbook, constructFilelist
import pfs.lam.imageAnalysis as imeas
# ================================================================================================================================= #


# ================================================================================================================================= #
# ======================================================== slitFindPeak =========================================================== #
# ================================================================================================================================= #
def slitFindPeak(data, radius=60, threshold=300,roi_size=150,fwhm_method='gaussian',com=False, doPrint=False): 
    """Find the maximum peak of brightness in a given image of slit through focus
        After the detection, add the detected object with his informations into a dataframe 
        
        Parameters
        ----------
        data : `dataframe`
            Path name from the LSST data butler. Besides the usual directory and
            filename with extension, this may include a suffix with additional
            characters added by the butler.
        radius : `int`
            Description :
        threshold : `int`
            Description :
        EE : `int`
            Description :
        roi_size : `int`
            Description :
        fmwh_method : `str`
            Description :
        com : `bool`
            Description :
        doPrint : `bool`
            Description :
            
        Returns
        -------
        out : `floats`
            Two values of maximum peak coordinates into all the detected objects.
        ------
        :Example:
 
        >>> cx, cy = slitFindPeak(data, doPrint=doPrint)
        
        -------
        .. seealso:: getPeakData
        .. warning:: There can be some anomalies in the data when you hide the different filter for the objlist parameter.
        .. note:: None
        .. todo:: None
    """
    
    doPlot = False
    maxi = 0
    calc = iqcalc.IQCalc(None)
    peaks = calc.find_bright_peaks(data, threshold=threshold, radius=radius)
    objlist = calc.evaluate_peaks(peaks, data, fwhm_radius=radius, cb_fn=None, ev_intr=None, fwhm_method=fwhm_method)
        
    if doPrint: 
        print(len(objlist))

    objlist = [elem for elem in objlist if (elem['fwhm'] > 15)
                                        and (elem['fwhm_x'] > 15) 
                                        and (elem['fwhm_y'] > 15)
                                        and (0 < elem['objx'] < 1940)
                                        and (0 < elem['objy'] < 1460)
                                        and (threshold < elem['brightness'] < 50000)]
        
    for obj in objlist :
        plt.scatter(obj['oid_x'],obj['oid_y'],s=80,c='red',marker='x',label='peak', edgecolors=None)
        
    if doPlot:
        plt.show()

    if doPrint: 
        print(f"Object detected after filtering: {len(objlist)}")
        
    if not objlist:
        print('peak has not been properly detected')
        obj ={"px": np.nan,
                "py": np.nan,
                "oid_x": np.nan,
                "oid_y": np.nan,
                "EE20" : np.nan
                }
        dict(obj)
        return np.nan, np.nan
    
    else :
        maxi = np.nanargmax([imeas.getEE(image=data,doBck=False, cx=peak['oid_x'], cy=peak['oid_y'], ee_size=20, roi_size=300)[0] for peak in objlist])
        
    return objlist[maxi]['oid_x'], objlist[maxi]['oid_y']
# ================================================================================================================================= #


# ================================================================================================================================= #
# ==================================================== getPeakData ================================================================ #
# ================================================================================================================================= #
def getPeakData(data, radius = 60, threshold = 300, roi_size=150, doPlot=False, com=False, doBck=False, fwhm_radius=60, fwhm_method='gaussian', doPrint=False, **kwargs):
    
    """Get the fwhm for a given peak
        TGet the corrdinates of a peak in a given image, and compute th fwhm for the given method
        
        Parameters
        ----------
        data : `DataFrame`
            Description :
        roi_size : `int`
            Description :
        doPlot : `bool`
            Description :
        com : `bool`
            Description :
        doBck : `bool`
            Description :
        fwhm_radius : `int`
            Description :
        fwhm_method : `bool`
            Description :
        doPrint : `bool`
            Description :
        **kwargs : `**`
            Description :
        
        Returns
        -------
        out : `float`
            The value of fwhm.
        ------
        :Example:
 
        >>> peak = getPeakData(data, com=com, doBck=doBck, doPlot=doPlot, doPrint=doPrint)

        -------
        .. seealso:: getSlitTF, slitFindPeak
        .. warning:: None
        .. note:: None
        .. todo:: None
    """
    
    cx, cy = slitFindPeak(data, radius=radius, threshold=threshold,roi_size=roi_size,fwhm_method=fwhm_method,com=com, doPrint=doPrint)
    
    if doPrint:
        print(f"cx: {cx:.2f}  cy: {cy:.2f}")
    if np.isnan(cx):
        obj ={"px": np.nan,
              "py": np.nan,
              "oid_x": np.nan,
              "oid_y": np.nan,
              "EE20" : np.nan
             }
        return dict(obj)
    
    return imeas.getPeakData(data, cx, cy, EE=[20], roi_size=roi_size, doBck=False, doPlot=doPlot, com=com, fwhm_radius=fwhm_radius, fwhm_method=fwhm_method, **kwargs)
# ================================================================================================================================= #



# ================================================================================================================================= #
# ================================================ getSlitPosFromMove ============================================================= #
# ================================================================================================================================= #
def getSlitPosFromMove(experimentId):
    """Get the values of the low and top positions given for the slit through focus analysis, with the number of positions between the two positions.
        Return a array with a size of (low + step * i for i in range(nbPos)) values.

        Parameters
        ----------
        experimentID : `list`
            Description :
        hdu : `int`
            Description :
        flags : `int`
            Description :
        Returns
        -------
        out : `np.array`
            Empty array
        ------
        :Example:
 
        >>> guessedPos = getSlitPosFromMove(experimentId)

        -------
        .. seealso:: getSlitTF
        .. warning:: None
        .. note:: None
        .. todo:: None
    
    """

    low, up, nbPos = Logbook.getParameter(experimentId=experimentId,param='position').split(',')
    low  = float(low)
    up = float(up)
    nbPos = int(nbPos)
    step = (up - low)/(nbPos - 1)
    
    return np.array([low + step * i for i in range(nbPos)])
# ================================================================================================================================= #


# ================================================================================================================================= #
# ================================================ stackedImage =================================================================== #
# ================================================================================================================================= #
def stackedImage(filelist, ind, duplicate, doBck=False):
    """Gives header and data for a given file list
        In option, the data can be given with the background retrived (if the artificial background corrector is set as True, or some background Ids is given)
        
        Parameters
        ----------
        filelist : `list`
            Description :
        ind : `list`
            Description :
        duplicate : `list`
            Description :
        Returns
        -------
        out : `pfs.drp.stella.FiberTraceSet`
            Traces read from FITS.
        ------    
        :Example:
 
        >>>  bck_hdr, bck_data = stackedImage(filelist=bck_filelist, ind=0, duplicate=bck_duplicate, bck=None)
        >>>  hdr, data = stackedImage(filelist=filelist, ind=i, duplicate=duplicate, doBck=doBck, bck = bck_data)

        -------
        .. seealso:: getSlitTF
        .. warning:: None
        .. note:: The previous version was considering res as float32 np.array, and was retriving the background after.
        .. todo:: None
    
    """
    
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
        
    if doBck is not None and type(doBck) is not bool:
        return hdr, res - bck
    
    elif doBck == True :
        res-=np.median(res)
        return hdr, res
    
    else :
        return hdr, res
# ================================================================================================================================= #


# ================================================================================================================================= #
# ================================================ getSlitTF ====================================================================== #
# ================================================================================================================================= #
def getSlitTF(experimentId, com=False, doBck=False, doPlot=False, doPrint=False, head=0, tail=0,threshold = 300 ,radius = 60  ,roi_size = 150,fwhm_radius = 60,fwhm_method = 'gaussian'):

    """Get the righ focus for a given fiber, with a given experimentID.
        With a given series
        
        Parameters
        ----------
        experimentID : `list`
            Description :
        com : `bool`
            Description :
        doBck : `bool` or `NoneType` or `list` 
            Description :
        doPlot : `bool`
            Description :
        doPrint : `bool`
            Description :
        head : `int`
            Description :
        tail : `int`
            Description :
        threshold : `int`
            Description :
        radius : `int`
            Description :
        roi_size : `int`
            Description :
        EE : `int`
            Description :
        fmwh_radius : `int`
            Description :
        fmwh_method : `str`
            Description :
        Returns
        -------
        out : `DataFrame`
            A dataframe of all the objects detected on the different images
        ------
        :Example:
 
        >>> for experimentId in experimentIds:
                dfs.append(getSlitTF(experimentId=experimentId, com=com, head=head, tail=tail, doBck=doBck,doPrint=False, doPlot=False))
            cube = pd.concat(dfs)
        >>> cube
                objx	objy	pos	oid_x	oid_y	fwhm_x	fwhm_y	fwhm	fwhm_radius	brightness	...	y	skylevel	background	px	py	EE20	TotEE20	experimentId	fiber	fca_x
        0	912.218014	818.003217	0.999102	908.352143	820.844504	82.226640	111.849192	98.161760	60.0	266.146016	...	820.0	40.0	0.0	908.352143	820.844504	0.047722	2219589.0	447	engtopend	-1.500000
        1	911.789595	818.875585	0.999084	911.691810	820.870137	82.227511	100.134550	91.619571	60.0	281.462116	...	820.0	40.0	0.0	911.691810	820.870137	0.050955	2206778.0	447	engtopend	-1.432203
        2	917.870636	819.142704	0.999078	915.094960	820.520054	74.969898	88.877204	82.218133	60.0	322.873802	...	820.0	40.0	0.0	915.094960	820.520054	0.057535	2238367.0	447	engtopend	-1.364407
        3	918.610294	818.798968	0.999086	918.835311	820.192194	74.508485	83.655742	79.214258	60.0	335.996653	...	819.0	40.0	0.0	918.835311	820.192194	0.062846	2196874.0	447	engtopend	-1.296610
        4	924.042145	819.123780	0.999079	922.065059	820.561443	66.678236	78.869550	73.028738	60.0	390.806635	...	820.0	40.0	0.0	922.065059	820.561443	0.069752	2248796.0	447	engtopend	-1.228814

        -------
        .. seealso:: getPeakData, stackedImage, getSlitPosFromMove, slitFindPeak
        .. warning:: if doBck = False, and there is no bck_Id given in the function : an error can appear because no objects was detected.
        .. note:: Working on scientif fibers take more time to compute data than engineering fibers
        .. todo:: Find a solution bout the do_bck/bck_ID problem
    """
        
    visitStart, visitEnd = Logbook.visitRange(experimentId=experimentId)
    filelist = constructFilelist(visitStart=visitStart, visitEnd=visitEnd)
    duplicate = int(Logbook.getParameter(experimentId=experimentId, param='duplicate'))
    filelist = filelist[head*duplicate:len(filelist)-tail*duplicate]
    guessedPos = getSlitPosFromMove(experimentId)
    fiberId = Logbook.getParameter(experimentId, 'fiber', doRaise=False)

    bck_data = None
    if doBck is not None and type(doBck) is not bool:
        bck_visitStart, bck_visitEnd = Logbook.visitRange(experimentId=bck_expId)
        bck_filelist = constructFilelist(visitStart=bck_visitStart, visitEnd=bck_visitEnd)
        bck_duplicate = int(Logbook.getParameter(experimentId=bck_expId, param='duplicate'))
        bck_hdr, bck_data = stackedImage(filelist=bck_filelist, ind=0, duplicate=bck_duplicate, doBck=None)
    
    res = []
    for i in range(len(filelist) // duplicate):
        
        if doBck is not None and type(doBck) is not bool:
            hdr, data = stackedImage(filelist=filelist, ind=i, duplicate=duplicate, doBck = bck_data)
        else :
            hdr, data = stackedImage(filelist=filelist, ind=i, duplicate=duplicate, doBck = doBck)

        if doPlot:
            plt.imshow(data, origin = 'lower', cmap = 'gray', vmin=300, vmax=800)
        
        peak = getPeakData(data, radius = radius, threshold = threshold, roi_size=roi_size, doPlot=doPlot, com=com, doBck=doBck, fwhm_radius=fwhm_radius, fwhm_method=fwhm_method, doPrint=doPrint)
        
        peak['experimentId'] = experimentId
        peak['fiber'] = fiberId

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
# ================================================================================================================================= #


# ================================================================================================================================= #
# ================================================ getFocus ============================================#========================== #
# ================================================================================================================================= #
def getFocus(series, index, corrector=False, doPrint=False):
    
    """Used in fitFocusData
        Sort data going out of getSlitTF with a fitted function (gauss 1D or parabola)
        
        Parameters
        ----------
        series : `dataFrame`
            Path name from the LSST data butler. Besides the usual directory and
            filename with extension, this may include a suffix with additional
            characters added by the butler.
        corrector : `bool`
            Description : Fit 1D Gauss
        doPrint : `bool`
            Description : 
        Returns
        -------
        out : `DataFrame`
            dataframe with two columns : index & criteria
        -------
        :Example:
 
        >>> thfoc = getFocus(series, criteria, index, corrector=corrector)

        -------
        .. seealso:: fitFocusData
        .. warning:: None
        .. note:: dedicated to fitFocusData
        .. todo:: None
    """
    
  
    if corrector:
        thfoc, parfit = imeas.fitgauss1D(series[index].values, series['EE20'].values)  
    else:
        thfoc = imeas.fitparabola(series[index].values, series['EE20'].values, deg=15)  
    
    return thfoc.rename(index=str, columns={"x": index, "y": 'EE20'})
# ================================================================================================================================= #


# ================================================================================================================================= #
# ================================================ fitFocusData =================================================================== #
# ================================================================================================================================= #
def fitFocusData(cube, corrector=False, doPlot=False, index='fca_x'):
    """Sort a given dataframe in order to have for each fiber, all the detected peaks
        This function catch a dataframe and sort all the peaks detected on each fibers (and associated experiment Id). We find also the position of the slit and the value of ensquarred energy.
        
        Parameters
        ----------
        Through focus data : `Data Frame`
            Path name from the LSST data butler. Besides the usual directory and
            filename with extension, this may include a suffix with additional
            characters added by the butler.
        corrector  : `bool`
            Part of the ``FitsCatalogStorage`` API, but not utilised.
        doPlot : `bool`
            Part of the ``FitsCatalogStorage`` API, but not utilised.
        index : `str`
            Part of the ``FitsCatalogStorage`` API, but not utilised.


        Returns
        -------
        out : `Dataframe`
            Data Frame sorted
        ------
        :Example:
 
        >>> thfocModel= fitFocusData(cube2[cube2.EE20_nbh_flag], corrector=corrector, doPlot=doPlot, criteria = criteria)
        >>> thfocModel
        0	-1.500000	0.047827	908.352143	820.844504	447	engtopend
        1	-1.499654	0.047811	908.369177	820.844634	447	engtopend
        2	-1.499308	0.047796	908.386211	820.844765	447	engtopend
        3	-1.498963	0.047781	908.403245	820.844896	447	engtopend
        4	-1.498617	0.047767	908.420279	820.845027	447	engtopend
        ...	...	...	...	...	...	...
        9995	2.498508	0.023265	963.793591	922.788602	450	engbotend
        9996	2.498881	0.023227	963.772601	922.790234	450	engbotend
        9997	2.499254	0.023189	963.751610	922.791865	450	engbotend
        9998	2.499627	0.023150	963.730620	922.793497	450	engbotend
        9999	2.500000	0.023110	963.709630	922.795129	450	engbotend

        ------
        .. seealso:: getSlitTF, getFocusModel
        .. warning:: This is a completly useless function. Use it only in a 
                      tutorial unless you want to look like a fool.
        .. note:: You may want to use a lambda function instead of this.
        .. todo:: Delete this function. Then masturbate with olive oil.
    
    """
    
    thfoc_data = []
    
    for experimentId, series in cube.groupby('experimentId'):
        series = series.dropna()
        thfoc = getFocus(series, index, corrector=corrector)

        thfoc['px'] = np.interp(thfoc[index], series[index], series['px'])
        thfoc['py'] = np.interp(thfoc[index], series[index], series['py'])
        thfoc['experimentId'] = experimentId
        thfoc['fiber'] = series['fiber'].unique()[0]

        thfoc_data.append(thfoc)
        
    thfoc_data = pd.concat(thfoc_data)

    if doPlot:
        kwargs = dict(grid=True, figsize=(14,10), legend=True, subplots=True)
        
        for experimentId, fit in thfoc_data.groupby('experimentId'):
            raw = cube.query("experimentId==%d"%(experimentId))
            axes = fit.set_index(index)[['EE20']].plot(**kwargs)
            for i, criteria in enumerate(['EE20']):
                axes[i].plot(raw[index].values, raw[criteria].values, 'o')
                
    return thfoc_data
# ================================================================================================================================= #


# ================================================================================================================================= #
# ================================================ getFocusModel ================================================================== #
# ================================================================================================================================= #
def getFocusModel(fitdata, index='fca_x'):
    
    """Get the righ focus for a given fiber.
        Resume in a dataframe with the right focus, for each fibers, and his posistion.
        Parameters
        ----------
        Fit focus data : `Dataframe`
            Description : 
        index : `str`
            Description : 
        Returns
        -------
        out : `DataFrame`
            For each fiber(with the associated ID, Criteria, positions on x & y), the compilated right focus fca_x and his coodinates.
        ------    
        :Example:
 
        >>> focusModel = getFocusModel(thfocModel)
        >>> focusModel
            experimentId	criteria	px	py	fca_x	fiber
        0	465	EE40	1038.690096	833.223406	0.542614	scitopend
        1	466	EE40	964.113287	856.918097	0.359796	scitopmid
        2	467	EE40	886.065193	846.060377	0.373393	scibotmid
        3	468	EE40	747.334953	840.743436	0.556201	scibotend

        ------
        .. seealso:: fitFocusData(cube, corrector=False, doPlot=False, index='fca_x', criteria = 'EE20'),getFocus(series, criteria, index, corrector=False, doPrint=False), mul()
        .. warning:: Efficient if used after fitFocusData
        .. note:: focusModel.fca_x.values allow to access directly the values of fca_x
        .. todo:: None
    """
    
    data = []
    for experimentId, series in fitdata.groupby('experimentId'):
        series = series.dropna()
        fiber = series['fiber'].unique()[0]
        for i in ['EE20']:

            ixmax = series[i].idxmax()
            focus = series[index][ixmax]
            px = series.px[ixmax]
            py = series.py[ixmax]
            mat = [experimentId, i, px, py, focus, fiber]
            data.append(tuple(mat))
    
    return pd.DataFrame(data, columns=['experimentId', 'EE20', 'px', 'py', index,'fiber'])
# ================================================================================================================================= #

