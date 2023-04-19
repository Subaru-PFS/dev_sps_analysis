from pfs.lam.imageAnalysis import *
import scipy.interpolate
import matplotlib.gridspec as gridspec
from scipy import spatial
from datetime import datetime
from pfs.lam.fileHandling import *
from pfs.lam.analysisPlot import plotRoiPeak
from pfs.lam.opdb import get_Visit_Set_Id, get_Visit_Set_Id_fromWeb
from pfs.lam.linePeaksList import removeClosePeak, removeFluxPeak
from pfs.lam.nir import isRoiBad
from pfs.lam.sep import *



def getImageQuality(image, peak_list, roi_size=20, EE=[3,5], seek_size=None,\
                    com=True, doPlot=False, scalePlot=False, doBck=False, doFit=False, doLSF=False,\
                    calexpMask=None):
    """
    Calulate Ensquared Energy in EExEE px for all peak given in peak_list for a given image
    Returns a pandas Dataframe
    """

    if type(image) is list and len(image) == 1 :
        image = image[0]

    if type(image) is str:     
        hdulist = fits.open(image, "readonly")
        prihdr = hdulist[0].header
        image = hdulist[1].data
 
    plist = pd.read_csv(peak_list) if type(peak_list) is str else peak_list
    objlist=[]
    for index, row in plist.iterrows():
        cx = row["X"]
        cy = row["Y"]
        wave = row["wavelength"] if "wavelength" in plist.columns else np.nan
        try:
            obj = getPeakData(image, cx,cy, EE=EE, roi_size=roi_size, seek_size=seek_size, com=com, doBck=doBck, doFit=doFit, doLSF=doLSF)
            obj["peak"] = row["peak"]
            obj["fiber"] = row["fiber"]
            obj["wavelength"] = wave
            obj["peak_element"] = row["element"] if "element" in plist.columns else np.nan
            obj["peak_exptime"] = row["exptime"] if "exptime" in plist.columns else np.nan
            obj["peak_lamp"] = row["lamp"] if "lamp" in plist.columns else np.nan
            obj["peak_dcb_wheel"] = row["dcb_wheel"] if "dcb_wheel" in plist.columns else np.nan
            
            # new 2023/01/26 add defect flag from DRP
            if calexpMask is not None:
                try:
                    obj["DRP_flag"] = isRoiBad(calexpMask, obj["px"], obj["py"], roi_size=roi_size, getMaskName=False)
                except Exception as e:
                    print("DRP FLAG FAILED: ",str(e), "cx:%i, cy:%i"%(cx,cy))
                
            objlist.append(obj)
        except Exception as e:
            print("getPeakData FAILED: ",str(e), "cx:%i, cy:%i"%(cx,cy))
            objlist.append(dict(peak=row["peak"], fiber=row["fiber"], wavelength=wave))
    #print(cx,cy)

    mdata = pd.DataFrame(objlist)
    if doPlot :
        plotRoiPeak(image, mdata, roi_size, scale=scalePlot)

    return mdata



# get Encercled Energy using SEP

def getImageEncerclEnergy(image, peak_list, roi_size=20, EE=[3,5], seek_size=None,\
                          mask_size=50, threshold= 50, subpix = 5 , maxPeakDist=80, maxPeakFlux=40000, minPeakFlux=2000, doPlot=False, scalePlot=False, doBck=False, doEE=False):

    if type(image) is list and len(image) == 1 :
        image = image[0]

    if type(image) is str:     
        hdulist = fits.open(image, "readonly")
        prihdr = hdulist[0].header
        image = hdulist[1].data
 
    if peak_list is not None :
        plist = pd.read_csv(peak_list) if type(peak_list) is str else peak_list
        objlist=[]

        for index, row in plist.iterrows():
            cx = row["X"]
            cy = row["Y"]
            try:
                #obj = getFluxPeakDataSep(image, cx,cy, EE=EE, roi_size=roi_size,mask_size=mask_size, subpix=subpix, doBck=doBck)
                obj = getPeakDataSep2(image, cx,cy, EE=EE, roi_size=roi_size, seek_size=seek_size, mask_size=mask_size,\
                                     subpix=subpix, doBck=doBck, doEE=doEE, threshold=threshold)
                obj["cx"] = cx
                obj["cy"] = cy
                obj["peak"] = row["peak"]
                obj["fiber"] = row["fiber"]
                obj["wavelength"] = row["wavelength"]

                objlist.append(obj)
            except Exception as e:
                print(str(e), "cx:%i, cy:%i"%(cx,cy))
                objlist.append(dict(peak=row["wavelength"], fiber=row["fiber"]))

        mdata = pd.concat(objlist)
    else : # do it on the whole image so every peak
        obj = sep.extract(image, threshold)
        df = pd.DataFrame(obj, columns=obj.dtype.names)
        df = df[["flux", "peak", "x", "y", "flag", "npix"]]
        df = df.rename(columns={'x': 'px','y': 'py', 'peak': 'brightness'})
        df = removeClosePeak(df, dist=maxPeakDist, doPlot=doPlot)
        df = removeFluxPeak(df, fmax=maxPeakFlux, fmin=minPeakFlux, doPlot=doPlot)
        if doEE :
            for (px,py),a in df.groupby(["px", "py"]):
                flux_tot, flux_tot_err, flux_tot_flag = sep.sum_circle(image, df['px'], df['py'],
                                                 roi_size/2., err=None, gain=1.0, subpix=subpix)
                EE = [3,5] if EE is None else EE
                for ee in EE:
                    df["ECE%d"%ee], df["fluxErr"], flag = sep.sum_circle(image, df['px'], df['py'],
                                                 ee/2., err=None, gain=1.0, subpix=subpix)
                    df["ECE%d"%ee] = df["ECE%d"%ee] / flux_tot
        mdata = df
        
        

    if doPlot :
        plt_data = mdata[["wavelength", "fiber", "px", "py"]]
        plt_data = plt_data.rename(columns={'px': 'X','py': 'Y'})

        plotRoiPeak(image, plt_data, roi_size, scale=scalePlot)


    return mdata



def getStatIM(dframe, par="EE5", thresold=0.9):
    z = dframe[par]
    zs = dframe[dframe[par]>thresold][par]
    print("%.f %% %s peaks >%s"%(100*len(zs)/len(z),par, thresold))


    
def getFullImageQuality(image, peaksList, roi_size=16, seek_size=None, imageInfo=None,\
                      com=True, doBck=True, EE=[3,5], doFit=True, doLSF=False, \
                      doSep=False, fullSep=False,\
                      doPlot=False, doPrint=False, csv_path=None, \
                      mask_size=50, threshold= 50, subpix = 5 , maxPeakDist=80,\
                      maxPeakFlux=40000, minPeakFlux=2000,\
                      calexpMask=None):
    """
    Calculate quality (EE, and sep info) for each peak from peaklist for a given image
    use getImagequality and getImageEncerclEnergy(sep) functions
    imageInfo must a dict: dict(arm="b", spectrograph=2, visit=1234, filename="visit_file_path")
    Returns a pandas Dataframe
    """
    data = getImageQuality(image, peaksList, roi_size=roi_size, seek_size=seek_size, EE=EE, com=com, \
                           doPlot=doPlot, doBck=doBck, doFit=doFit, doLSF=doLSF, calexpMask=calexpMask)
    if doSep:
        dsep = getImageEncerclEnergy(image, peaksList, roi_size=roi_size, EE=EE,\
        mask_size=mask_size, threshold= threshold, subpix = subpix ,\
        maxPeakDist=maxPeakDist, maxPeakFlux=maxPeakFlux, minPeakFlux=minPeakFlux,\
        doPlot=doPlot, doBck=doBck, doEE=fullSep)
        dsep = dsep.add_prefix("sep_")
        dsep = dsep.rename(columns={'sep_wavelength': 'wavelength','sep_fiber': 'fiber'})
        
        data = data.merge(dsep, on=["wavelength","fiber"])

        
    if imageInfo is not None:
        visitfilepath = imageInfo["filename"] 
        visit = imageInfo["visit"] 
        cam = f"{imageInfo['arm']}{imageInfo['spectrograph']}"
        if "experimentId" in imageInfo.keys():
            experimentId = imageInfo["experimentId"]
        else:
            try:
                experimentId = get_Visit_Set_Id(visit)
            except:
                try:
                    experimentId = get_Visit_Set_Id_fromWeb(visit)
                except:
                    experimentId = np.nan
                    raise(f"Unable to get experimentId from logbook for visit: {visit}")



            # populate dataframe with keywords information
        __, fname = os.path.split(visitfilepath)

        data["filename"] = fname
        data["visit"] = visit

        if type(peaksList) is str: 
            data["peaklist"] = peaksList        
        
        if getFitsKey(visitfilepath, 'W_XM1POS', doRaise=False) is not np.nan :
            data["xm1pos"] = np.float(getFitsKey(visitfilepath, 'W_XM1POS'))
            data["xm2pos"] = np.float(getFitsKey(visitfilepath, 'W_XM2POS'))
            data["xm3pos"] = np.float(getFitsKey(visitfilepath, 'W_XM3POS'))
        else :
            data["xm1pos"] = np.float(getFitsKey(visitfilepath, 'HIERARCH W_XCU_MOTOR1_MICRONS'))
            data["xm2pos"] = np.float(getFitsKey(visitfilepath, 'HIERARCH W_XCU_MOTOR2_MICRONS'))
            data["xm3pos"] = np.float(getFitsKey(visitfilepath, 'HIERARCH W_XCU_MOTOR3_MICRONS'))
           
        
        data["motor1"] = data["xm1pos"] 
        data["motor2"] = data["xm2pos"]
        data["motor3"] = data["xm3pos"]

        fcax = np.float(getFitsKey(visitfilepath, 'W_ENFCAX', doRaise=False))
        fcay = np.float(getFitsKey(visitfilepath, 'W_ENFCAY', doRaise=False))
        fcaz = np.float(getFitsKey(visitfilepath, 'W_ENFCAZ', doRaise=False))

        # OneChannel back compatiblity
        fcax = np.float(getFitsKey(visitfilepath, 'HIERARCH W_FCA_FOCUS', doRaise=False)) if np.isnan(fcax) else fcax
        data['fcaFocus'] = fcax
        data['fcaX'] = fcax
        data['fcaY'] = fcay
        data['fcaZ'] = fcaz

        ccdTemp = np.float(getFitsKey(visitfilepath, 'W_XTDET1', doRaise=False))
        data['ccdTemp'] = ccdTemp
        detBoxTemp = np.float(getFitsKey(visitfilepath, 'W_XTDBOX', doRaise=False))
        data['detBoxTemp'] = detBoxTemp

        data['cam'] = cam
        data['obsdate'] = getFitsKey(visitfilepath, 'DATE-AVG')
        data['experimentId'] = experimentId
        
        data["exptime"] = round(getFitsKey(visitfilepath, 'exptime'),1)
        data["lamp"] = getArcLampForNist(None,fitsfile=visitfilepath, strict=False)
        data["dcb_wheel"] = getFitsKey(visitfilepath, 'W_AITLWH')
        
    return data


def ImageQualityToCsv(butler, dataId, peaksList, csv_path=".",\
                      roi_size=16, EE=[3,5], seek_size=None,\
                      com=True, doBck=True, doFit=True, doLSF=False, doSep=False,fullSep=False,\
                      doPlot=False, doPrint=False, \
                      mask_size=50, threshold= 50, subpix = 5 , maxPeakDist=80,\
                      maxPeakFlux=40000, minPeakFlux=2000, experimentId=None, bkgId=None):
    """
    Calculate quality (EE, and sep info) for each peak from peaklist for a given visit defined in dataId dict
    butler is required to access the data
    use getImagequality and getImageEncerclEnergy(sep) functions
    imageInfo must a dict: dict(arm="b", spectrograph=2, visit=1234, filename="visit_file_path")
    csv file is <csv_path>\Imquality_{cam}_Exp{experimentId}_{visit}_{date_time}.csv
    Returns a pandas Dataframe and write csv file 
    
    """
    
    exp = butler.get("calexp", dataId)
    calexfilePath = butler.getUri("calexp", dataId)
    visit = dataId["visit"] 
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
    imageInfo.update(filename=calexfilePath)
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
    if bkgId is not None:
        bkg = butler.get("calexp", visit=bkgId, arm=dataId["arm"], spectrograph=dataId["spectrograph"])
        data = exp.image.array - bkg.image.array
    else:
        data = exp.image.array
    
    data = getFullImageQuality(data, peaksList, imageInfo=imageInfo,\
                      roi_size=roi_size, EE=EE, seek_size=seek_size,\
                      com=com, doBck=doBck, doFit=doFit, doLSF=doLSF, doSep=doSep,fullSep=fullSep,\
                      doPlot=doPlot, doPrint=doPrint, \
                      mask_size=mask_size, threshold= threshold, subpix = subpix , maxPeakDist=maxPeakDist,\
                      maxPeakFlux=maxPeakFlux, minPeakFlux=minPeakFlux, calexpMask=exp)
    
    now = datetime.now() # current date and time\n",
    date_time = now.strftime("%Y%m%dT%Hh%M")

    csvName = f"Imquality_{cam}_Exp{experimentId}_{visit}_{date_time}.csv"
    if not os.path.exists(csv_path):
        os.makedirs(csv_path,exist_ok =True)
    data.to_csv(os.path.join(csv_path, csvName))
    if doPrint:
        print(os.path.join(csv_path, csvName))

    return data


# this function is keep for compatibility but ImageQualityToCsv should be used instead

def VisitImageQualityToCsv(visit, \
                           peak_list, roi_size, com=True, doBck=True, EE=[3,5], doFit=True, doLSF=False,\
                           cam=None, repo="cluster",rerun="sm1-dither", cluster=False, doPlot=False,\
                           doSep=False, mask_size=50, threshold= 50, subpix = 5 , maxPeakDist=80,\
                           maxPeakFlux=40000, minPeakFlux=2000,\
                          csv_path = "", drpImage=None, experimentId = None, doPrint=False):
    """
    this function is keep for compatibility but ImageQualityToCsv should be used instead
    """
    if drpImage is None:
        experimentId = Logbook.visitExperimentId(visit=visit)
        
        doRaise = True
        try:
            filepath, date = getvisitfilepath(visit, rerun, None, cam, repo, cluster=cluster)
            visitfilepath = filepath
        except IOError:
            if doRaise:
                raise
    else :
        filepath = drpImage[0]
        visitfilepath = drpImage[1]
        visit = drpImage[2]
        try:
            #experimentId = Logbook.visitExperimentId(visit=visit)
            experimentId = get_Visit_Set_Id(visit)
        except:
            experimentId = np.nan
            print("Unable to get experimentId from logbook")

    data = getImageQuality(filepath, peak_list, roi_size=roi_size, EE=EE, com=com, doPlot=doPlot, doBck=doBck, doFit=doFit, doLSF=doLSF)
    if doSep:
        dsep = getImageEncerclEnergy(filepath, peak_list, roi_size=roi_size, EE=EE,\
        mask_size=mask_size, threshold= threshold, subpix = subpix ,\
        maxPeakDist=maxPeakDist, maxPeakFlux=maxPeakFlux, minPeakFlux=minPeakFlux,\
        doPlot=doPlot, doBck=doBck)
        dsep = dsep.add_prefix("sep_")
        dsep = dsep.rename(columns={'sep_peak': 'peak','sep_fiber': 'fiber'})
        data = data.merge(dsep, on=["peak","fiber"]) 
    __, fname = os.path.split(visitfilepath)

    data["filename"] = fname
    data["visit"] = visit

#    data["peaklist"] = peak_list        

#    data["motor1"] = np.float(getFitsKey(visitfilepath, 'W_XCU_MOTOR1_MICRONS'))
#    data["motor2"] = np.float(getFitsKey(visitfilepath, 'W_XCU_MOTOR2_MICRONS'))
#    data["motor3"] = np.float(getFitsKey(visitfilepath, 'W_XCU_MOTOR3_MICRONS'))
    data["xm1pos"] = np.float(getFitsKey(visitfilepath, 'W_XM1POS'))
    data["xm2pos"] = np.float(getFitsKey(visitfilepath, 'W_XM2POS'))
    data["xm3pos"] = np.float(getFitsKey(visitfilepath, 'W_XM3POS'))
    data["motor1"] = data["xm1pos"] 
    data["motor2"] = data["xm2pos"]
    data["motor3"] = data["xm3pos"]   
    fcax = np.float(getFitsKey(visitfilepath, 'W_ENFCAX', doRaise=False))
    fcay = np.float(getFitsKey(visitfilepath, 'W_ENFCAY', doRaise=False))
    fcaz = np.float(getFitsKey(visitfilepath, 'W_ENFCAZ', doRaise=False))

    # OneChannel back compatiblity
    fcax = np.float(getFitsKey(visitfilepath, 'HIERARCH W_FCA_FOCUS', doRaise=False)) if np.isnan(fcax) else fcax
    data['fcaFocus'] = fcax
    data['fcaY'] = fcay
    data['fcaZ'] = fcaz

#    feeTemp = np.float(getFitsKey(visitfilepath, 'temps.FEE'))
#    feeTemp = np.float(getFitsKey(visitfilepath, 'HIERARCH temps.FEE', doRaise=False)) if np.isnan(feeTemp) else feeTemp

#    data['feeTemp'] = feeTemp
#    ccdTemp = np.float(getFitsKey(visitfilepath, 'temps.CCD0'))
#    ccdTemp = np.float(getFitsKey(visitfilepath, 'HIERARCH temps.CCD0', doRaise=False)) if np.isnan(ccdTemp) else ccdTemp
#    data['ccdTemp'] = ccdTemp
    ccdTemp = np.float(getFitsKey(visitfilepath, 'W_XTDET1', doRaise=False))
    data['ccdTemp'] = ccdTemp
    detBoxTemp = np.float(getFitsKey(visitfilepath, 'W_XTDBOX', doRaise=False))
    data['detBoxTemp'] = detBoxTemp

    data['cam'] = cam
    data['obsdate'] = getFitsKey(visitfilepath, 'DATE-AVG')
    data['experimentId'] = experimentId

    now = datetime.now() # current date and time\n",
    date_time = now.strftime("%Y%m%dT%Hh%M")

    csvName = f"Imquality_{cam}_Exp{experimentId}_{visit}_{date_time}.csv"
    if not os.path.exists(csv_path):
        os.makedirs(csv_path)
    data.to_csv(csv_path+csvName)
    if doPrint:
        print(csv_path+csvName)

    return
