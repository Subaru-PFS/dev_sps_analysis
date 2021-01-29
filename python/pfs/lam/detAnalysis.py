from pfs.imageAnalysis import *
import scipy.interpolate
import matplotlib.gridspec as gridspec
from scipy import spatial
from datetime import datetime
from pfs.fileHandling import *

from pfs.sep import *


def getImageQuality(image, peak_list, roi_size=20, EE=[3,5], com=False, doPlot=False, scalePlot=False, doBck=False, doFit=True, doLSF=False):

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
            obj = getPeakData(image, cx,cy, EE=EE, roi_size=roi_size, com=com, doBck=doBck, doFit=doFit, doLSF=doLSF)
            obj["peak"] = row["peak"]
            obj["fiber"] = row["fiber"]
            obj["wavelength"] = wave
            obj["element"] = row["element"] if "element" in plist.columns else np.nan
            objlist.append(obj)
        except Exception as e:
            print(str(e), "cx:%i, cy:%i"%(cx,cy))
            objlist.append(dict(peak=row["peak"], fiber=row["fiber"], wavelength=wave))

    mdata = pd.DataFrame(objlist)
    if doPlot :
        plt_data = mdata[["peak", "fiber", "px", "py"]]
        plt_data = plt_data.rename(columns={'px': 'X','py': 'Y'})

        plotRoiPeak(image, plt_data, roi_size, scale=scalePlot)


    return mdata

def VisitImageQualityToCsv(visit, \
                           peak_list, roi_size, com=True, doBck=True, EE=[3,5], doFit=True, doLSF=False,\
                           cam=None, repo="cluster",rerun="sm1-dither", cluster=False, doPlot=False,\
                           doSep=False, mask_size=50, threshold= 50, subpix = 5 , maxPeakDist=80,\
                           maxPeakFlux=40000, minPeakFlux=2000,\
                          csv_path = "", drpImage=None, experimentId = None, doPrint=False):
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


# get Encercled Energy using SEP

def getImageEncerclEnergy(image, peak_list, roi_size=20, EE=[3,5], mask_size=50, threshold= 50, subpix = 5 , maxPeakDist=80, maxPeakFlux=40000, minPeakFlux=2000, doPlot=False, scalePlot=False, doBck=False, noEE=False):

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
                obj = getPeakDataSep(image, cx,cy, EE=EE, roi_size=roi_size,mask_size=mask_size,\
                                     subpix=subpix, doBck=doBck)
                obj["cx"] = cx
                obj["cy"] = cy
                obj["peak"] = row["peak"]
                obj["fiber"] = row["fiber"]              
                objlist.append(obj)
            except Exception as e:
                print(str(e), "cx:%i, cy:%i"%(cx,cy))
                objlist.append(dict(peak=row["peak"], fiber=row["fiber"]))

        mdata = pd.concat(objlist)
    else : # do it on the whole image so every peak
        obj = sep.extract(image, threshold)
        df = pd.DataFrame(obj, columns=obj.dtype.names)
        df = df[["flux", "peak", "x", "y", "flag", "npix"]]
        df = df.rename(columns={'x': 'px','y': 'py', 'peak': 'brightness'})
        df = removeClosePeak(df, dist=maxPeakDist, doPlot=doPlot)
        df = removeFluxPeak(df, fmax=maxPeakFlux, fmin=minPeakFlux, doPlot=doPlot)
        if noEE != True:
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
        plt_data = mdata[["peak", "fiber", "px", "py"]]
        plt_data = plt_data.rename(columns={'px': 'X','py': 'Y'})

        plotRoiPeak(image, plt_data, roi_size, scale=scalePlot)


    return mdata



def removeClosePeak(df, dist=80, doPlot=False):
    #from scipy import spatial
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_point.html#scipy.spatial.KDTree.query_ball_point
    
    points = np.c_[df.px.ravel(), df.py.ravel()]
    tree = spatial.KDTree(points)
    
    badpoints = []
    for point in points:
        result = tree.query_ball_point(point, dist)
        # remove the point itself 
        if len(result) > 1 : badpoints.append(result)
    if doPlot:
        points = np.asarray(points)
        plt.plot(points[:,0], points[:,1], '.')
        for badpoint in badpoints:
            nearby_point = points[badpoint]
            plt.plot(nearby_point[:,0], nearby_point[:,1], 'x')
        plt.show()
    
    # set a flag "ignore" in the dataframe to identify the peak to be ignored

    df["ignore"] = 0
    for result in badpoints:
        nearby_points = points[result]
        for nearby_point in nearby_points:
            df.loc[(df.px == nearby_point[0]) & (df.py == nearby_point[1]), "ignore"] = 1

#    return df.where(df.ignore<1).dropna()
    return df


# def removeBrightPeak(df, threshold=45000, doPlot=False):
#    if doPlot:
#        ax= df.plot.scatter(x="px", y="py")
#        df.where(df.brightness> threshold).plot.scatter(x="px", y="py", c="red", ax=ax)
#    
#    return df.where(df.brightness< threshold).dropna()

def removeFluxPeak(df, fmin=2000, fmax=40000, doPlot=False):
    if doPlot:
        ax= df.plot.scatter(x="px", y="py")
        df.where(df.brightness> fmax).plot.scatter(x="px", y="py", c="red", ax=ax)
    
    df.loc[(df.brightness< fmax) & (df.brightness> fmin), "ignore"] = 1

#    return df.where((df.brightness< fmax) & (df.brightness> fmin)).dropna()
    return df



def plotRoiPeak(image, peak_list, roi_size=20, raw=False, scale=True):
    if type(image) is str:     
        hdulist = fits.open(image, "readonly")
        image = hdulist[1].data
    
    plist = pd.read_csv(peak_list) if type(peak_list) is str else peak_list

    plist = plist.sort_values(["X","Y"], ascending=[True,False])    
    #peak_data = peak_data.reset_index()

    nbfiber = len(plist.fiber.unique())
    nbpeak = len(plist.peak.unique())
    # have a list of fiber from 0 to nbfiber 
    listfibertoplot = pd.DataFrame(plist.fiber.unique(), columns=["fiber"])
    
    print("#Fiber= %d and #peak= %d"%(nbfiber, nbpeak))
    f, axarr = plt.subplots(nbpeak, nbfiber,  sharex='col', sharey='row',figsize=(12,8))
#    print(axarr.shape)
    vmin=None
    vmax=None
    
    for (peak, fiber), group in plist.groupby(['peak','fiber']):
        k = int(round(peak))
        i = listfibertoplot[listfibertoplot.fiber == fiber].index.tolist()[0]
#        print(k,i)
        indy = group["X"]
        indx = group["Y"]
        cut_data = image[int(indx-roi_size/2):int(indx+roi_size/2), int(indy-roi_size/2):int(indy+roi_size/2)]
        if nbpeak == 1 and nbfiber == 1:
            axarr.imshow(cut_data,interpolation="none", origin="lower", vmin=vmin, vmax=vmax)
        else:
            axarr[nbpeak -1 -k, nbfiber - i -1].imshow(cut_data,interpolation="none", origin="lower", vmin=vmin, vmax=vmax)
    f.subplots_adjust(hspace=0.5,wspace=0.5)
    plt.show()


def getStatIM(dframe, par="EE5", thresold=0.9):
    z = dframe[par]
    zs = dframe[dframe[par]>thresold][par]
    print("%.f %% %s peaks >%s"%(100*len(zs)/len(z),par, thresold))

def plotImageQuality(dframe, vmin=-1,vmax=-1, par="EE3", hist=None, filelist=None, com=False, doSave=False, imgPath="/home/pfs/shared/Pictures/" ):
    
    # select peak center 
    # default x , y are objx and objy, but if it is the center of Mass it is oid_x and oid_y
#    x = dframe["objx"]
#    y = dframe["objy"]
#    if com :
#        x = dframe["oid_x"]
#        y = dframe["oid_y"]

    ## should now be px, py which are affected during calculation according com or not
    x = dframe["px"]
    y = dframe["py"]    
    
    z = dframe[par]
#    stat = 'ECE5' if par == 'ECE5' else 'EE5'
    stat = par
    xs = dframe[dframe[stat]>0.90].px
    ys = dframe[dframe[stat]>0.90].py
    zs = dframe[dframe[stat]>0.90][stat]

    print("%.f %% %s peaks >0.90"%(100*len(zs)/len(z),stat))
    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    # Interpolate
    rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    zi = rbf(xi, yi)

    if vmin == -1 :
        vmin = z.min()
    if vmax == -1 :
        vmax = z.max()
    if hist is not None:
        #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
        fig = plt.figure(figsize=(16, 8))
        gs = gridspec.GridSpec(1, 3,
 #                      width_ratios=[3,1],
 #                      height_ratios=[1,1]
                       )
        ax1 = plt.subplot(gs[0,:2])
        ax2 = plt.subplot(gs[0,2])
        im = ax1.imshow(zi, vmin=vmin, vmax=vmax, origin='lower',
                   extent=[x.min(), x.max(), y.min(), y.max()])
        ax1.scatter(x, y, c=z, vmin=vmin, vmax=vmax)
        ax1.set_title(par)
        ax1.set_xlim([0,4095])
        ax1.set_ylim([0,4175])
        fig.colorbar(im, ax=ax1, shrink=1)
        dframe[hist].plot.hist(ax=ax2, bins=20)
        ax2.set_title("Histogram")
        if filelist is not(None):
            fig.suptitle(getFileInfo(filelist))
        
    else : 
        fig = plt.figure(figsize=(8, 8))

        plt.imshow(zi, vmin=vmin, vmax=vmax, origin='lower',
                   extent=[x.min(), x.max(), y.min(), y.max()])
        plt.scatter(x, y, c=z, vmin=vmin, vmax=vmax)
        
        plt.colorbar(shrink=0.65)
        plt.scatter(xs, ys, c="r", marker="x")

        plt.xlim(0,4095)
        plt.ylim(0,4175)
        if filelist is not(None):
            plt.title(getFileInfo(filelist))
          
    if doSave:
        fig.patch.set_alpha(0.5)
        plt.savefig(os.path.join(imgPath,\
            f"SM1_ImQual_{filelist.cam[0]}_EXP{filelist.experimentId[0]}_{par}_{filelist.obsdate[0]}.png"))  
        plt.show()


def plotImageQualityScatter(dframe, par="EE3", vmin=-1,vmax=-1, hist=None, savePlotFile=None, com=False, doSave=False, title=None ):
    # select peak center 
    # default x , y are objx and objy, but if it is the center of Mass it is oid_x and oid_y
#    x = dframe["objx"]
#    y = dframe["objy"]
#    if com :
#        x = dframe["oid_x"]
#        y = dframe["oid_y"]

    ## should now be px, py which are affected during calculation according com or not
    x = dframe["px"]
    y = dframe["py"]   
    
    z = dframe[par]
    
    #    stat = 'ECE5' if par == 'ECE5' else 'EE5'
    val = 0.5 if par == "EE3" else 0.9

    stat = par
    xs = dframe[dframe[stat]>val].px
    ys = dframe[dframe[stat]>val].py
    zs = dframe[dframe[stat]>val][stat]

    statEE = f"{100*len(zs)/len(z):.1f}% of peak have a {par} > {val}"
        
    if vmin == -1 :
        vmin = z.min()
    if vmax == -1 :
        vmax = z.max()
    fact = 100
    if par == "brightness":
        fact = 1./100
    
    if hist is not None:
        #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
        fig = plt.figure(figsize=(16, 8))
        gs = gridspec.GridSpec(1, 3,
    #                      width_ratios=[3,1],
    #                      height_ratios=[1,1]
                       )
        ax1 = plt.subplot(gs[0,:2])
        ax2 = plt.subplot(gs[0,2])
        im = ax1.scatter(x, y, c=z, s= z*fact, vmin=vmin, vmax=vmax)
        ax1.set_title(par)
        ax1.set_xlim([0,4095])
        ax1.set_ylim([0,4175])
        plt.colorbar(im,ax=ax1,shrink=1)

        dframe[hist].plot.hist(ax=ax2, bins=20)
        val = 0.5 if par == "EE3" else 0.9
        ax2.axvline(x=val, color='k')
        ax2.set_title("Histogram")
        ax2.set_xlim(vmin,vmax)
        if title is not(None):
            fig.suptitle(title+"\n"+statEE)
    else:
        fig = plt.figure(figsize=(10, 8))

        plt.scatter(x, y, c= z, s= z*fact, vmin=vmin, vmax=vmax)
        plt.colorbar(shrink=1)
        plt.xlim(0,4095)
        plt.ylim(0,4175)
        plt.xlabel('X')
        plt.ylabel('Y')

        plt.show()
        if title is not(None):
            fig.suptitle(title+"\n"+statEE)
          
    if doSave:
        fig.patch.set_alpha(0.5)
        plt.savefig(savePlotFile+f"_{par}.png", bbox_inches = "tight" )
        plt.show()


def getStatIM(dframe, par="EE5", thresold=0.9):
    z = dframe[par]
    zs = dframe[dframe[par]>thresold][par]
    print("%.f %% %s peaks >%s"%(100*len(zs)/len(z),par, thresold))

