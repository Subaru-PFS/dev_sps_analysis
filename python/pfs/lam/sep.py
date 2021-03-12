# Functions that use sep

import sep

def createSepMask(shape, cx, cy, mask_size=50):
    indx = cy
    indy = cx
    mask = np.ones(shape)
    xmin = int(indx-mask_size/2) if int(indx-mask_size/2) >0 else 0
    xmax = int(indx+mask_size/2) if int(indx+mask_size/2) < shape[1] else shape[1]
    ymin = int(indy-mask_size/2) if int(indy-mask_size/2) >0 else 0
    ymax = int(indy+mask_size/2) if int(indy+mask_size/2) < shape[0] else shape[0]
    
    mask[xmin:xmax, ymin: ymax]= 0.
    
    return mask

def getPeakDataSep2(image, cx, cy, EE=None, roi_size=30, mask_size=50, threshold= 500, subpix=5, \
                            doPlot=False, doBck=True, nRows=5, nRound=0, **kwargs):
    if type(image) is str:
        hdulist = fits.open(image, "readonly")
        image = hdulist[1].data
    
    objlist=[]
    
    outer_data, inner_data = getRois(image, cx, cy, inner_size=5, outer_size=roi_size, doBck=doBck, nRows=nRows, nRound=nRound)
    outer_data = np.ascontiguousarray(outer_data, dtype=np.float32)


#    mask = createSepMask(image.shape, cx, cy, mask_size=mask_size)
    obj = sep.extract(outer_data, threshold)
    if len(obj) == 0 :
        df = pd.DataFrame({'cx': [cx], 'cy': [cy]})
        return df
    
    df = pd.DataFrame(obj, columns=obj.dtype.names)
    df = df[["flux", "peak", "x", "y", "flag", "npix", "theta"]]
    df = df.rename(columns={'x': 'px','y': 'py', 'peak': 'brightness'})
    
    flux_tot, flux_tot_err, flux_tot_flag = sep.sum_circle(outer_data, df['px'], df['py'],
                                     roi_size/2., err=None, gain=1.0, subpix=subpix)
    df["flux_tot"] = flux_tot

    EE = [3,5] if EE is None else EE
    for ee in EE:
        df["ECE%d"%ee], df["fluxErr"], flag = sep.sum_circle(outer_data, df['px'], df['py'],
                                     ee/2., err=None, gain=1.0, subpix=subpix)
        df["ECE%d"%ee] = df["ECE%d"%ee] / flux_tot
    
    df.px = df.px + roi_size/2
    df.py = df.py + roi_size/2
    
    return df



def getPeakDataSep(image, cx, cy, EE=None, roi_size=30, mask_size=50,  seek_size=None, threshold= 500, subpix=5, \
                            doPlot=False, doBck=True, nRows=5, **kwargs):
    if type(image) is str:
        hdulist = fits.open(image, "readonly")
        image = hdulist[1].data
    
    objlist=[]
    if seek_size is not None:
        (cx,cy) = estimateCOM(image, cx, cy, roi_size=seek_size, doBck=doBck, nRows=nRows)     
        
    mask = createSepMask(image.shape, cx, cy, mask_size=mask_size)
    obj = sep.extract(image, threshold, mask=mask)
    if len(obj) == 0 :
        df = pd.DataFrame({'cx': [cx], 'cy': [cy]})
        return df
    
    df = pd.DataFrame(obj, columns=obj.dtype.names)
#    df = df[["flux", "peak", "x", "y", "flag", "npix", "theta"]]
    df = df.rename(columns={'x': 'px','y': 'py', 'peak': 'brightness'})
    
    if doBck:
        flux_tot, flux_tot_err, flux_tot_flag = sep.sum_circle(image, df['px'], df['py'],
                                         roi_size/2., err=None, gain=1.0, subpix=subpix,bkgann=(roi_size/2., nRows+roi_size/2.))
        df["flux_tot"] = flux_tot

        EE = [3,5] if EE is None else EE
        for ee in EE:
            df["flux_EC%d"%ee], df["fluxErr"], flag = sep.sum_circle(image, df['px'], df['py'],
                                                         ee/2., err=None, gain=1.0, subpix=subpix, 
                                                         bkgann=(roi_size/2., nRows+roi_size/2.))
            df["flux_CR%d"%ee], df["fluxErrR%d"%ee], flag = sep.sum_circle(image, df['px'], df['py'],
                                                         np.sqrt(2)*ee/2., err=None, gain=1.0, subpix=subpix,
                                                         bkgann=(roi_size/2., nRows+roi_size/2.))  
            df["ECE%d"%ee] = df["flux_EC%d"%ee] / flux_tot
            df["ECR%d"%ee] = df["flux_CR%d"%ee] / flux_tot

    else:
        flux_tot, flux_tot_err, flux_tot_flag = sep.sum_circle(image, df['px'], df['py'],
                                         roi_size/2., err=None, gain=1.0, subpix=subpix)
        df["flux_tot"] = flux_tot

        EE = [3,5] if EE is None else EE
        for ee in EE:
            df["flux_EC%d"%ee], df["fluxErr"], flag = sep.sum_circle(image, df['px'], df['py'],
                                         ee/2., err=None, gain=1.0, subpix=subpix)
            df["ECE%d"%ee] = df["flux_EC%d"%ee] / flux_tot
            df["flux_CR%d"%ee], df["fluxErrR%d"%ee], flag = sep.sum_circle(image, df['px'], df['py'],
                                                         np.sqrt(2)*ee/2., err=None, gain=1.0, 
                                                        subpix=subpix)  
            df["ECE%d"%ee] = df["flux_EC%d"%ee] / flux_tot
            df["ECR%d"%ee] = df["flux_CR%d"%ee] / flux_tot
         
    return df

# only calculate flux
def getFluxPeakDataSep(image, cx, cy, EE=None, roi_size=30, mask_size=50, threshold= 50, subpix=5, \
                            doPlot=False, doBck=True, nRows=5, **kwargs):
    if type(image) is str:

        hdulist = fits.open(image, "readonly")
        image = hdulist[1].data
    
    objlist=[]

#    mask = createSepMask(image.shape, cx, cy, mask_size=mask_size)
    outer_data, inner_data = getRois(image, cx, cy, inner_size=5, outer_size=50, doBck=doBck)
    data = outer_data.copy(order='C')
    
    obj = sep.extract(data, threshold ) #, mask=mask)
    if len(obj) == 0 :
        df = pd.DataFrame({'cx': [cx], 'cy': [cy]})
        return df
    
    df = pd.DataFrame(obj, columns=obj.dtype.names)
#    df = df[["flux", "peak", "x", "y", "flag", "npix", "theta"]]
    df = df.rename(columns={'x': 'px','y': 'py', 'peak': 'brightness'})

    
    
    flux_tot, flux_tot_err, flux_tot_flag = sep.sum_circle(data, df['px'], df['py'],
                                         roi_size/2., err=None, gain=1.0, subpix=subpix)
    df["flux_tot"] = flux_tot

    EE = [3,5] if EE is None else EE
    for ee in EE:
        df["EC%d"%ee], df["fluxErr%d"%ee], flag = sep.sum_circle(data, df['px'], df['py'],
                                     ee/2., err=None, gain=1.0, subpix=subpix)
        df["ECR%d"%ee], df["fluxErrR%d"%ee], flag = sep.sum_circle(data, df['px'], df['py'],
                                     np.sqrt(2)*ee/2., err=None, gain=1.0, subpix=subpix)         
    indx = cy
    indy = cx
    df.px = int(indy-roi_size/2) + df.px + 0.5
    df.py = int(indx-roi_size/2) + df.py + 0.5
    
    if doPlot:
        # plot background-subtracted image
        fig, ax = plt.subplots()
        m, s = np.mean(data), np.std(data)
        s = s
        im = ax.imshow(data, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')
    
    return df


def getAllSepFlux(filelist, peak_list, roi_size, com=False, doBck=False, doPlot=False, EE=[3,5], trackPeak=False):
   
    imdata = []
    for cam, camfilelist in filelist.groupby('cam'):
        for j, row in camfilelist.iterrows():
        # keep only filename
            data = getImageEncerclEnergy(row.filepath, peak_list, roi_size=roi_size, EE=EE, doBck=doBck)
            __, fname = os.path.split(row.filepath)

            data["filename"] = fname        
            data["motor1"] = np.float(getFitsKey(row.filepath, 'W_XCU_MOTOR1_MICRONS'))
            data["motor2"] = np.float(getFitsKey(row.filepath, 'W_XCU_MOTOR2_MICRONS'))
            data["motor3"] = np.float(getFitsKey(row.filepath, 'W_XCU_MOTOR3_MICRONS'))
            fcax = np.float(getFitsKey(row.filepath, 'W_ENFCAX', doRaise=False))
            # OneChannel back compatiblity
            fcax = np.float(getFitsKey(row.filepath, 'HIERARCH W_FCA_FOCUS', doRaise=False)) if np.isnan(fcax) else fcax
            data['fcaFocus'] = fcax
            feeTemp = np.float(getFitsKey(row.filepath, 'temps.FEE'))
            feeTemp = np.float(getFitsKey(row.filepath, 'HIERARCH temps.FEE', doRaise=False)) if np.isnan(feeTemp) else feeTemp
            
            data['feeTemp'] = feeTemp
            ccdTemp = np.float(getFitsKey(row.filepath, 'temps.CCD0'))
            ccdTemp = np.float(getFitsKey(row.filepath, 'HIERARCH temps.CCD0', doRaise=False)) if np.isnan(ccdTemp) else ccdTemp
            data['ccdTemp'] = ccdTemp
            data['cam'] = cam
            data['obsdate'] = row.obsdate
            data['experimentId'] = row.experimentId
            imdata.append(data)
            # to udpate peak centroid using the previous image
            if trackPeak == True:
                plist = pd.read_csv(peak_list) if type(peak_list) is str else peak_list
                plist.X[0] = data.px
                plist.Y[0] = data.py
                peak_list = plist
        
    imdata = pd.concat(imdata, ignore_index=True)
    minPos = imdata[['motor1','motor2','motor3']]
    minPos = minPos - minPos.min()
    imdata['relPos'] = minPos['motor1']
    
    return imdata

