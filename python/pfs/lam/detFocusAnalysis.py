# -*- coding: utf-8 -*-
from pfs.detAnalysis import *
from pfs.fileHandling import *

from pfs.sep import *

def getAllImageQuality(filelist, peak_list, roi_size, com=False, doBck=False, doPlot=False, EE=[3,5], trackPeak=False, doFit=True, doLSF=False):
   
    imdata = []
    for cam, camfilelist in filelist.groupby('cam'):
        for j, row in camfilelist.iterrows():
        # keep only filename
            data = getImageQuality(row.filepath, peak_list, roi_size=roi_size, EE=EE, com=com, doPlot=doPlot, doBck=doBck, doFit=doFit, doLSF=doLSF)
            __, fname = os.path.split(row.filepath)

            data["filename"] = fname        
            data["motor1"] = np.float(getFitsKey(row.filepath, 'W_XCU_MOTOR1_MICRONS'))
            data["motor2"] = np.float(getFitsKey(row.filepath, 'W_XCU_MOTOR2_MICRONS'))
            data["motor3"] = np.float(getFitsKey(row.filepath, 'W_XCU_MOTOR3_MICRONS'))
            fcax = np.float(getFitsKey(row.filepath, 'W_ENFCAX', doRaise=False))
            fcay = np.float(getFitsKey(row.filepath, 'W_ENFCAY', doRaise=False))
            fcaz = np.float(getFitsKey(row.filepath, 'W_ENFCAZ', doRaise=False))

            # OneChannel back compatiblity
            fcax = np.float(getFitsKey(row.filepath, 'HIERARCH W_FCA_FOCUS', doRaise=False)) if np.isnan(fcax) else fcax
            data['fcaFocus'] = fcax
            data['fcaY'] = fcay
            data['fcaZ'] = fcaz

            feeTemp = np.float(getFitsKey(row.filepath, 'temps.FEE'))
            feeTemp = np.float(getFitsKey(row.filepath, 'HIERARCH temps.FEE', doRaise=False)) if np.isnan(feeTemp) else feeTemp
            
            data['feeTemp'] = feeTemp
            ccdTemp = np.float(getFitsKey(visitfilepath, 'W_XTDET1', doRaise=False))
            data['ccdTemp'] = ccdTemp
            detBoxTemp = np.float(getFitsKey(visitfilepath, 'W_XTDBOX', doRaise=False))
            data['detBoxTemp'] = detBoxTemp
            #ccdTemp = np.float(getFitsKey(row.filepath, 'temps.CCD0'))
            #ccdTemp = np.float(getFitsKey(row.filepath, 'HIERARCH temps.CCD0', doRaise=False)) if np.isnan(ccdTemp) else ccdTemp
            data['cam'] = cam
            data['obsdate'] = getFitsKey(row.filepath, 'DATE-AVG')
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



# get best focus using fit gauss 1D
# Return the the best focus found and the corresponding peak centroid 

def getFocus(series, criteria="EE5", index='relPos', doPrint=False, doRaise=False):
    
    try:
        if ('fwhm' in criteria) or ('2ndM' in criteria):
            thfoc = fitparabola(series[index].values, series[criteria].values, deg=3)  
        else:
            thfoc, parfit = fitgauss1D(series[index].values, series[criteria].values)
    except RuntimeError:
        if doRaise:
            raise
        else:
            print('could not fit data for %s'%criteria)
            thfoc = interpdata(series[index].values, series[criteria].values, criteria)

    if doPrint:
        print("Best Focus(%s) = %.2f µm  %.2f"%(criteria, thfoc.focus, thfoc.y.max()))
    
    return thfoc.rename(index=str, columns={"x": index, "y": criteria})

def fitFocusData(imdata, index='relPos', criterias=['EE5','EE3', 'brightness', 'fwhm'], doPlot=False, doPrint=False, head=0, tail=0):
    thfoc_data = []
    tmpdata = imdata[head:imdata.count()[0]-tail]

#    criterias = critierias if criterias is None else ['EE5', 'EE3', 'brightness', 'fwhm']

    for (wavelength, fiber), series in tmpdata.groupby(['wavelength','fiber']):
        #series = series.dropna()
        thfoc = getFocus(series, criteria=criterias[0], index=index, doPrint=doPrint)
        for criteria in criterias:
            thfoc[criteria] = getFocus(series, criteria=criteria, index=index, doPrint=doPrint)[criteria]

        thfoc['peak'] =  series.peak.unique()[0]
        thfoc['wavelength'] = wavelength
        thfoc['fiber'] = fiber
        thfoc['px'] = np.interp(thfoc[index], series[index], series['px'])
        thfoc['py'] = np.interp(thfoc[index], series[index], series['py'])
        # re-create motor value with the higer sampling given by gaussfit
        thfoc['motor1'] = thfoc[index] #+ series.motor1.min()
        thfoc['motor2'] = thfoc[index] #+ series.motor2.min()
        thfoc['motor3'] = thfoc[index] #+ series.motor3.min()

        thfoc_data.append(thfoc)

    thfoc_data = pd.concat(thfoc_data)
    
    if doPlot:
        kwargs = dict(grid=True, figsize=(14,10), legend=True, subplots=True)
        criterias = ['EE5','EE3', 'brightness', 'fwhm']
        for (wavelength, fiber), fit in thfoc_data.groupby(['wavelength','fiber']):
### WARINIG: CA VA PAS MARCHER !! ??
            imdata['waveStr'] = imdata['wavelength'].map('{:,.5f}'.format).astype("str")
            raw = imdata.query("waveStr==%d and fiber==%d"%(wavelength, fiber))
            axes = fit.set_index('motor1')[criterias].plot(**kwargs)
            for i, criteria in enumerate(criterias):
                axes[i].plot(raw['motor1'].values, raw[criteria].values, 'o')
                
    return thfoc_data

def getFocusMap(fitdata, index='relPos', criterias=['EE5','EE3', 'brightness', 'fwhm']): 
    data = []
    for (wavelength, fiber), series in fitdata.groupby(['wavelength','fiber']):
        #series = series.dropna()
        peak = series.peak.unique()[0]
        for criteria in criterias:
            if ('fwhm' in criteria) or ('2ndM' in criteria) :
                ixmax = series[criteria].idxmin()
                value = series[criteria].min()
            else:
                ixmax = series[criteria].idxmax()
                value = series[criteria].max()
            focus = series[index][ixmax]                
            px = series.px[ixmax]
            py = series.py[ixmax]
            mat = [peak,wavelength,fiber,criteria, px, py, focus, focus+series.motor1.min(), focus+series.motor2.min(), focus+series.motor3.min(), value]
            data.append(tuple(mat))
    
    columns = ['peak','wavelength', 'fiber', 'criteria', 'px', 'py', 'relPos', 'motor1', 'motor2', 'motor3', 'value']
    return pd.DataFrame(data, columns=columns)

def getBestPlane(data, order=1, doPlot=False, plot_path=None, exp=None):
    # data should contains x,y,value
    import scipy.linalg
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt    

    lowerBound_x = 0
    upperBound_x = 4200

    lowerBound_y = 0
    upperBound_y = 4200


    # regular grid covering the domain of the data
    X,Y = np.meshgrid(np.arange(lowerBound_x, upperBound_x, 100), np.arange(lowerBound_y, upperBound_y, 100))
    XX = X.flatten()
    YY = Y.flatten()

    # 1: linear, 2: quadratic
    if order == 1:
        # best-fit linear plane
        A = np.c_[data.px, data.py, np.ones(data.px.shape[0])]
        C,_,_,_ = scipy.linalg.lstsq(A, data.relPos)    # coefficients

        # evaluate it on grid
        Z = C[0]*X + C[1]*Y + C[2]

        # or expressed using matrix/vector product
        #Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)

    elif order == 2:
        # best-fit quadratic curve
#        A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
#        C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
        A = np.c_[np.ones(data.shape[0]), data[["px","py"]], \
                  np.prod(data[["px","py"]].values, axis=1),\
                  data[["px","py"]].values**2]
        C,_,_,_ = scipy.linalg.lstsq(A, data.relPos)
        # evaluate it on a grid
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)

    if doPlot:
    # plot points and fitted surface
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
        ax.scatter(data.px, data.py, data.relPos, c='r', s=50)
        plt.xlabel('X')
        plt.ylabel('Y')
        ax.set_zlabel('Z')
#        ax.axis('equal')
#        ax.axis('tight')

        plt.suptitle(f"Exp{exp}")

        plt.savefig(plot_path+f"Focus_plane_Exp{exp}.png")
        plt.show()
    
    return C


def findMotorPos(plane, inv_mat=None, doPrint=False):
    if inv_mat is not None:
        inv_mat = np.load(inv_mat, allow_pickle=True)
    if doPrint:
        print(np.dot(inv_mat, plane))
    return np.dot(inv_mat, plane)

def findInvMatrix(dfPlanes, tilt0=None, doPrint=False, saveMat=None):
    calTilt = round(dfPlanes[["motor1", "motor2", "motor3"]].max().iloc[0]-dfPlanes[["motor1", "motor2", "motor3"]].min().iloc[0])
    tilt0 = calTilt if tilt0 is None else tilt0
    planes = dfPlanes[["a","b","c"]].values
    
    mat = np.zeros((3,3))
    mat[:,0] = planes[0,:] - planes[1,:] 
    mat[:,1] = planes[0,:] - planes[2,:]
    mat[:,2] = planes[0,:] - planes[3,:]

    mat = mat /tilt0
    print(mat)
    inv_mat = np.linalg.inv(mat)
    print(inv_mat)
    if doPrint:
        print(f'tilt: {tilt0}')
        print(np.dot(inv_mat,mat))
        print(f'plane {np.dot(inv_mat, planes[0,:])}')
    if saveMat is not None:
        inv_mat.dump(saveMat)
    return inv_mat
