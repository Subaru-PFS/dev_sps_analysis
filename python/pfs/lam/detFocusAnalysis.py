# -*- coding: utf-8 -*-
from pfs.lam.detAnalysis import *
from pfs.lam.fileHandling import *
from pfs.lam.imageAnalysis import parabola

from pfs.lam.sep import *


# Extending Pandas Dataframe for thFocus fit data

@pd.api.extensions.register_dataframe_accessor("thFocus")
class thFocusAccessor:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj
        
    @staticmethod
    def _validate(obj):
        # verify there is a column latitude and a column longitude
        if "focus" not in obj.columns or "value" not in obj.columns:
            raise AttributeError("Must have 'focus' and 'value'.")
    @property
    def vline(self):
        return dict(x=self._obj.focus, ymin=0, ymax=self._obj.value)

    @property
    def fitdata(self):
        fit_function = self._obj.fit_function.values[0]
        xmin = self._obj.fit_xmin
        xmax = self._obj.fit_xmax
        newx = np.linspace(xmin, xmax, 100)
        best = self._obj.focus.values[0]
        value = self._obj.value.values[0]
        width = self._obj.width.values[0]
        
        if fit_function == "thFocus_Gaussian":
            return newx, thFocus_Gaussian(newx, best, value, width)
        elif fit_function == "thFocus_Parabola":
            return newx, thFocus_Parabola(newx, best, value, width)

# Class to create DataFrame that contains thfocus fit data
        
class ThFocusDF(pd.DataFrame):
    def __init__(self, focus, value, width, fit_function, fit_xmin, fit_xmax):
        pd.DataFrame.__init__(self, data=[(focus, value, width, fit_function, fit_xmin, fit_xmax)],\
                              columns=["focus", "value", "width", "fit_function", "fit_xmin", "fit_xmax"])
        
def thFocus_Gaussian(motor_piston, motor_peak, ee_peak, width):
    """
    Gaussian curve used for Through Focus analysis
    """
    return ee_peak * np.exp(-np.power( (motor_piston -  motor_peak) / width, 2.))


def fit_thFocus_Gaussian(x, y, width=None, max_value=None):
    """
    Run the fit of the thFocus_Gaussian
    Estimate initial parameter
    if width and/or max_value is given they are pass as a fixed value and thus not part of optimisation parameters
    run curve_fit
    return a DataFrame as defined by ThFocusDF class
    """
    ee0 = np.max(y)
    motor0 = x[np.argmax(y)]
    
 
    # width
    half_pos = x[(np.abs(y - ee0/2)).argmin()]
    width0 = 2 * np.abs((motor0 - half_pos))
    
    if width is None and max_value is None:
        [focus_fit, value_fit, width_fit], pcov = curve_fit(thFocus_Gaussian, x, y, p0=[motor0,ee0, width0], maxfev=10000, sigma=None, absolute_sigma= True)
    elif width is None and max_value is not(None):
        value_fit = max_value
        [focus_fit, width_fit], pcov = curve_fit(lambda x, best, width: thFocus_Gaussian(x, best, max_value, width), x, y, p0=[motor0, width0], maxfev=10000, sigma=None, absolute_sigma= True) 
    elif width is not(None) and max_value is None:
        width_fit = width
        [focus_fit, value_fit], pcov = curve_fit(lambda x, best, value: thFocus_Gaussian(x, best, value, width), x, y, p0=[motor0,ee0], maxfev=10000, sigma=None, absolute_sigma= True)
    else:
        width_fit = width
        value_fit = max_value
        [focus_fit], pcov = curve_fit(lambda x, best: thFocus_Gaussian(x, best, max_value, width), x, y, p0=[motor0], maxfev=10000, sigma=None, absolute_sigma= True)

    
    return ThFocusDF(focus_fit, value_fit, width_fit, "thFocus_Gaussian", np.min(x), np.max(x))


def thFocus_Parabola(motor_piston, motor_peak, value_peak, width):
    """
    Parabola curve used for Through Focus analysis
    f(x) = (1/width) * (x - motor_peak) + EE_peak
    with motor_peak = -b/2a
    and value_peak = (4ac - b*b)/4a
    latus / width = 1/a
    """
    
    
    return (1/width) * np.power( (motor_piston - motor_peak), 2.) + value_peak

def fit_thFocus_parabola(x, y):
    """
    Run Parabola fit for Through Focus analysis
    
    f(x) = a(x - motor_peak) + EE_peak
    with motor_peak = -b/2a
    and EE_peak = (4ac - b*b)/4a
    latus / width = 1/a
    """
    [a, b, c], pcov = scipy.optimize.curve_fit(parabola, x, y)

  
    min_radius = (4*a*c-b*b)/(4*a)
    min_pos = -b/(2*a)
    latus = 1/a
    
    #return ThFocusDF([min_pos, min_radius, latus]), np.sqrt(np.diag(pcov))
    return ThFocusDF(min_pos, min_radius, latus, "thFocus_Parabola", np.min(x), np.max(x))

def thF_interpdata(x, y, criteria):
    f = interpolate.interp1d(x, y, kind='cubic')
    newx = np.linspace(np.min(x), np.max(x), 1000)
    
    if criteria == 'min':
            focus_fit = np.argmin(f(newx))
            value_fit = np.min(f(newx))
    else:
            focus_fit = np.argmax(f(newx))
            value_fit = np.max(f(newx))        
        
    return ThFocusDF(*[focus_fit, value_fit, None])


def getBestFocus(series, criteria="EE5", index='relPos', width=None, max_value=None, doPrint=False, doRaise=False):
    
    try:
        if ('fwhm' in criteria) or ('sep_2ndM' in criteria):
            thfoc = fit_thFocus_parabola(series[index].values, series[criteria].values)  
        else:
            thfoc = fit_thFocus_Gaussian(series[index].values, series[criteria].values, width=width, max_value=max_value)
    except RuntimeError:
        if doRaise:
            raise
        else:
            print('could not fit data for %s'%criteria)
            thfoc = ThFocusDF(*(np.nan * np.ones(6)))
            #thfoc = thF_interpdata(series[index].values, series[criteria].values, criteria)

    if doPrint:
        print("Best Focus(%s) = %.2f µm  %.2f"%(criteria, thfoc.focus, thfoc.value))
    
    return thfoc


def getAllBestFocus(piston, index="relPos", criterias=["EE5", "EE3", "2ndM"], doPlot=False, doPrint=False, head=0, tail=0, \
                   savePlot=False, plot_path="", plot_title=None, width=None, max_value=None):
    thfoc_data = []
    tmpdata = piston[head:piston.count()[0]-tail]
    for criteria in criterias:
        if criteria in piston.columns:
            for (wavelength, fiber), series in tmpdata.groupby(['wavelength','fiber']):
                thfoc = getBestFocus(series, criteria, index=index, width=width, max_value=max_value)
                thfoc['peak'] =  series.peak.unique()[0]
                thfoc['wavelength'] = wavelength
                thfoc['fiber'] = fiber
                thfoc['criteria'] = criteria
                thfoc['axis'] = index

                #thfoc['visit'] = series.visit.unique()[0]
                thfoc['experimentId'] = series.experimentId.unique()[0]


                thfoc['px'] = np.interp(thfoc.focus, series[index], series['px'])
                thfoc['py'] = np.interp(thfoc.focus, series[index], series['py'])
                thfoc_data.append(thfoc)
        else:
            if doPrint:
                print(f"{criteria} not in DataFrame")
    thfoc_data = pd.concat(thfoc_data, ignore_index=True)
    
    if doPlot:
        
        grouped = piston.groupby(['wavelength','fiber'])
        grouped_focus = thfoc_data.groupby(['wavelength','fiber'])
        # Need to work on a better criteria that match good focus-able data and others....
        width_limit = 2*thfoc_data.width.mean()

        ncols = len(piston.fiber.unique())
        nrows = len(piston.wavelength.unique())

        newx = np.linspace(np.min(piston[index].values), np.max(piston[index].values), 100)
        for criteria in criterias:
            if criteria in piston.columns:
                fig, axs = plt.subplots(nrows,ncols, figsize=(20,12), sharey=True, sharex=True)
                visit_info = f"{piston.visit.values.min()} to {piston.visit.values.max()}"
                temp_info = f"detBox: {piston.detBoxTemp.mean():.1f}K  ccd: {piston.ccdTemp.mean():.1f}K"
                cam_info = piston.cam.unique()[0]
                date_info = piston.obsdate[0].split('T')[0]
                fca_focus = piston.fcaFocus.mean()
                fig.suptitle(f"{cam_info.upper()} ExpId {str(int(piston.experimentId.unique()[0]))} - {visit_info} - {criteria} - {temp_info} - FCA Focus {fca_focus:.1f}mm - {date_info}")
                plt.subplots_adjust(top=0.93)
                for (name, df), ax, (f, focus) in zip(grouped, axs.flat, grouped_focus):
                    ax.set_title(f"{name[0]:.2f}, {name[1]}")
                    df.plot.scatter(x=index,y=criteria, ax=ax)
                    if np.isnan(focus.focus.values) != True :
                        ax.plot(*focus[focus.criteria == criteria].thFocus.fitdata, "r")
                        if focus.width.values < width_limit:
                            ax.vlines(**focus[focus.criteria == criteria].thFocus.vline)
        if savePlot:
            if plot_title is None:
                plot_title = f"{cam_info.upper()}_ExpId_{str(int(piston.experimentId.unique()[0]))}_{criteria}_thFocusPlot{date_info}.png"
            plt.savefig(plot_path+plot_title)
                    
    return thfoc_data
    
    
def fit3dPlane(df, coords=["x","y","z"], order=1, x_bound=None, y_bound=None, \
               doPlot=False, plot_path=None, plot_title=None, plot_name=None, savePlot=False):
#def fit3dPlane(df, coords=["x","y","z"], order=1, doPlot=False, plot_path=None, exp=None):
    
    import scipy.linalg
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    
    x = df[coords[0]]
    y = df[coords[1]]
    z = df[coords[2]]
    x_bound = [x.min(),x.max()] if x_bound is None else x_bound
    y_bound = [y.min(),y.max()] if y_bound is None else y_bound

    # regular grid covering the domain of the data
    X,Y = np.meshgrid(np.arange(x_bound[0], x_bound[1], 100), np.arange(y_bound[0], y_bound[1], 100))
    XX = X.flatten()
    YY = Y.flatten()

    # 1: linear, 2: quadratic
    if order == 1:
        # best-fit linear plane
        A = np.c_[x, y, np.ones(x.shape[0])]
        C,_,_,_ = scipy.linalg.lstsq(A, z)    # coefficients

        # evaluate it on grid
        Z = C[0]*X + C[1]*Y + C[2]

        # or expressed using matrix/vector product
        #Z = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)

    elif order == 2:
        # best-fit quadratic curve
        A = np.c_[np.ones(x.shape[0]), df[[coords[0],coords[1]]], \
                  np.prod(df[[coords[0],coords[1]]].values, axis=1),df[[coords[0],coords[1]]].values**2]
        C,_,_,_ = scipy.linalg.lstsq(A, z)
        # evaluate it on a grid
        Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)

    if doPlot:
    # plot points and fitted surface
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
        ax.scatter(x, y, z, c='r', s=50)
        plt.xlabel('X')
        plt.ylabel('Y')
        ax.set_zlabel('Z')
#        ax.axis('equal')
#        ax.axis('tight')

        if plot_title is not None:
            plt.suptitle(plot_title)

        if savePlot :
            if plot_path is None or plot_title is None:
                  raise Exception("you must specify a plot path and a title")
            plt.savefig(plot_path+plot_title)
        plt.show()
#        plt.cla()
    
    return C

def getBestPlane(data, order=1, doPlot=False, plot_path=None, exp=None, coords=["px", "py", "relPos"]):
    
    #coords = ["px", "py", "relPos"]
    x_bound = [0,4096]
    y_bound = [0,4176]
    
    if exp is not None:
        plot_title = f"Focus_plane_Exp{exp}.png"
        plot_name = f"Exp{exp}"
    else :
        plot_title = None
        plot_name = None

    return fit3dPlane(data, coords=coords, order=1, x_bound=x_bound, y_bound=y_bound, \
                      doPlot=doPlot, plot_path=plot_path, plot_title=plot_title, plot_name=plot_name)


def getFocusInvMat(cam):
    if cam is not None:
        arm = cam[0]
        inv_mat_path = os.path.join(os.environ['LAM_SPS_ANALYSIS_DIR'],"notebooks/optical/CamUnitAlignement/invMat")
        if arm == "r" or arm =="m":
            invMatName = "InvMat_sm1_R1_17sept2020.mat"
        elif arm == "b":
            invMatName = "InvMat_sm1_R1_17sept2020.mat"
        else:
            raise Exception("arm must be b, r or m")
    else:
        raise Exception("Either inv_mat or arm must be provided")
    inv_mat = os.path.join(inv_mat_path, invMatName)
    invMat = np.load(inv_mat, allow_pickle=True)
    return invMat, invMatName
    


def findMotorPos(plane, inv_mat=None, doPrint=False, cam=None):
    """
    Determine focus motor position given the best plane.
    it used a inversion matrix that gives each motor effect
    motor_pos  = np.dot(invMat, best plane)
    if no inv_mat is specify you have to specify the cam argument
    an inv_mat will be automatically used
    """
    if inv_mat is not None:
        invMat = np.load(inv_mat, allow_pickle=True)
    elif cam is not None:
        arm = cam[0]
        inv_mat_path = os.path.join(os.environ['LAM_SPS_ANALYSIS_DIR'],"notebooks/optical/CamUnitAlignement/invMat")
        if arm == "r" or arm =="m":
            inv_mat = os.path.join(inv_mat_path, "InvMat_sm1_R1_17sept2020.mat")
            invMat = np.load(inv_mat, allow_pickle=True)
        elif arm == "b":
            inv_mat = os.path.join(inv_mat_path, "InvMat_sm1_B1_02oct2020.mat")
            invMat = np.load(inv_mat, allow_pickle=True)
        else:
            raise Exception("arm must be b, r or m")
    else:
        raise Exception("Either inv_mat or arm must be provided")
    
    motors_pos = np.dot(invMat, plane)
    if doPrint:
        print(inv_mat)
        print("xcu_%s motors moveCcd a=%.2f b=%.2f c=%.2f microns abs"%(cam,*motors_pos))
    return motors_pos


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


def getPlaneTilt(plane, doPrint=False):
    # PFS detector size in px
    x_size = 4096
    y_size = 4176
    pixel_size = 15. #microns
    
    tip = plane[0]/pixel_size
    tilt = plane[1]/pixel_size
    
    tip_mic = tip*x_size*pixel_size
    tilt_mic = tilt*y_size*pixel_size
    if doPrint:
        print(f"Tip {tip:.3e} rad => {tip_mic:.1f} microns")
        print(f"Tilt {tilt:.3e} rad => {tilt_mic:.1f} microns")
    
    return tip, tilt, tip_mic, tilt_mic 

def getPlaneFocus(x,y, plane):
    a = plane[0]
    b = plane[1]
    foc0 = plane[2]
    
    return foc0 - (a*x + b*y) 

def getPlaneDeFocus(plane, doPrint=False):
    # PFS detector size in px
    x_size = 4096
    y_size = 4176    
    defoc = getPlaneFocus((x_size/2),(y_size/2), plane)
    if doPrint:
        print(f"defoc {defoc:.1f} microms")
    return defoc

def getPlaneInfo(plane, doPrint=False):
    return getPlaneTilt(plane, doPrint=doPrint), getPlaneDeFocus(plane, doPrint=doPrint)




#####
#####

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

