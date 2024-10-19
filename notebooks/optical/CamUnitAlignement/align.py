import math
import matplotlib.pyplot as plt
import numpy as np



RfrontBell=372.5 #rayon en mm
VPHG_Set = 0.867 # arsec / µm


# Define center wavelenght to the corresponding arm
switchWave = {
    "r":[795.0362,826.6794],
    "b":[577.121, 579.2276],
    "n":[1087.788, 1119.017]
    }

switchWaveRange = {
    "r":[630,970],
    "b":[380, 650],
    "n":[940, 1260]
    }


def getRollCu(df, extremFibers=(2,650), doPrint=False, arm=None):
    # This shift is calculated using position of a center wavelength for the 2 extreme fibers (2 and 650)
    # TBO: for visible camera 
    # sign is given by doing differnce between 2 and 650
    # for nir it is between 650 and 2 
    i0 = 0
    i1 = 1
    if arm == "n":
        i0 = 1
        i1 = 0
    cu_roll_shift = df[df.fiber == extremFibers[i0]].y.values - df[df.fiber == extremFibers[i1]].y.values
    distance =  df[df.fiber == extremFibers[0]].x.values - df[df.fiber == extremFibers[1]].x.values

    # Calculation of shift
    roll_angle = math.atan(np.mean(cu_roll_shift/distance))
    cu_roll_shift = cu_roll_shift.mean()

    roll_shim = RfrontBell * math.tan(roll_angle)
    annot = "roll_shift %.1f px\n"%cu_roll_shift
    annot += "roll_angle %.1f arsec \n"%(math.degrees(roll_angle)*3600)
    annot += "roll_shim %.2f mm\n"%roll_shim
    
    if doPrint:
        print(distance)
        print(roll_angle)
        print(annot)
        print()
        
    if abs(math.degrees(roll_angle)*3600) < 50:
        if doPrint:
            print(f"{abs(math.degrees(roll_angle)*3600)} < 50 arsec : alignment is ok")
    return annot, roll_shim, roll_angle

def plot_CU_roll(df, dataId, midWaves, fname=None, doSavePlot=False, site='LAM', doPrint=False, px2fiber=None, fiber2px=None):
    
    annot, roll_shim, roll_angle = getRollCu(df, doPrint=doPrint, arm=dataId['arm'])
    
    fig, ax = plt.subplots(figsize=(12,8))
    plt.grid()
    plt.title(f"{site} - SM{dataId['spectrograph']} - {dataId['arm'].upper()}{dataId['spectrograph']} - CU Roll adjustment \n visitId = {dataId['visit']}")

    lns1 = ax.plot("x", "y", "*-", label=f"{midWaves[0]}nm", data=df[df.wavelength == midWaves[0]])
    secax = ax.secondary_xaxis(1, functions=(px2fiber, fiber2px))
    secax.set_xlabel('Fiber')
    lns = lns1
    if len(midWaves) == 2 :
        axy = ax.twinx()
        lns2 = axy.plot("x", "y", "*-", label=f"{midWaves[1]}nm", data=df[df.wavelength == midWaves[1]], color="r")

        # added these three lines
        lns = lns1+lns2
    
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)

    #ax = df.groupby("wavelength")[["x","y"]].plot(kind="scatter")
    #ax.plot("x", "y", data=df.groupby("wavelength")) #  "*-", label=f"{midWave}nm")
    color=None
    if abs(math.degrees(roll_angle)*3600) <= 50 :
        color = "green"
    ax.annotate(annot,  xy=(.45,.45),xycoords="figure fraction", color=color)
    if doSavePlot:
    #fig.patch.set_alpha(0.5)
        plt.savefig(fname)
#        plt.savefig(fname, transparent=True)


def plot_vphg_roll(df, dataId, fibers, fname=None, doSavePlot=False, site='LAM'):
    
    fig, ax = plt.subplots(figsize=(12,8))
    plt.grid()

    plt.xlabel("X pixel")
    plt.ylabel("Y pixel")
    plt.title(f"{site} - SM{dataId['spectrograph']} - {dataId['arm'].upper()}{dataId['spectrograph']} - VPHG Roll adjustment \n visitId = {dataId['visit']}")


    x0 = df.x.values[0] 

    lns1 = ax.plot( df.x, df.y, 'o', label=f"{*fibers,}", color="black")

    # fit data 
    popt = np.polyfit( df.x, df.y, 1)
    print(popt)
    t = np.arange(df.x.min(), df.x.max()+1)
    lns2 = ax.plot(t,np.polyval(popt, t), 'b--', label="linear fit")

    annot = f"Shift { df.x.max() - df.x.min():.1f} px over {df.y.max() - df.y.min():.1f} px \n" 

    Roll_VPHG_shim_fit = np.rad2deg(1/popt[0])*3600*VPHG_Set
    annot += f"VPHG Roll Shim: {Roll_VPHG_shim_fit:.0f} um \n"

    print(np.rad2deg(popt[0])*3600*VPHG_Set)

    print(annot)

    plt.annotate(annot, xy=(.58,.2),xycoords="figure fraction")

    lns = lns1+lns2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)


    if doSavePlot:
    #fig.patch.set_alpha(0.5)
        plt.savefig(fname)
#        plt.savefig(fname, transparent=True)
    return Roll_VPHG_shim_fit



def plotDetAlign(detMap, specId=None, waves=None):
    
    dMapInfo = detMap.metadata.toDict()
    detMapvisitId = dMapInfo['CALIB_INPUT_0']
    cam = f"{dMapInfo['ARM']}{specId}"
    
    pix_size=15
    waveLow, waveMid, waveHigh = waves
    
    ft = 2 + (specId-1)*651
    xMed_fib2 = int(detMap.findPoint(ft,waveMid).getX()) #ftFunction.yCenter + ftFunction.yLow

    ft = 650 + (specId-1)*651
    xMed_fib650 = int(detMap.findPoint(ft,waveMid).getX()) #ftFunction.yCenter + ftFunction.yLow

    naxis1, naxis2 = detMap.bbox.getDimensions()
    print(naxis1 - xMed_fib2, xMed_fib650)

    # decenter is given by comparing the margin there between the side of the detector and spectrum of the 2 extermes sciences fibers 2 and 650
    # xdim so naxis1 - x of fiber #2 for center wave 
    # 

    spatial_decenter = (naxis1 - xMed_fib2) - xMed_fib650
    print(f"spatial decenter: {spatial_decenter}px")

    ft = 650 + (specId-1)*651
    yHigh_fib339 = int(detMap.findPoint(ft,waveLow).getY())
    yLow_fib339 = int(detMap.findPoint(ft,waveHigh).getY())

    spectral_decenter_ft650 = (naxis2 - yHigh_fib339) - yLow_fib339
    print(f"spectral decenter: {spectral_decenter_ft650}px")

    ft = 2 + (specId-1)*651
    yHigh_fib339 = int(detMap.findPoint(ft,waveLow).getY())
    yLow_fib339 = int(detMap.findPoint(ft,waveHigh).getY())

    spectral_decenter_ft2 = (naxis2 - yHigh_fib339) - yLow_fib339
    print(f"spectral decenter: {spectral_decenter_ft2}px")

    
    
    
    
    
    
    
    
    
    fig, ax = plt.subplots(figsize=(12,8))
    ft = 2 + (specId-1)*651
    xCenters = detMap.getXCenter(ft)
    yLow = int(detMap.findPoint(ft,waveLow).getY()) #ftFunction.yCenter + ftFunction.yLow
    yHigh = int(detMap.findPoint(ft,waveHigh).getY())
    print(yLow, yHigh)
    offset = 0
    print(detMap.findWavelength(ft,yLow+offset ))
    ax.plot(xCenters[yLow+offset:yHigh], np.arange(yLow, yHigh), color="red")
    ax.annotate(f"fiber {ft}",  xy=(.85,.5),xycoords="figure fraction", color="red")
    ax2 =ax.twinx()
    xCenters = detMap.getXCenter(ft)
    w = detMap.getWavelength(ft)
    ax2.plot(xCenters[yLow+offset:yHigh], w[yLow+offset:yHigh], color="red")


    ft = 650 + (specId-1)*651
    xCenters = detMap.getXCenter(ft)
    yLow = int(detMap.findPoint(ft,waveLow).getY()) #ftFunction.yCenter + ftFunction.yLow
    yHigh = int(detMap.findPoint(ft,waveHigh).getY())
    print(yLow, yHigh)
    offset = 0
    print(detMap.findWavelength(ft,yLow+offset ))
    ax.plot(xCenters[yLow+offset:yHigh], np.arange(yLow, yHigh), color="blue")
    ax.annotate(f"fiber {ft}",  xy=(.1,.5),xycoords="figure fraction", color="blue")
    xCenters = detMap.getXCenter(ft)
    w = detMap.getWavelength(ft)
    ax2.plot(xCenters[yLow+offset:yHigh], w[yLow+offset:yHigh], color="blue")

    edgeX = []
    edgeY = []

    for f in detMap.getFiberId():
        edgeX.append(detMap.findPoint(f,waveLow).getX())
        edgeY.append(detMap.findPoint(f,waveLow).getY()) 
    ax.plot(edgeX, edgeY, color="green")
    ax.annotate(f"{waveLow}nm",  xy=(.5,.1),xycoords="figure fraction", color="green")

    #ax3 = ax.twiny()
    #ax3.plot(detMap.getFiberId(), edgeY)
    edgeX = []
    edgeY = []

    for f in detMap.getFiberId():
        edgeX.append(detMap.findPoint(f,waveHigh).getX())
        edgeY.append(detMap.findPoint(f,waveHigh).getY()) 
    ax.plot(edgeX, edgeY, color="purple")
    ax.annotate(f"{waveHigh}nm",  xy=(.5,.9),xycoords="figure fraction", color="purple")

    #ax3.plot(detMap.getFiberId(), edgeY)

    ax.set_xlim(1,4095 )
    ax.set_ylim(1,4176 )
    #ax3.set_xlim(detMap.getFiberId()[-1], detMap.getFiberId()[0])
    ax.set_xlabel("pixel x")
    ax.set_ylabel("pixel y")
    ax2.set_ylabel("Wavelength (nm)")
    plt.annotate(f"Spatial decenter at mid wavelength {waveMid:.1f}nm = {spatial_decenter}px ({spatial_decenter*pix_size:.1f}µm)", (0.3, 0.5), xycoords="figure fraction")
    plt.annotate(f"Spectral decenter at \n(fiber 2, {waveLow:.1f}nm) = {spectral_decenter_ft2}px ({spectral_decenter_ft2*pix_size:.1f}µm)", (0.65, 0.2), xycoords="figure fraction")
    plt.annotate(f"Spectral decenter at \n(fiber 650, {waveLow:.1f}nm) = {spectral_decenter_ft650}px ({spectral_decenter_ft650*pix_size:.1f}µm)", (0.15, 0.2), xycoords="figure fraction")

    plt.title(f"{cam.upper()} centering - detectorMap {detMapvisitId}")
    #plt.savefig(imgPath+f"{cam.upper()}centering-detectorMap{detMapvisitId}"+".png", bbox_inches = "tight")
    
    

def getXY(fiber, wave, detMap, specId=None):
    '''
    fiber in DCB unit so from 1 to 651
    '''
    ft = fiber + (specId-1)*651
    x = detMap.findPoint(ft,wave).getX()
    y = detMap.findPoint(ft,wave).getY()
    return x, y

def getShifts(detMap, specId=None,fibers=[2,339,650], simMap=None, waves=None):
    '''
    get X and Y shift in px 
    between simulatated detectorMap (so the goal) and the current detectorMap
    waves: [waveLow, waveMid, waveHigh]
    '''
    
    pos = []
    for fib in fibers:
        for wave in waves:
            xs, ys = getXY(fib,wave, simMap, specId=specId )
            x, y = getXY(fib,wave, detMap, specId=specId )
            pos.append({"fiber": fib, "wave": wave, "x":x, "y":y, "dx":xs-x, "dy":ys-y})
    
    return pos

def plotDetMapContour(detMap, specId=None, ax=None, sim=False, colors=["red", "blue"], fibers=[1,651], waves=None):
    if sim is False:
        ls = '-'
        alpha = 1
    else:
        ls = '--'
        alpha = 0.3
    
    xxx = float('NaN')
    
    infos = [
        {"fiber": fibers[0], "color": "red", "annot": (.85,.5)},
        {"fiber": fibers[1], "color": "blue", "annot": (.1,.5)}
        ]
    
    fibers = [d['fiber'] for d in infos if 'fiber' in d]
    
    waveLow, waveMid, waveHigh = waves
    for fib in fibers:
        inf = [item for item in infos if item["fiber"] == fib][0]
        color = inf["color"]
        annot = inf["annot"]
        ft = fib + (specId-1)*651
        print(f"fiber: {ft}")
        xCenters = detMap.getXCenter(ft)
        yLow = detMap.findPoint(ft,waveLow).getY()
        yHigh = detMap.findPoint(ft,waveHigh).getY()      
        yLow = int(yLow) if yLow != xxx else yLow
        yHigh = int(yHigh) if yHigh != xxx else yHigh
     
        ax.plot(xCenters[yLow:yHigh], np.arange(yLow, yHigh), color=color, linestyle=ls, alpha=alpha)
            
        if sim is False:
            ax.annotate(f"fiber {ft}",  xy=annot,xycoords="figure fraction", color=color)
            
    edgeX = []
    edgeY = []
    for f in detMap.getFiberId():
        edgeX.append(detMap.findPoint(f,waveLow).getX())
        edgeY.append(detMap.findPoint(f,waveLow).getY()) 
    ax.plot(edgeX, edgeY, color="green", linestyle=ls, alpha=alpha)
    if sim is False:
        ax.annotate(f"{waveLow}nm",  xy=(.5,.15),xycoords="figure fraction", color="green")


    edgeX = []
    edgeY = []
    for f in detMap.getFiberId():
        edgeX.append(detMap.findPoint(f,waveHigh).getX())
        edgeY.append(detMap.findPoint(f,waveHigh).getY()) 
    ax.plot(edgeX, edgeY, color="purple", linestyle=ls, alpha=alpha)
    if sim is False:
        ax.annotate(f"{waveHigh}nm",  xy=(.5,.8),xycoords="figure fraction", color="purple")

    dim = detMap.bbox.getDimensions()
    ax.set_xlim(0,dim.getX() )
    ax.set_ylim(0,dim.getY() )
    ax.set_xlabel("pixel x")
    out = ax.set_ylabel("pixel y")
    return ax
    

def plotDetAlignAxe(detMap, specId=None, ax=None, simMap=None,fibers=[1,651],  waves=None, fname=None, doSavePlot=False, site='LAM'):
    pix_size = 15

    dMapInfo = detMap.metadata.toDict()
    detMapvisitId = dMapInfo['CALIB_INPUT_0']
    cam = f"{dMapInfo['ARM']}{specId}"
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12,8))
    
    ax = plotDetMapContour(detMap, specId=specId, ax=ax, waves=waves, fibers=fibers)
    if simMap is not None:
        plotDetMapContour(simMap, specId=specId, ax=ax, sim=True, waves=waves, fibers=fibers)

        offsets = getShifts(detMap, specId=specId, simMap=simMap, waves=waves)
        #fibers = list(set([d['fiber'] for d in offsets if 'fiber' in d]))
        #waves = list(set([d['wave'] for d in offsets if 'wave' in d]))
        for offset in offsets:
            fiber = offset["fiber"]
            wave = offset["wave"]
            x = offset["x"]
            y = offset["y"]
            dx = offset["dx"]
            dy = offset["dy"]
            ax.annotate(f"{fiber} {wave:.1f}nm  \n dx={dx:.1f}px \n dy={dy:.1f}px", (int(x), int(y)), xycoords="data")
            #print(f"{fiber} {wave:.1f}nm  dx={dx:.1f}px dy={dy:.1f}px ({int(x)},{int(y)}) xycoords='axes pixels'")

    ax.annotate(dMapInfo["CALIB_ID"], (0.2,0.02), xycoords="figure fraction")
    out = plt.title(f"{cam.upper()} centering - detectorMap {detMapvisitId}", pad=30)
    if doSavePlot:
        plt.savefig(fname, bbox_inches = "tight")
    return out

def getTiltCUShim(shift, doPrint=True):
    '''
    shift in pixel
    '''
    F_CU = 300. #mm
    pix_size = 15e-3
    e = 523.639 #mm entraxe en sperolinders de la CU
    shift_mm = shift * pix_size
    shim = e*shift_mm/F_CU
    if doPrint:
        annot = "CU tilt \n"
        annot += "Shift %.2f px\n"%shift
        annot += "Shift %.2f microns\n"%(shift*15)
        annot += "angle %.3f deg\n"%(math.degrees(math.atan(shift_mm/F_CU)))
        annot += "angle %.3f arsec\n"%(math.degrees(math.atan(shift_mm/F_CU))*3600)
        annot += "shim %.3f mm"%(e*shift_mm/F_CU)
        print(annot)
    return 

def getAbsOffsets(detMap, specId=None,fibers=[2,339,650], waves=None):
    dim = detMap.bbox.getDimensions()
    xl = dim.getX()
    yl = dim.getY()
    infos = [
        {"fiber": fibers[0], "wave":waves[0], "color": "lightblue", "limit":(xl,0)},
        {"fiber": fibers[1], "wave":waves[0], "color": "lightblue", "limit":(xl/2,0)},
        {"fiber": fibers[2], "wave":waves[0], "color": "lightblue", "limit":(0,0)},
        {"fiber": fibers[0], "wave":waves[1], "color": "lightblue", "limit":(xl,yl/2)},
        {"fiber": fibers[2], "wave":waves[1], "color": "lightblue", "limit":(0,yl/2)},
        {"fiber": fibers[0], "wave":waves[2], "color": "lightblue", "limit":(xl,yl)},
        {"fiber": fibers[1], "wave":waves[2], "color": "lightblue", "limit":(xl/2,yl)},
        {"fiber": fibers[2], "wave":waves[2], "color": "lightblue", "limit":(0,yl)},
    ]
        
    pos = []
    for p in infos:
        x, y = getXY(p["fiber"],p["wave"], detMap, specId=specId )
        xs, ys = p["limit"]
        p.update({"x":x, "y":y, "dx":x-xs, "dy":y-ys})
#        pos.append({"fiber": p["fiber"], "wave": p["wave"], "x":x, "y":y, "dx":x-xs, "dy":y-ys})
        pos.append(p)

    
    return pos

def plotReqDetAlignAxe(detMap, specId=None, ax=None, simMap=None,  waves=None, fname=None, doSavePlot=False, site='LAM'):
    pix_size = 15

    dMapInfo = detMap.metadata.toDict()
    detMapvisitId = dMapInfo['CALIB_INPUT_0']
    cam = f"{dMapInfo['ARM']}{specId}"
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12,8))
    
    
    ax = plotDetMapContour(detMap, specId=specId, ax=ax, waves=waves)
    if simMap is not None:
        plotDetMapContour(simMap, specId=specId, ax=ax, sim=True, waves=waves)

        offsets = getShifts(detMap, specId=specId, simMap=simMap, waves=waves)
        #fibers = list(set([d['fiber'] for d in offsets if 'fiber' in d]))
        #waves = list(set([d['wave'] for d in offsets if 'wave' in d]))
        for offset in offsets:
            fiber = offset["fiber"]
            wave = offset["wave"]
            x = offset["x"]
            y = offset["y"]
            dx = offset["dx"]
            dy = offset["dy"]
            #print(f"{fiber} {wave:.1f}nm  dx={dx:.1f}px dy={dy:.1f}px ({int(x)},{int(y)}) xycoords='axes pixels'")
            
    absOffsets = getAbsOffsets(detMap, specId=specId, waves=waves)
    for offset in absOffsets:
        fiber = offset["fiber"]
        wave = offset["wave"]
        x = offset["x"]
        y = offset["y"]
        dx = offset["dx"]
        dy = offset["dy"]
        mx = np.min((x-0, detMap.getBBox().width-x))
        my = np.min((y-0, detMap.getBBox().height-y))
        ax.annotate(f"margin: x {mx:.1f}px, \n  y {my:.1f}px", (int(x), int(y)-5), xycoords="data", size=8)
        ax.annotate(f"{fiber} {wave:.1f}nm  \n dx={dx:.1f}px \n dy={dy:.1f}px", (int(x), int(y)+120), size=9,
                    xycoords="data", annotation_clip=False)
        #arrow = mpatches.FancyArrowPatch(offset["limit"], (x, y),
        #                     mutation_scale=1)
        #ax.add_patch(arrow)
        #ax.annotate(f"{fiber} {wave:.1f}nm  \n dx={dx:.0f}µm \n dy={dy:.0f}µm", (int(x), int(y)), xycoords="data")    

    ax.annotate(dMapInfo["CALIB_ID"], (0.2,0.02), xycoords="figure fraction")
    out = plt.title(f"{cam.upper()} centering - detectorMap {detMapvisitId}", pad=20)
    if doSavePlot:
        plt.savefig(fname, bbox_inches = "tight")
    return out
