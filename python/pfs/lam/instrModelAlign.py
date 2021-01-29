import numpy as np
import math
import pandas as pd

class zemax():
    def __init__(self, cam, filemodel=None, test=False):
        self.cam = cam
        path = "/home/pfs/dev/lib/pfs/input/"

        if filemodel is None:
            if cam == "r1":
#                filemodel = path+"PFS-Design-Red-allfiber.pkl"
                filemodel = path+"PFS-Design-Red-allfiber_Nov2019_R1.pkl"
            elif cam == "b1":
#                filemodel = path+"PFS-DesignCoeff-Blue-allfiber.pkl"
                filemodel = path+"PFS-DesignCoeff-Blue-allfiber_B1.pkl"
            elif cam == "m1":
                filemodel = path+"PFS-Design-Med-allfiber.pkl"
            else:
                print("Error should have model file...")
        if cam == "r1":
            raw = path+"redZemaxModel.csv"
            pixRes = 0.0847 #nm/px
            waveband = [630,970]
        elif cam == "b1":
            raw = path+"blueZemaxModel.csv"
            pixRes = 0.0671 #nm/px
            waveband = [380,650]
        elif cam == "m1":
            raw = path+"medZemaxModel.csv"          
            pixRes = 0.0477 #nm/px
            waveband = [0,0]

        else:
            print("Error should have Med Res file...")
                
        #        pd.DataFrame.__init__(self, data=data, columns=['x', 'y'])
        self.filemodel = filemodel
        self.model = pd.read_pickle(self.filemodel)
        self.rawmodel =  pd.read_csv(raw)
        self.pixRes = pixRes
        self.waveband = waveband

    @property
    def getmodelname(self):
        return(self.filemodel)
    
    def getmodel(self):
        return(self.model)
    
    def getrawmodel(self):
        return(self.rawmodel)   
    
    def getpixres(self):
        return(self.pixRes)   
    
    def wave2pix(self,wave,fiber):
        x = np.polyval(self.model.xs(fiber).cx, wave)
        y = np.polyval(self.model.xs(fiber).cy, wave)
        return x, y
    def pix2wave(self,y,fiber):
        wave = np.polyval(self.model.xs(fiber).p2w_cy, y)
        return wave


def wave2pix(wave,fiber,zemodel):
    x = np.polyval(zemodel.xs(fiber).cx, wave)
    y = np.polyval(zemodel.xs(fiber).cy, wave)
    return x, y

def pix2wave(y,fiber,zemodel):
    wave = np.polyval(zemodel.xs(fiber).p2w_cy, y)
    return wave

def pix2wave2(y,fiber,model_file):
    zemodel = pd.read_pickle(model_file)
    wave = np.polyval(zemodel.xs(fiber).p2w_cy, y)
    return wave


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

# def find_nearest(X, value):
#    return X[np.unravel_index(np.argmin(np.abs(X - value)), X.shape)]


def GetRollCU(df, doPrint=False):
    y1 = df.sort_values(["px"]).py.values[0]
    y2 = df.sort_values(["px"]).py.values[-1]
    x1 = df.sort_values(["px"]).px.values[0]
    x2 = df.sort_values(["px"]).px.values[-1]

    roll_shift = y2 - y1

    dist = x2 - x1
    RfrontBell=372.5 #rayon en mm

    roll_angle = math.atan(roll_shift/dist)
    roll_shim = RfrontBell * math.tan(roll_angle)
    annot = "roll_shift %.2f px\n"%roll_shift
    annot += "roll_angle %.2f arsec \n"%(math.degrees(roll_angle)*3600)
    annot += "roll_shim %.3f mm\n"%roll_shim
    
    if doPrint:
        print(dist)
        print((math.atan(roll_shift/dist)))    
        print(annot)
    
    return roll_shim



def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = np.arange(len(v))
    
    v = np.asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not np.isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN
    
    lookformax = True
    
    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab), np.array(mintab)
