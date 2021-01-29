import sqlite3
import os
from skimage.feature import register_translation
import pandas as pd
from astropy.io import fits
from datetime import datetime as dt
import numpy as np
from scipy.optimize import curve_fit


# Exposure SAC logbook class

class Logbook:
    engine = '///data/ait/experimentLog-sac.db'
    
    @staticmethod
    def lastExperimentId():

        conn = sqlite3.connect(Logbook.engine)
        c = conn.cursor()
        c.execute("""SELECT MAX(experimentId) FROM Experiment""")
        (experimentId,) = c.fetchone()
        return experimentId

    @staticmethod
    def visitRange(experimentId=False):
        experimentId = Logbook.lastExperimentId() if not experimentId else experimentId
        
        conn = sqlite3.connect(Logbook.engine)
        c = conn.cursor()

        c.execute('''select visitStart,visitEnd from Experiment where ExperimentId=%i''' % experimentId)
        (visitStart, visitEnd) = c.fetchone()

        return visitStart, visitEnd
    
    @staticmethod
    def cmdStr(experimentId=False):
        experimentId = Logbook.lastExperimentId() if not experimentId else experimentId
        
        conn = sqlite3.connect(Logbook.engine)
        c = conn.cursor()

        c.execute('''select cmdStr from Experiment where ExperimentId=%i''' % experimentId)
        (cmdStr,) = c.fetchone()

        return cmdStr
    
    @staticmethod
    def getParameter(experimentId, param, doRaise=True):
        cmdStr = Logbook.cmdStr(experimentId=experimentId)
        cmdStr = cmdStr.replace(', ', ',').replace(' ,', ',')
        res = cmdStr.split('%s='%param)
        if len(res)==1:
            if doRaise:
                raise ValueError('parameter %s in not in the command'%param)
            else:
                return ''
            
        return res[1].split(" ")[0]
    
    @staticmethod
    def getDuplicate(experimentId):
        try:
            duplicate = int(Logbook.getParameter(experimentId=experimentId, param='duplicate'))
        except ValueError:
            duplicate = 1
            
        return duplicate


def constructFilelist(visitStart, visitEnd):
    visits = {}
    directory = "/data/ait/sac"

    datefolders = os.listdir(directory)
    datefolders.remove('nextSeqno')

    for datefolder in datefolders:
        files = os.listdir(os.path.join(directory, datefolder))
        for file in files:
            try: 
                visits[int(file[3:9])] = os.path.join(directory, datefolder, file)

            except:
                pass
    
    filelist = []
    for visit in range(visitStart, visitEnd+1):
        if visit not in visits.keys():
            print ('visit %i is missing !'%visit)
        else:
            filelist.append(visits[visit])
    
    filelist.sort()
    
    return filelist

def stackedImage(filelist, ind, duplicate, bck=False, doMeanBck=False):
    sublist = filelist[ind*duplicate:ind*duplicate+duplicate]
    first = sublist[0]
    img = fits.open(first)
    hdr =  img[0].header
    data = img[0].data
    res = np.zeros(data.shape, dtype=np.float32)
    
    for filepath in sublist:
        img = fits.open(filepath)
        res += np.copy(img[0].data)
        
    res = res/duplicate
    
    if bck:
        back_im = fits.open(bck)
        back_data = back_im[0].data
        return hdr, res - back_data
        
    if doMeanBck:
        res -= np.median(res)

    return hdr, res

def getDateObs(experimentId):
    visitStart, visitEnd = Logbook.visitRange(experimentId=experimentId)
    start = constructFilelist(visitStart=visitStart, visitEnd=visitEnd)[0]
    img = fits.open(start)
    hdr =  img[0].header
    return hdr['DATE-OBS']

    
