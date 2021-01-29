import sqlite3
import os
from astropy.io import fits
from datetime import datetime as dt
import numpy as np
import pandas as pd
from datetime import timedelta

# from pfs.utils.dummyCableB import *


# Exposure logbook class

class Logbook:
    import os
    try:
        engine=os.environ['ENGINE_PFS_PATH']
    except:
#        engine = '///data/ait/experimentLog.db'
        engine = '///data/drp/fmadec/experimentLog_sm2.db'

    
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
    def date(experimentId=False):
        experimentId = Logbook.lastExperimentId() if not experimentId else experimentId
        
        conn = sqlite3.connect(Logbook.engine)
        c = conn.cursor()

        c.execute('''select startdate from Experiment where ExperimentId=%i''' % experimentId)
        (startdate) = c.fetchone()
        date = startdate[0].split('T')[0]

        return date
    
    @staticmethod
    def cmdStr(experimentId=False):
        experimentId = Logbook.lastExperimentId() if not experimentId else experimentId
        
        conn = sqlite3.connect(Logbook.engine)
        c = conn.cursor()

        c.execute('''select cmdStr from Experiment where ExperimentId=%i''' % experimentId)
        (cmdStr,) = c.fetchone()

        return cmdStr
    
    @staticmethod
    def getParameter(experimentId, param):
        cmdStr = Logbook.cmdStr(experimentId=experimentId)
        res = cmdStr.split('%s='%param)
        if len(res)==1:
            raise ValueError('parameter %s in not in the command'%param)
            
        return res[1].split(" ")[0]
    
    @staticmethod
    def newAnomalies(experimentId, anomalies):
        sqlRequest = 'UPDATE Experiment SET anomalies = "%s" WHERE experimentId=%i' % (anomalies.replace('"', ""),
                                                                                       experimentId)
        Logbook.newRow(sqlRequest=sqlRequest)
        
        
    @staticmethod
    def newRow(sqlRequest):
        conn = sqlite3.connect(Logbook.engine)
        c = conn.cursor()
        try:
            c.execute(sqlRequest)
            conn.commit()

        except sqlite3.IntegrityError:
            pass
        
    @staticmethod
    def getDuplicate(experimentId):
        try:
            duplicate = int(Logbook.getParameter(experimentId=experimentId, param='duplicate'))
        except ValueError:
            duplicate = 1
            
        return duplicate
    
    @staticmethod
    def getCams(visit):
        conn = sqlite3.connect(Logbook.engine)
        c = conn.cursor()
        c.execute('''select spectrograph,arm from Exposure inner join CamExposure on Exposure.exposureId=CamExposure.exposureId where visit=%i'''%visit)

        return ['%s%d'%(arm, specId) for specId,arm in c.fetchall()]
    
    @staticmethod    
    def visitDate(visit):
        conn = sqlite3.connect(Logbook.engine)
        c = conn.cursor()

        c.execute('''select obsdate from Exposure where visit=%i''' % visit)
        (obsdate) = c.fetchone()
        date = obsdate[0].split('T')[0]

        return date

    @staticmethod    
    def visitExperimentId(visit):
        conn = sqlite3.connect(Logbook.engine)
        c = conn.cursor()
        c.execute('''select experimentId from Experiment where %i BETWEEN Experiment.visitStart AND Experiment.visitEnd'''%visit)

        (experimentId,) = c.fetchone()
        return experimentId



def getfilepath(visit, rerun, date, cam, repo, doRaise=False, cluster=False, basePath=None, subaru=False, doPrint=False, verbose=False):
    repo = cam if repo is None else repo
    basePath = "/drp" if basePath is None else basePath
    site = "L" if subaru == False else "S"
    path = "%s/%s/rerun/%s/detrend/calExp/%s" %(basePath,repo, rerun, date)
    if cluster : 
        path = "/net/SRVSTK20C/drp/%s/rerun/%s/detrend/calExp/%s" %(repo, rerun, date)
    if site == "L":
        filepath = '%s/v%s/calExp-%sA%s%s.fits' % (path, str(visit).zfill(7), site, str(visit).zfill(6), cam)
    else :
        filepath = '%s/v%s/calExp-%sA%s%s.fits' % (path, str(visit).zfill(6), site, str(visit).zfill(6), cam)

    
    if doPrint:
        print("filepath = ", filepath)
        print("site = ",site)

    
    if not os.path.isfile(filepath):
        if doRaise:
            raise IOError("%s does not exist" % filepath)
        else:
            if verbose: print("checking the day after")
            datetime_object = dt.strptime(date, '%Y-%m-%d')
            newdate = datetime_object + timedelta(days=1)
            return getfilepath(visit, rerun, newdate.strftime('%Y-%m-%d'), cam, repo, doRaise=True,basePath=basePath, subaru= subaru)
    return filepath, date

def constructFilelist(experimentId, doPrint=False, rerun="pfs", repo=None, doRaise=True, basePath=None, subaru=False):
    filelist = []
    
    visitStart, visitEnd = Logbook.visitRange(experimentId=experimentId)
    date = Logbook.date(experimentId=experimentId)

    for visit in range(visitStart, visitEnd + 1):
        for cam in Logbook.getCams(visit=visit):
            try:
                filepath, date = getfilepath(visit, rerun, date, cam, repo, basePath=basePath, subaru=subaru, doPrint=doPrint)
                filelist.append((filepath, experimentId, cam, date, visit))
            except IOError:
                if doRaise:
                    raise
    if doPrint:
        print("date = ", date)
        print("visitStart = ", visitStart)
        print("visitEnd = ", visitEnd)
#        print("cam = ", cam)
#        print("filepath = ", filepath)
        
    filelist = pd.DataFrame(filelist, columns=['filepath', 'experimentId', 'cam', 'obsdate',"visitId"])
    
    return filelist

def stackedImage(filelist, ind, duplicate, bck=False):
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
        
    return hdr, res

def getFitsKey(filename, key, doRaise=True):
    hdulist = fits.open(filename, "readonly")
    prihdr = hdulist[0].header
    
    try:
        return prihdr[key]
    except KeyError:
        if doRaise:
            raise
        else:
            keyval = np.nan
    
    return keyval

def getFiberIds(fitsfile):
    pfsdsgn = getFitsKey(fitsfile,"W_PFDSGN")
    dcb = DummyCableBDatabase()
    setup = dcb.interpret(pfsdsgn)
    fibers = dcb.getFiberIds(*setup)
    return fibers, setup

def getSourcesUsed(image_file):
    sources = dict(
        neon = "W_AITNEO",
        argon = "W_AITARG",
        hgar = "W_AITHGA",
        krypton = "W_AITKRY",
        xenon = "W_AITXEN",
        deuterium = "W_AITDEU",
        qth = "W_AITQTH"
        )
       
    listsource = []
    for key, fitskey in sources.items():
        try:
            state =  getFitsKey(image_file,fitskey)
            if state in ['off','undef']:
                raise KeyError
        except KeyError:
            state = False
            
        state = True if state=='on' else state
        
        if state==True:
            listsource.append(key)

    return ",".join(listsource)

def getArcLampForNist(lamp, fitsfile=None):
    lamp = lamp if fitsfile is None else getSourcesUsed(fitsfile)
    if lamp == "neon" :
        arclamp = "Ne"
    elif lamp == "argon":
        arclamp = "Ar"
    elif lamp == "hgar":
        arclamp = "Hg | Ar"
    elif lamp == "krypton":
        arclamp = "Kr"
    elif lamp == "xenon":
        arclamp = "Xe"
    
    return arclamp


def getFileInfo(filelist):
    return f"ExpId {filelist.experimentId[0]} -- {(filelist.cam[0].upper())} -- {filelist.obsdate[0]}"


def getvisitfilepath(visit, rerun, date, cam, repo, doRaise=False, cluster=False, doPrint=False):
    repo = cam if repo is None else repo
    date = Logbook.visitDate(visit) if date is None else date
    cam = Logbook.getCams(visit=visit) if cam is None else cam

    path = "/drp/%s/rerun/%s/detrend/calExp/%s" %(repo, rerun, date)
    if cluster : 
        path = "/net/SRVSTK20C/drp/%s/rerun/%s/detrend/calExp/%s" %(repo,rerun, date)
    filepath = '%s/v%s/calExp-LA%s%s.fits' % (path, str(visit).zfill(7), str(visit).zfill(6), cam)
    if doPrint:
        print(filepath)
        print(date)
    if not os.path.isfile(filepath):
        if doRaise:
            raise IOError("%s does not exist" % filepath)
        else:
            datetime_object = dt.strptime(date, '%Y-%m-%d')
            newdate = datetime_object + timedelta(days=1)
            return getvisitfilepath(visit, rerun, newdate.strftime('%Y-%m-%d'), cam, repo, doRaise=True, cluster=cluster)
    
    return filepath, date


def constructVisitFilelist(visit, doPrint=False, rerun="ginga", repo=None, cluster=False, doRaise=True):
    filelist = []
    for cam in Logbook.getCams(visit=visit):
            try:
                filepath, date = getvisitfilepath(visit, rerun, None, cam, repo, cluster=cluster)
                filelist.append((filepath, visit, cam, date))
            except IOError:
                if doRaise:
                    raise
    filelist = pd.DataFrame(filelist, columns=['filepath', 'visitId', 'cam', 'obsdate'])
    
    return filelist

