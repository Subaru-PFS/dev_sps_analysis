import numpy as np
from sps_engineering_Lib_dataQuery.databasemanager import DatabaseManager
from sps_engineering_Lib_dataQuery.confighandler import loadConf
import pandas as pd
from matplotlib.dates import num2date
colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf']

def getTables(start):
    conf = loadConf(start)
    d = dict([(device.tablename, device) for device in conf])
    print('\n'.join(d.keys()))
    return d


def getConf(start, actors=None, doPrint=True, dbname='archiver'):
    db = DatabaseManager('tron', 5432, '', dbname=dbname)
    db.init()
#    db = DatabaseManager('localhost', 5432, dbname=dbname)
#    db = DatabaseManager('localhost', dbname=dbname )

    conf = dict([(device.tablename, device) for device in loadConf(start)])
    tables = dict([(device.tablename, device) for device in  db.pollDbConf(start)])

    for key in tables.keys():
        tables[key] = conf[key] if key in conf.keys() else tables[key]

    if actors is not None:
        t=[]
        for actor in actors:
            t+=[(key, val) for key,val in tables.items() if actor in key]
        tables = dict(t)

    if doPrint:
        print('\n'.join(['%s : %s' %(key, ','.join(table.labels)) for key, table in sorted(tables.items())]))

    db.close()

    return tables


def extractData(tables, start, end=False, interpolate=True, dbname='archiver'):
    datas = []
    db = DatabaseManager('tron', 5432, '', dbname=dbname)
    db.init()
    d = getConf(start, doPrint=False, dbname=dbname)

    for i, device in enumerate(tables):
        devConf = d[device]
        df = db.dataBetween(devConf.tablename, ','.join(devConf.keys), start=start, end=end)
#        df = df.drop(columns=['id'])
        df.columns = ['id','tai'] + devConf.labels
        datas.append(df)
        print(f'{100*(i+1)/len(tables)}%')

    if interpolate:
        datas = interp(datas)

    return datas


def getColname(col, dff, i=0):
    colname = col if not i else f'{col}{i}'
    if colname in dff.columns:
        return getColname(col, dff, i+1)
    return colname

def interp(dfs, sec=15, tai=None):
    samp_start = np.max([df.tai.values[0] for df in dfs])
    samp_end = np.min([df.tai.values[-1] for df in dfs])
    step = sec * (1/24/3600)
    tai = np.arange(samp_start, samp_end, step) if tai is None else tai
    dff = pd.DataFrame(dict(date=num2date(tai), tai=tai))
    for df in dfs:
        for col in df.columns:
            if col=='tai':
                continue
            colname = getColname(col, dff)
            dff[colname] = np.interp(tai, df.tai, df[col])

    return dff


