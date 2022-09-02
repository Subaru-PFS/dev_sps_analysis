from opdb import utils, opdb
import pandas as pd
import re

# def getVisitRange(visit_set_id):
#    sql_all = f"select sps_sequence.visit_set_id, sps_exposure.pfs_visit_id from sps_exposure #\
#    inner join visit_set on sps_exposure.pfs_visit_id=visit_set.pfs_visit_id \
#    inner join sps_sequence on visit_set.visit_set_id=sps_sequence.visit_set_id \
#    inner join sps_visit on sps_exposure.pfs_visit_id=sps_visit.pfs_visit_id \
#    WHERE sps_sequence.visit_set_id = {visit_set_id}"
#    
#    df = utils.fetch_query(opdb.OpDB.url, sql_all)
#    visit_min = df.pfs_visit_id.min()
#    visit_max = df.pfs_visit_id.max()
#    return visit_min, visit_max

def fetch_sps_sequence(visit_set_id):
    sql_all = f"select iic_sequence.visit_set_id, sequence_type, name, comments,sps_exposure.pfs_visit_id, cmd_str, exp_type, sps_module_id, arm, notes, data_flag,time_exp_start,status from sps_exposure \
    inner join visit_set on sps_exposure.pfs_visit_id=visit_set.pfs_visit_id \
    inner join iic_sequence on visit_set.visit_set_id=iic_sequence.visit_set_id \
    inner join sps_visit on sps_exposure.pfs_visit_id=sps_visit.pfs_visit_id \
    inner join sps_camera on sps_exposure.sps_camera_id = sps_camera.sps_camera_id \
    left outer join sps_annotation on sps_exposure.pfs_visit_id=sps_annotation.pfs_visit_id \
    WHERE iic_sequence.visit_set_id = {visit_set_id}"
    df = utils.fetch_query(opdb.OpDB.url, sql_all)
    return df

def getVisitRange(visit_set_id):
    query = f"select min(pfs_visit_id) as visitStart, max(pfs_visit_id) as visitEnd from visit_set WHERE visit_set_id  = {visit_set_id}"

    df = utils.fetch_query(opdb.OpDB.url, query)
    visit_min = df.visitstart.values[0]
    visit_max = df.visitend.values[0]
    return int(visit_min), int(visit_max)

def get_Visit_Set_Id(visit):
    query = f'select visit_set_id from visit_set WHERE pfs_visit_id = {visit}'
    df = utils.fetch_query(opdb.OpDB.url, query)
    return df.values[0][0]


#####
# Alternative function using the logbook web page to retrieve informations
#####

def get_Visit_Set_Id_fromWeb(visit, url = "https://people.lam.fr/madec.fabrice/pfs/spsLogbook.html"):
    df = pd.read_html(url, header=0, index_col=0, keep_default_na=True, skiprows=0)[0]
    df.drop(index=df.index[0], inplace=True)
    df.index.name = 'experimentId'
    df.reset_index(inplace=True)
    df = df.astype({"visitStart": int, "visitEnd": int,"experimentId": int})

    return df.query(f"visitStart <= {int(visit)} <= visitEnd").experimentId.values[0]


def getVisitRange_fromWeb(expId, url = "https://people.lam.fr/madec.fabrice/pfs/spsLogbook.html"):
    df = pd.read_html(url, header=0, index_col=0, keep_default_na=True, skiprows=0)[0]
    df.drop(index=df.index[0], inplace=True)
    df.index.name = 'experimentId'
    df.reset_index(inplace=True)
    df = df.astype({"visitStart": int, "visitEnd": int,"experimentId": int})
    visitStart, visitEnd = df[df.experimentId == expId][["visitStart","visitEnd"]].values[0]

    return visitStart, visitEnd

def getDuplicate_fromWeb(expId, url = "https://people.lam.fr/madec.fabrice/pfs/spsLogbook.html"):
    df = pd.read_html(url, header=0, index_col=0, keep_default_na=True, skiprows=0)[0]
    df.drop(index=df.index[0], inplace=True)
    df.reset_index(inplace=True)
    df.rename(columns={"index": "experimentId"}, inplace=True)
    df = df.astype({"experimentId": int})
    df = df[df['experimentId'] == expId]
    res = re.findall('duplicate=\d+', df.cmd_str.values[0])

    try :
        res = int(re.findall('\d+', res[0])[0])
    except :
        res = 1

    return res

def getMultiTimeExp_fromWeb(expId, url = "https://people.lam.fr/madec.fabrice/pfs/spsLogbook.html"):
    df = pd.read_html(url, header=0, index_col=0, keep_default_na=True, skiprows=0)[0]
    df.drop(index=df.index[0], inplace=True)
    df.reset_index(inplace=True)
    df.rename(columns={"index": "experimentId"}, inplace=True)
    df = df.astype({"experimentId": int})
    df = df[df['experimentId'] == expId]
    res = re.findall('exptime=([^\s]+)', df.cmd_str.values[0])

    try :
        res = len(re.findall('\d+', res[0]))
    except :
        res = 1

    return res

