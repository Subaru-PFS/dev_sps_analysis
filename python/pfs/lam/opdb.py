from opdb import utils, opdb

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
    sql_all = f"select sps_sequence.visit_set_id, sequence_type, name, comments,           sps_exposure.pfs_visit_id, cmd_str, exp_type, sps_module_id, arm, notes, data_flag,time_exp_start,status from sps_exposure \
    inner join visit_set on sps_exposure.pfs_visit_id=visit_set.pfs_visit_id \
    inner join sps_sequence on visit_set.visit_set_id=sps_sequence.visit_set_id \
    inner join sps_visit on sps_exposure.pfs_visit_id=sps_visit.pfs_visit_id \
    inner join sps_camera on sps_exposure.sps_camera_id = sps_camera.sps_camera_id \
    left outer join sps_annotation on sps_exposure.pfs_visit_id=sps_annotation.pfs_visit_id \
    WHERE sps_sequence.visit_set_id = {expID}"
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
