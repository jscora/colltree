import numpy as np
from astropy.table import vstack,Table,Column,setdiff,join
from tqdm import tqdm
from tqdm import tqdm_notebook as tqdm
import astropy.units as u
import astropy.constants as c
import astropy.io 

#set constants
au = 1.495978707e8 #in km
munit = 5.0428958e31 #in grams
mearth = 5.972e24*1e3 #in g
msol = 1.989e30*1e3 #in g

def get_vel(base_dir,fdir):
    '''Get a table of giant collisions that combines info from follow.maxcorecollisions and important_velcoities.out
    
    Input:
    base_dir = directory where runs are
    fdir = run directory
    
    Output:
    colls'''

    #read in data
    coll = Table.read(base_dir+fdir+'/follow.maxcorecollisions-nograzefa',format='ascii.no_header',
                      names=('time','a','iac','iap','tmass','CMFt','il','ilp','pmass','CMFp','itype','iLR',
                             'LRMass','CMFLR','iSLR','SLRMass','CMFSLR','inew','ideb','mdeb'))
    collmin = Table.read(base_dir+fdir+'/follow.mincorecollisions-nograzefa',format='ascii.no_header',
                      names=('time','a','iac','iap','tmass','CMFt','il','ilp','pmass','CMFp','itype','iLR',
                             'LRMass','CMFLR','iSLR','SLRMass','CMFSLR','inew','ideb','mdeb'))

    if len(coll) != len(collmin):
        print('Error: collision tables not same length')
        print(len(coll))
        print(len(collmin))
        print(fdir)
        #exit code here

    #combine min and max info into one table
    coll['CMFt_min'] = collmin['CMFt'] 
    coll['CMFp_min'] = collmin['CMFp']
    coll['CMFLR_min'] = collmin['CMFLR']
    coll['CMFSLR_min'] = collmin['CMFSLR']


    #read in velocity data
    vel = Table.read(base_dir+fdir+'/important_velocities.out',format='ascii.no_header',
                 names=('id1','id2','ctype','gamma','b','bcrit','time','vesc','vimp','vescalpha','veros',
                       'vcat','vsup','vhr','reveros','revsup'))

    #match up velocity info with collision info
    vround = np.vectorize(fround)
    vel['time'] = vround(vel['time'])
    tab = join(coll,vel,keys=['time'],join_type='left')

    #read in mass to get giant collisions from
    mtiny = np.genfromtxt(base_dir+fdir+'/continue.in')
    mask8 = tab['pmass'] >= mtiny[2]

    colls = tab[mask8]

    return(colls)


def get_allvel(base_dir,dir_tab,fname,ovw):
    '''Get a table of giant collision info, including impact parameters, for a list of directories
    
    Input:
    base_dir = directory where runs are
    dir_tab = table with runs to use
    fname = name of output file of giant collision info
    ovw = True if you want to overwrite any pre-existing files with that name
    
    Output:
    tot_colls = table with all collision data'''

    #define table
    tot_colls = Table(names=('dir','dloss','ecc','inc','slope','time','a','iac','iap','tmass','CMFt','CMFt_min','il','ilp','pmass','CMFp','CMFp_min','itype','iLR','LRMass','CMFLR','CMFLR_min','iSLR','SLRMass','CMFSLR','CMFSLR_min','inew','ideb','mdeb','id1','id2','ctype','gamma','b','bcrit','vesc','vimp','vescalpha','veros','vcat','vsup','vhr','reveros','revsup'),
                      dtype=('U100','U12','U12','U12','U12','float32','float32','int32','int32','float32','float32','float32','int32','int32','float32','float32','float32','int32','int32','float32','float32','float32','int32','float32','float32','float32',
                             'int32','int32','float32','int32','int32','int32','float32','float32','float32','float32','float32','float32','float32','float32','float32','float32','float32','float32'))

    for i in range (0,len(dir_tab['dirs'])):
        colls = get_vel(base_dir,dir_tab['dirs'][i])

        #add columns with init params
        dirs_c = Column([dir_tab['dirs'][i]]*len(colls))
        dloss = Column([dir_tab['dloss'][i]]*len(colls))
        ecc = Column([dir_tab['ecc'][i]]*len(colls))
        inc = Column([dir_tab['inc'][i]]*len(colls))
        slope =Column([dir_tab['slope'][i]]*len(colls))
        colls.add_column(dirs_c,name='dir',index=0)
        colls.add_column(dloss,name='dloss',index=1)
        colls.add_column(ecc,name='ecc',index=2)
        colls.add_column(inc,name='inc',index=3)
        colls.add_column(slope,name='slope',index=4)

        tot_colls = vstack([tot_colls,colls])

    #write table to a csv file in base_dir
    if ovw:
        tot_colls.write(base_dir+fname,format='ascii.csv',overwrite=True)

    return(tot_colls)


def iterate_clist(pcoll,clist,master_ids,new_ids,iname):
    """Iterate through collisions list and add to collision history if not duplicate
    Input:
    pcoll = collision history table for this planet
    clist = list of collisions to search through
    master_ids = embryos that make up this planet
    new_ids = new planet ids in this list 
    iname = name of ids to add to new_ids

    Output:
    pcoll = updated version
    new_ids = updated version"""

    #add in all collisions that are not already in the collision table
    #add pid column
    pid_col = Column(p*np.ones(len(clist)))
    clist.add_column(pid_col, name='pid',index=0)

    #iterate through collision list
    for j in range(0,len(clist)):
        #search through only collisions for this planet
        maskp = pcoll['pid'] == p
        mask_collpid = []
        for r in pcoll[maskp]:
            mask_row = []
            for keys in pcoll.keys():
                mask_row += [r[keys] == clist[keys][j]]
            mask_collpid += [np.prod(mask_row)]
        duplicate = np.sum(mask_collpid)

        if not duplicate: #if the collision is not already in the history
            maskt = pcoll['time'] == clist[j]['time']
            pcoll.add_row(clist[j]) #add to pcoll

            if clist['itype'][j] == 4 or clist['itype'][j] >= 8:
                #add in new ids, otherwise don't care about other id coll history
                #if the secondary id is new, add to list
                if clist[iname][j] not in master_ids['id']:
                    new_ids.add_row([clist[iname][j],clist['time'][j]])
                #if secondary id has later time, add to list
                else:
                    mask_master = master_ids['id'] == clist[iname][j]
                    if master_ids[mask_master]['time'] < clist['time'][j]:
                        new_ids.add_row([clist[iname][j],clist['time'][j]])
    return(pcoll,new_ids)



def find_prev_coll(pcoll,ids,master_ids,collnew,p):
    """Function that finds previous collisions for a given list of planet ids.
    Inputs:
    pcoll = the table of collisions that are already in the collision history table
    ids = table of ids and their collision times to find previous collisions for
    master_ids = table of ids that already have previous collisions found
    collnew = list of collisions for the whole run
    p = the final planet id

    Output:
    pcoll = updated version of the pcoll table"""


    #table to put in ids that need to find previous collisions for
    ids_to_find = Table(names=('id','time'))
    #table of ids that are new in this iteration, used to make ids_to_find
    new_ids = Table(names=('id','time'))

    #iterate through ids
    for i in range(0,len(ids)):

        mask1 = collnew['iap'] == ids['id'][i] #check if they are targets
        mask2 = collnew['ilp'] == ids['id'][i] #check if they are projectiles
        maskt = collnew['time'] < ids['time'][i] #make sure it happens before the child collision

        #get list of new collisions to consider
        clist1 = collnew[mask1&maskt]
        clist2 = collnew[mask2&maskt]

        if len(clist1) > 0:

            pcoll, new_ids = iterate_clist(pcoll,clist1,master_ids,new_ids,'ilp')

        if len(clist2) > 0:

            pcoll, new_ids = iterate_clist(pcoll,clist2,master_ids,new_ids,'iap')

    #cut this down to unique values and earliest time
    if len(new_ids) > 1:

        ids_group = new_ids.group_by('id')
        #pick time that is largest for each group
        for key, group in zip(ids_group.groups.keys,ids_group.groups):
            new_time = np.max(group['time'])

            #if it's not already there, add to master list
            if key['id'] not in master_ids['id']:
                master_ids.add_row([key['id'],new_time])

                #if it's not already there, need to look for new collisions
                ids_to_find.add_row([key['id'],new_time])
            elif  key['id'] in master_ids['id']:
                mask_master = master_ids['id'] == key['id']
                #if this collisions is at an earlier time than the one you have, update time
                if master_ids[mask_master]['time'] < new_time:
                    master_ids['time'][np.where(master_ids['id'] == key['id'])] = new_time

                    #add to collisions to look for later collisions
                    ids_to_find.add_row([key['id'],new_time])



    #if there's only one new id, just check it against master list
    elif len(new_ids) == 1:
        if new_ids['id'][0] not in master_ids['id']:
            master_ids.add_row([new_ids['id'][0],new_ids['time'][0]])

            #need to find new collisions for this id
            ids_to_find.add_row([new_ids['id'][0],new_ids['time'][0]])
        elif  new_ids['id'][0] in master_ids['id']:
                mask_master = master_ids['id'] == new_ids['id'][0]
                if master_ids[mask_master]['time'] < new_ids['time'][0]:
                    master_ids['time'][np.where(master_ids['id'] == new_ids['id'][0])] = new_ids['time'][0]

                    #add to collisions to look for later collisions
                    ids_to_find.add_row([new_ids['id'][0],new_ids['time'][0]])

    #find new collisions if there are any new ids
    if len(ids_to_find) == 0:
        return(pcoll,master_ids)
    elif len(ids_to_find) > 0:
        return(find_prev_coll(pcoll,ids_to_find,master_ids,collnew,p))


def get_small_coll(base_dir,dirn,scoll,emb_list,mtiny,p):
    """Get the number of small collisions that occur over a planet's formation

    Input:
    base_dir = directory that folders are in
    dirn = folder name for this run
    scoll = table of the small collisions that happen to the planet
    emb_list = a table of the embryo ids that make up a planet's history, and time they entered that history
    mtiny = minimum embryo mass
    p = planet id

    Output:
    scoll = table of collision data """
    
    #read in collision info
    colls = Table.read(base_dir+dirn+'/follow.maxcorecollisions-nograzefa',format='ascii.no_header',
                      names=('time','a','iac','iap','tmass','CMFt','il','ilp','pmass','CMFp','itype','iLR',
                             'LRMass','CMFLR','iSLR','SLRMass','CMFSLR','inew','ideb','mdeb'))
    colls_min = Table.read(base_dir+dirn+'/follow.mincorecollisions-nograzefa',format='ascii.no_header',
                      names=('time','a','iac','iap','tmass','CMFt','il','ilp','pmass','CMFp','itype','iLR',
                             'LRMass','CMFLR','iSLR','SLRMass','CMFSLR','inew','ideb','mdeb'))
    #add min CMF to main coll table
    colls.add_column(colls_min['CMFt'],name='CMFt_min',index=6)
    colls.add_column(colls_min['CMFp'],name='CMFp_min',index=11)
    colls.add_column(colls_min['CMFLR'],name='CMFLR_min',index=16)

    #get table with only small collisions
    maskm = colls['pmass'] < mtiny
    colls_small = colls[maskm]
    
    #get debris origin here
    debris = Table.read(base_dir+dirn+'/debris.origin',format='ascii.no_header',
                        names=('time','d_iinit','origin'))
    
    
    for m in emb_list['id']:
        #only need ones where target is the same, we're looking at debris projectiles
        tmask = colls_small['iap'] == m
        slist = colls_small[tmask]
        
        if len(slist)>0:
        
            origin_time = [] ; origin_id = []
            #add debris origin and time
            for i in range(0,len(slist)):
                maskdi = debris['d_iinit'] == slist['ilp'][i]
                
                if len(debris[maskdi]) <= 0:
                    origin_time.append(0)
                    origin_id.append(0)
                else:
                    origin_time.append(debris[maskdi]['time'][0])
                    origin_id.append(debris[maskdi]['origin'][0])
        

            #add debris origin info
            slist.add_column(origin_time,name='origin_time')
            slist.add_column(origin_id,name='origin_id')
        
        
            #add pid
            slist.add_column([p]*len(slist),name='pid',index=0)
            slist.add_column(['other']*len(slist),name='origin_type')
            scoll = vstack([scoll,slist]) 
        
    #origin type
    #sort into three groups: re-accretion, accretion from other bodies, 
    #or specifically accretion from planets destroyed by supercatastrophic disruption
    mask = np.isin(scoll['origin_id'],emb_list['id'])
    scoll['origin_type'][mask] = 'reacc'
    
    mask1 = colls['itype'] == 1
    ids1 = np.concatenate((np.array(colls[mask1]['iap']),np.array(colls[mask1]['ilp'])),axis=0)
    maskin1 = np.isin(scoll['origin_id'],ids1)
    scoll['origin_type'][maskin1] = 'supercat'
        
    
    return(scoll)

def calc_coll_all(base_dir,fdir,coll,cparam,minemb):
    """Calculating collisional history for all final planets. Calls find_prev_coll
    Inputs:
    base_dir = directory where runs are
    fdir = directory of simulation
    cparam = collision size to do collhist for. Options are: 'giant','small', or 'all'.
    coll = table of all giant collisions 
    minemb = minimum embryo mass
    
    Outputs:
    pcoll = table of all giant collision histories
    scoll = table of all small collision histories"""
    
    print(fdir)

    dtsyn = Table(names=('dt','tloss1','tloss2'))
    DT = []
    pcoll = Table(names=('pid','dir','dloss','ecc','inc','slope','time','a','iac','iap','tmass',
                         'CMFt','CMFt_min','il','ilp','pmass','CMFp','CMFp_min','itype','iLR','LRMass',
                         'CMFLR','CMFLR_min','iSLR','SLRMass','CMFSLR','CMFSLR_min','inew','ideb','mdeb',
                         'id1','id2','ctype','gamma','b','bcrit','vesc','vimp','vescalpha','veros','vcat',
                         'vsup','vhr','reveros','revsup'),
                 dtype=('int64','U100','U12','U12','U12','U12','float64','float64','int64','int64',
                        'float64','float64','float64','int64','int64','float64','float64','float64','int64',
                        'int64','float64','float64','float64','int64','float64','float64','float64','int64',
                        'int64','float64','int64','int64','int64','float64','float64','float64','float64',
                        'float64','float64','float64','float64','float64','float64','float64','float64'))

    scoll = Table(names=('pid','time','a','iac','iap','tmass',
                         'CMFt','CMFt_min','il','ilp','pmass','CMFp','CMFp_min','itype','iLR','LRMass',
                         'CMFLR','CMFLR_min','iSLR','SLRMass','CMFSLR','inew','ideb','mdeb',
                         'origin_time','origin_id','origin_type'),
                  dtype=('int64','float64','float64','int64','int64','float64','float64',
                         'float64','int64','int64','float64','float64','float64','int64',
                         'int64','float64','float64','float64','int64','float64',
                         'float64','int64','int64','float64','float64','int64','U24'))

    #read in comp file to get final planets
    comp = Table.read(base_dir+fdir+'/pl.maxcorecompositions-nograzefa',format='ascii.no_header',
                  names=('time','iinit','a','e','inc','mass','inew','mcore','mmant','mtot'))

    print(len(coll))


    mtiny = np.genfromtxt(base_dir+fdir+'/continue.in')

    #make list of final planets
    plnts = comp.group_by('time')
    l = len(plnts.groups)

    #planets here are all embryos above some min embryos mass
    allplanets = plnts.groups[l-1]
    if minemb == 'embryo':
        maskmp = allplanets['mass'] > mtiny[2]*munit/mearth
    else:
        maskmp = allplanets['mass'] >= minemb
    planets = allplanets[maskmp]


    for p in planets['iinit']:

        #make table for the ids you're searching for
        ids = Table(names=('id','time'))
        ids.add_row([p,1e8])
        #make table for the master list of ids and add planet id
        master_ids = Table(names=('id','time'))
        master_ids.add_row([p,1e8])

        lp = len(pcoll)
        lm = len(master_ids)
        ls = len(scoll)

        if cparam == 'all' or cparam == 'giant':
            #find all collisions
            maskm = coll['tmass'] >= mtiny[2]
            maskm2 = coll['pmass'] >= mtiny[2]
            collnew = coll[maskm&maskm2]
            pcoll,master_ids = find_prev_coll(pcoll,ids,master_ids,collnew,p)

        if cparam == 'all' or cparam == 'small':
            #here is where we figure out how many small collisions have happened to a body
            scoll = get_small_coll(base_dir,fdir,scoll,master_ids,mtiny[2],p)

        if cparam == 'all' or cparam == 'giant':
            if len(pcoll) == lp:
                if len(master_ids) == lm:
                    #add a line for this id that has no collisions
                    dloss = collnew['dloss'][0]
                    ecc = collnew['ecc'][0]
                    inc = collnew['inc'][0]
                    slope = collnew['slope'][0]
                    time = planets['time'][0]
                    ip = np.where(planets['iinit'] == p)
                    a = planets['a'][ip]
                    mass = planets['mtot'][ip]
                    pcoll.add_row([p,fdir,dloss,ecc,inc,slope,time,a,0,p,mass,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
                elif len(master_ids) != lm:
                    print('Error: pcoll same length but not master_ids')
                    print(lm,len(master_ids))
                    print(lp,len(pcoll))
                    print(master_ids[lm-1])
                    print(pcoll[lp-1])

        if cparam == 'all' or cparam == 'small':
            if len(scoll) == ls:
                #add a line for this id that has no collisions
                ip = np.where(planets['iinit'] == p)
                scoll.add_row([p,planets['time'][0],planets['a'][ip],0,p,planets['mtot'][ip],
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,p,'none'])


    return(pcoll,scoll)

def get_collhist(base_dir,dir_tab,vtab,cparam,minemb,fname,fname_s,ovw):
    """Calls calc_coll_all for dirs. Creates and writes a table for small collisions and giant collisions.
        
        Input:
        base_dir = directory where data is
        dir_tab = table of directory names and parameters that you want your table generated from
        vtab = table of all giant collisions
        cparam = collision size to do collhist for. Options are: 'giant','small', or 'all'
        minemb = minimum embryo mass. Options are 'embryo' or the minimum embryo mass you want.
        fname = name for giant collision history table
        fname_s = name for small collision history table
        ovw = if True, overwrite any existing files of that name
        
        Output:
        collhist = giant collision history table
        collhist_small = small collision history table"""
    
    collhist = Table(names=('pid','dir','dloss','ecc','inc','slope','time','a','iac','iap','tmass','CMFt','CMFt_min','il','ilp','pmass','CMFp','CMFp_min',
                        'itype','iLR','LRMass','CMFLR','CMFLR_min','iSLR','SLRMass','CMFSLR','CMFSLR_min','inew','ideb','mdeb',
                            'id1','id2','ctype','gamma','b','bcrit','vesc','vimp','vescalpha','veros','vcat','vsup','vhr','reveros','revsup'),
                    dtype=('int64','U100','U12','U12','U12','U12','float64','float64','int64','int64','float64','float64','float64',
                           'int64','int64','float64','float64','float64','int64','int64','float64','float64','float64','int64','float64','float64','float64','int64',
                           'int64','float64','int64','int64','int64','float64','float64','float64','float64','float64','float64','float64','float64',
                           'float64','float64','float64','float64'))
    
    collhist_small = Table(names=('pid','dir','dloss','ecc','inc','slope','time','a','iac','iap','tmass','CMFt','CMFt_min',
                                  'il','ilp','pmass','CMFp','CMFp_min','itype','iLR','LRMass','CMFLR','CMFLR_min',
                                  'iSLR','SLRMass','CMFSLR','inew','ideb','mdeb','origin_time','origin_id'),
                    dtype=('int64','U100','U12','U12','U12','U12','float64','float64','int64','int64','float64','float64',
                           'float64','int64','int64','float64','float64','float64','int64','int64','float64','float64',
                           'float64','int64','float64','float64','int64','int64','float64','float64','int64'))
    
    
    for i in range (0,len(dir_tab)):
        #iterate through directories in the table
        dmask = vtab['dir'] == dir_tab['dirs'][i] #get velocity table for that run
        pcoll, scoll = calc_coll_all(base_dir,dir_tab['dirs'][i],vtab[dmask],cparam,minemb)
        collhist = vstack([collhist,pcoll])
        
        #add columns with init params from dir_tab
        dirs_c = Column([dir_tab['dirs'][i]]*len(scoll))
        dloss = Column([dir_tab['dloss'][i]]*len(scoll))
        ecc = Column([dir_tab['ecc'][i]]*len(scoll))
        inc = Column([dir_tab['inc'][i]]*len(scoll))
        slope =Column([dir_tab['slope'][i]]*len(scoll))
        scoll.add_column(dirs_c,name='dir',index=1)
        scoll.add_column(dloss,name='dloss',index=2)
        scoll.add_column(ecc,name='ecc',index=3)
        scoll.add_column(inc,name='inc',index=4)
        scoll.add_column(slope,name='slope',index=5)
        collhist_small = vstack([collhist_small,scoll])
    
    if cparam == 'giant' or cparam == 'all':
        collhist.write(base_dir+fname,format='ascii.csv',overwrite=ovw)
    if cparam == 'small' or cparam == 'all':
        collhist_small.write(base_dir+fname_s,format='ascii.csv',overwrite=ovw)
    
    return(collhist,collhist_small)
