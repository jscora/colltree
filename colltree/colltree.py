import pandas as pd

def find_children(pid,tloss,tab):
    """Find the children of a planet.
        Input:
        pid = planet id
        tloss = time of collision
        tab = table of collision data
        
        Output:
        child = name of child embryo"""

    maskt = tab['iap'] == pid
    maskp = tab['ilp'] == pid
    targ_colls = tab[maskt]
    proj_colls = tab[maskp]

    #pick the earliest collision for each one
    if len(targ_colls) >= 1:
        child_t = targ_colls[targ_colls['time']> tloss] #all collisions after the collision time
        if len(child_t) <= 0:
            child_tname = None
        else:
            child_tname = 'targ{0}_time{1}_c{2}'.format(child_t['iap'].iloc[0],child_t['time'].iloc[0],child_t['itype'].iloc[0])
    else:
        child_tname = None

    if len(proj_colls) >= 1:
        child_p = proj_colls[proj_colls['time'] > tloss]
        if len(child_p) <= 0:
            child_pname = None
        else:
            child_pname = 'proj{0}_time{1}_c{2}'.format(child_p['ilp'].iloc[0],child_p['time'].iloc[0],child_p['itype'].iloc[0])
    else:
        child_pname = None

    if all([child_tname,child_pname]):
        #then you have to compare them and pick the earliest one
        if child_p['time'].iloc[0] > child_t['time'].iloc[0]:
            child = child_tname
        elif child_p['time'].iloc[0] < child_t['time'].iloc[0]:
            child = child_pname
        elif child_p['time'].iloc[0] == child_t['time'].iloc[0]:
            print('Error: collisions happen at same time')
            print(child_p,child_t)
    elif child_pname is not None:
        child = child_pname
    elif child_tname is not None:
        child = child_tname
    else:
        #there are no children for this id
        child = None

    return(child)

def find_parents(pid,tloss,tab):
    """Find parents of a planet.
        Input:
        pid = planet id
        tloss = collision time
        tab = table of collision data
        
        Output:
        parent = name of parent"""
    
    #get collisions with this planet id
    maskt = tab['iap'] == pid
    maskp = tab['ilp'] == pid
    targ_colls = tab[maskt]
    proj_colls = tab[maskp]

    #pick the earliest collision for each one
    if len(targ_colls) >= 1:
        parent_t = targ_colls[targ_colls['time'] < tloss] #all collisions before this one
        lt = len(parent_t)
        if lt <= 0:
            parent_tname = None
        else:
            parent_tname = 'LR{0}_time{1}_c{2}'.format(parent_t['iap'].iloc[lt-1],parent_t['time'].iloc[lt-1],parent_t['itype'].iloc[lt-1])
    else:
        parent_tname = None

    if len(proj_colls) >= 1:
        parent_p = proj_colls[proj_colls['time'] < tloss]
        lp = len(parent_p)
        if lp <= 0:
            parent_pname = None
        else:
            parent_pname = 'SLR{0}_time{1}_c{2}'.format(parent_p['ilp'].iloc[lp-1],parent_p['time'].iloc[lp-1],parent_p['itype'].iloc[lp-1])
    else:
        parent_pname = None

    #pick which parent id is the most recent
    if all([parent_tname,parent_pname]):
        if parent_p['time'].iloc[lp-1] < parent_t['time'].iloc[lt-1]:
            #then LR is more recent, pick that as parent
            parent = [parent_tname]
        elif parent_t['time'].iloc[lt-1] < parent_p['time'].iloc[lp-1]:
            #SLR is more recent, pick that one
            parent = [parent_pname]
        elif parent_p['time'].iloc[lp-1] == parent_t['time'].iloc[lt-1]:
            print('Error: collisions happen at same time')
            print(parent_p,parent_LR)
    elif parent_tname is not None:
        parent = [parent_tname]
    elif parent_pname is not None:
        parent = [parent_pname]
    else:
        parent = [None]

    return(parent)

def add_collision(new,CollTree,i,collhistg_sorted):
    """Adds the bodies in one collision row into the colltree table with their parents and children
        Input:
        new = the collision being added
        CollTree = CollTree table for this planet
        i = current line of the table
        collhistg_sorted = collision info to search for parents and children
        
        Output:
        CollTree = updated CollTree table for this planet
        i = current line of table"""

    #start with the target
    ##########################################################################################
    name_t = 'targ{0}_time{1}_c{2}'.format(new['iap'],new['time'],new['itype'])

    #children name
    name_LR = 'LR{0}_time{1}_c{2}'.format(new['iap'],new['time'],new['itype'])
    if new['SLRMass'] > 0:
        name_SLR = 'SLR{0}_time{1}_c{2}'.format(new['ilp'],new['time'],new['itype'])
        children = [name_LR,name_SLR]
    else:
        children = [name_LR]

    if i > 0:
        #find parents
        targ_parents = find_parents(new['iap'],new['time'],collhistg_sorted)
        #or somehow assign parents when finding children?
    else:
        targ_parents = [None]

    CollTree.loc[i] = [new['time'],name_t,new['iap'],new['tmass'],(new['CMFt']+new['CMFt_min'])/2,
                       new['a'],new['itype'],children,targ_parents]
    #add to i
    i = i+1

    #do the projectile
    ###########################################################################################
    name_p = 'proj{0}_time{1}_c{2}'.format(new['ilp'],new['time'],new['itype'])

    if i > 1:
        #find parents
        proj_parents = find_parents(new['ilp'],new['time'],collhistg_sorted)
    else:
        proj_parents = [None]

    CollTree.loc[i] = [new['time'],name_p,new['ilp'],new['pmass'],(new['CMFp']+new['CMFp_min'])/2,new['a'],new['itype'],children,proj_parents]
    #add to i
    i = i+1

    #do the LR
    ###########################################################################################
    #find children
    child = find_children(new['iap'],new['time'],collhistg_sorted)

    #set parents
    parents = [name_t,name_p]

    CollTree.loc[i] =[new['time'],name_LR,new['iap'],new['LRMass'],(new['CMFLR']+new['CMFLR_min'])/2,new['a'],new['itype'],[child],parents]
    i = i+1 #add to i

    #check if there's a second largest remnant
    #########################################################################################
    if new['SLRMass'] > 0.:
        #see if there's a child for this one
        child_2 = find_children(new['ilp'],new['time'],collhistg_sorted)

        CollTree.loc[i] = [new['time'],name_SLR,new['ilp'],new['SLRMass'],(new['CMFSLR']+new['CMFSLR'])/2,new['a'],new['itype'],[child_2],parents]
        #add to i
        i = i+1

    return(CollTree,i)

def get_CollTree(d,p_id,collhists):
    """Make a collision tree table for one planet.
        
        Input:
        d = directory
        p_id = planet id
        collhists = collision history table"""

    #set up table
    CollTree = pd.DataFrame(columns=['time','body','id','mass','cmf','a','ctype','children','parents'])

    #mask collisions directory to get collision history for this planet
    maskd = collhists['dir'] == d
    maskp = collhists['pid'] == p_id
    collhist_giant = collhists[maskp&maskd]

    #sort table by time
    collhistg_sorted = collhist_giant.sort_values('time')

    i=0
    #iterate through collisions to add to CollTree
    for j in range(0,len(collhistg_sorted)):
        new = collhistg_sorted.iloc[j]
        CollTree, i = add_collision(new,CollTree,i,collhistg_sorted)

    return(CollTree)
