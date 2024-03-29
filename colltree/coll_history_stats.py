import numpy as np
from astropy.table import vstack,Table,Column,setdiff,join
#from tqdm import tqdm
from tqdm import tqdm_notebook as tqdm
import astropy.units as u
import astropy.constants as c
import astropy.io
import colltree.util as util

#set constants
au = 1.495978707e8 #in km
munit = 5.0428958e31 #in grams
mearth = 5.972e24*1e3 #in g
msol = 1.989e30*1e3 #in g

def colls_per_planet(base_dir,collhist,good_dirs,sparam,minemb,pl_param,cname_max,cname_min):
    """A function that takes the collision history table and counts the total number of collisions that occurred in a planet's history.
        Inputs:
        base_dir = directory where data directories are
        collhist = collision history table
        good_dirs = table of directories and parameters to read data from
        sparam = types of collisions to count (i.e. small, giant)
        minemb = minimum planet mass. Options are 'embryo', where it takes the minimum embryo mass from continue.in, or the minimum planet mass you want.
        pl_param = if 'all', using final planets, if 'some', could be other bodies
        cname_max = name of composition file (if not provided, default is used)
        cname_min = name of composition file (if not provided, default is used)
        
        Output:
        numcoll_tot = table with total number of collisions per planet"""
    
    if sparam == 'giant':
        #set up table to tally number of each collision type
        numcoll_tot = Table(names=('dirs','pid','c1','c2','c3','c4','c5','c6','c7','c8','c9','tot_giant',
                                   'mass','a','e','inc','cmf_max','cmf_min','cmf_av','av_ecc',
                                   'dloss','av_init','slope'),
                    dtype=('U32','int64','int32','int32','int32','int32','int32','int32','int32','int32',
               'int64','int64','float64','float64','float32','float32','float64','float64',
               'float64','float64','U12','float64','U12'))
        colls = [1,2,3,4,5,6,7,8,9]
    
    elif sparam == 'small':
        #set up table
        numcoll_tot = Table(names=('dir','pid','tot_small','mass','a','e','inc','cmf_max','cmf_min','cmf_av',
                                   'reacc','supercat','other','m_reacc','m_super','m_other','mcore_min',
                                   'mcore_max','mdebris','av_ecc','dloss','av_inc','slope'),
                            dtype=('U32','int32','int32','float64','float64','float64','float64',
                                   'float64','float64','float64','int32','int32','int32','float64',
                                   'float64','float64','float64','float64','float64','float64','U12',
                                   'float64','U12'))

    new = collhist.group_by('dir')

    for key, group in zip(new.groups.keys, new.groups):
        print(key['dir'])
        p_mask = good_dirs['dir'] == key['dir']
        params = good_dirs[p_mask]
        if len(params['av_ecc']) > 1:
            print('error in masking')
            print(params)

        coll = group.group_by('pid')

        #read in pl file for that directory to match final planet info to pid
        comp = util.read_comp(base_dir,key['dir'],'astropy',cname_min,cname_max)

        mtiny = np.genfromtxt(base_dir+key['dir']+'/continue.in')
        
        #mask for planets above a certain mass
        if minemb == 'embryo':
            maskm = comp['mtot'] > mtiny[2]*munit/mearth
        else:
            maskm = comp['mtot'] >= minemb

        embryos = comp[maskm]
            
        if pl_param == 'all':
            #make table of final bodies for this run
            pl_group = embryos.group_by('time')
            lc = len(pl_group.groups)
            planets = pl_group.groups[lc-1]
        else:
            planets = embryos
        
        if len(planets)==0:
            raise Exception('no final planets found')

        for key, group in zip(coll.groups.keys, coll.groups):

            if pl_param == 'all':
                #get planet semi-major axis, cmf, and final mass
                maskpid = planets['iinit'] == key['pid']
                n_planet = planets[maskpid]
            elif pl_param == 'some':
                maskpid = planets['iinit'] == key['pid']
                if len(planets[maskpid]) > 1:
                    planets[maskpid].sort('time')
                    n_planet = planets[maskpid][-1]
                else:
                    n_planet = planets[maskpid]
            else:
                print('something wrong, pl_param =',pl_param)
                    
            if len(n_planet) < 1:
                #something is wrong, couldn't find planet to match
                print('no matching planet for: ',group['dir'][0],key['pid'])
                cmf = -1
                cmfm = -1
                cmf_av = -1
                mass = -1
                a = -1
                e = -1
                inc = -1
            else:
                cmfm = n_planet['cmf_max']
                cmf = n_planet['cmf_min']
                cmf_av = n_planet['cmf']
                mass = n_planet['mtot']
                a = n_planet['a']
                e = n_planet['e']
                inc = n_planet['inc']


            if sparam == 'giant' or sparam=='all':
                #total number of each type of collision
                ncolls = []
            
                for c in colls:
                    #calculate number of each type of collision that happened
                    maskc=group['itype']==c
                    l = len(group['itype'][maskc])
                    ncolls.append(l)
                    if len(group['itype']) == 1:
                        if group['itype'][0] == 0:
                            #if there is only one entry and it has an itype of 0, it had no collisions
                            totl_giant = 0
                        else:
                            totl_giant = len(group['itype'])
                    elif len(group['itype']) > 1:
                        totl_giant = len(group['itype'])
                
                numcoll_tot.add_row([group['dir'][0],key['pid'],ncolls[0],ncolls[1],ncolls[2],ncolls[3],ncolls[4],ncolls[5],ncolls[6],ncolls[7],ncolls[8],totl_giant,mass,a,e,inc,cmfm,cmf,cmf_av,params['av_ecc'],params['dloss'],params['av_inc'],params['slope']])

            elif sparam == 'small' or sparam=='all':
                if group['itype'][0] == 0:
                    tot = 0
                else:
                    tot = len(group['itype']) #total number of small collisions
                                                     
                counts = Counter(group['origin_type'])
                #get mass fraction of debris for each origin type
                origins = group.group_by('origin_type')
                #use aggregate to get total mass of debris
                massfrac = origins['origin_type','pmass'].groups.aggregate(np.sum)

                #assign massfraction of each origin type to a variable
                maskr = massfrac['origin_type'] == 'reacc'
                if len(massfrac[maskr]['pmass']) == 0:
                    m_reacc = 0.
                else:
                    m_reacc = massfrac[maskr]['pmass'][0]

                masks = massfrac['origin_type'] == 'supercat'
                if len(massfrac[masks]['pmass']) == 0:
                    m_super = 0.
                else:
                    m_super = massfrac[masks]['pmass'][0]

                masko = massfrac['origin_type'] == 'other'
                if len(massfrac[masko]['pmass']) == 0:
                    m_other = 0.
                else:
                    m_other = massfrac[masko]['pmass'][0]

                #get planet parameters to add to numcoll
                mcore_min = np.sum(group['CMFp_min']*group['pmass']*munit/mearth)
                mcore_max = np.sum(group['CMFp']*group['pmass']*munit/mearth)
                mtot = np.sum(group['pmass'])*munit/mearth

                #add to numcoll
                numcoll_tot.add_row([group['dir'][0],key['pid'],tot,mass,a,e,inc,cmfm,cmf,cmf_av,counts['reacc'],
                                     count['supercat'],counts['other'],m_reacc,m_super,m_other,mcore_min,mcore_max,mtot,
                                     params['av_ecc'],params['dloss'],params['av_inc'],params['slope']])


                                                                                                         
    return(numcoll_tot)


def get_numcoll(base_dir,dirtable,collhists,cparam,minemb,numcoll_name,pl_param,fwrite,cname_min='pl.mincorecompositions-nograzefa',cname_max='pl.maxcorecompositions-nograzefa'):
    """gets numcoll table for a specified set of runs
        
        Input:
        base_dir = base directory where all runs can be found
        dirtable = csv file with directories to use
        collhists = collhist file to use 
        cparam = what type of collisions to do (default = giant)
        minemb = minimum mass for planets. If you want all embryos, use 'embryo'
        numcoll_name = name of numcoll file
        pl_param = if 'all', then you give the code collhist file made with 'all'. Otherwise, use 'some'. 
        fwrite = if True, overwrites pre-existing files
        cname_min = name of composition file, default used if not given
        cname_max = name of composition file, default used if not given
        
        Output:
        numcoll = table of collision stats"""
    
    dirs = Table.read(base_dir+dirtable)
    
    numcoll = colls_per_planet(base_dir,collhists,dirs,cparam,minemb,pl_param,cname_min,cname_max)
    if fwrite:
        numcoll.write(base_dir+numcoll_name,format='ascii.csv',overwrite=True)

    return(numcoll)
