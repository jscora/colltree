#set constants
au = 1.495978707e8 #in km
munit = 5.0428958e31 #in grams
mearth = 5.972e24*1e3 #in g
msol = 1.989e30*1e3 #in g

def munit_to_mearth(x):
    munit = 5.0428958e31 #in grams
    """Function that converts code units to earth masses"""
    return(x*munit/mearth)

def read_comp(base_dir,pdir,fname_min='pl.mincorecompositions-nograzefa',fname_max='pl.maxcorecompositions-nograzefa'):
    """Reads in output files for one simulation"""

    try:
        comp = pd.read_csv(base_dir+pdir+fname_max,
                        names=('time','iinit','a','e','inc','mass','inew','mcore','mmant','mtot'),
                        delim_whitespace=True)
        comp_min = pd.read_csv(base_dir+pdir+fname_min,
                        names=('time','iinit','a','e','inc','mass','inew','mcore','mmant','mtot'),
                        delim_whitespace=True)
    except:
        print('cannot open files from: ',base_dir+pdir)
        
    #combine comp_min and comp so there is an average CMF, min CMF and max CMF
    if len(comp) != len(comp_min):
        raise Exception('output files are not the same length. Min file is {0} and max file is {1}'.format(len(comp_min),len(comp)))
    
    comp['cmf_max'] = comp['mcore']/comp['mtot']
    comp_min['cmf_min'] = comp_min['mcore']/comp_min['mtot']
    comp['cmf'] = (comp['cmf_max']+comp_min['cmf_min'])/2

    return(comp)

def read_follow(base_dir,pdir,fname_min='follow.mincorecollisions-nograzefa',fname_max='follow.maxcorecollisions-nograzefa'):
    
    try:
        coll = Table.read(base_dir+pdir+fname_max,format='ascii.no_header',
                      names=('time','a','iac','iap','tmass','CMFt','il','ilp','pmass','CMFp','itype','iLR',
                             'LRMass','CMFLR','iSLR','SLRMass','CMFSLR','inew','ideb','mdeb'))
        collmin = Table.read(base_dir+pdir+fname_min,format='ascii.no_header',
                      names=('time','a','iac','iap','tmass','CMFt','il','ilp','pmass','CMFp','itype','iLR',
                             'LRMass','CMFLR','iSLR','SLRMass','CMFSLR','inew','ideb','mdeb'))
    except:
        print('Could not open files in directory ',base_dir+pdir)

    if len(coll) != len(collmin):
        raise Exception('Error: collision tables not same length for directory: ',pdir)

    #combine min and max info into one table
    coll['CMFt_min'] = collmin['CMFt']
    coll['CMFp_min'] = collmin['CMFp']
    coll['CMFLR_min'] = collmin['CMFLR']
    coll['CMFSLR_min'] = collmin['CMFSLR']

    return(coll)
