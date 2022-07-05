#set constants
au = 1.495978707e8 #in km
munit = 5.0428958e31 #in grams
mearth = 5.972e24*1e3 #in g
msol = 1.989e30*1e3 #in g

def munit_to_mearth(x):
    munit = 5.0428958e31 #in grams
    """Function that converts code units to earth masses"""
    return(x*munit/mearth)
