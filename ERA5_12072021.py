'''
12 Jul 2021 
Loading ERA-5 data 
@vikki.thompson
'''

# Load neccessary libraries
import iris
import iris.coord_categorisation as icc
from iris.coord_categorisation import add_season_membership
import numpy as np
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy as cart
import glob
import matplotlib.cm as mpl_cm
import sys
sys.path.append('/home/hh21501/unseen') 
import core_functions as mine
from iris.experimental.equalise_cubes import equalise_attributes
import scipy.stats as sps

def load_ERA5_tasmax(constr):
    file_list = glob.glob('/bp1store/geog-tropical/data/ERA-5/day/tasmax/*')
    cubes = iris.load(file_list, constr)
    equalise_attributes(cubes) # to allow merge to one cube 
    return cubes.concatenate_cube()


############
' Dictionary of individual locations for assessing' 
lookup_lon = {'PacNW':[245, 250], 'Lytton':[236, 241], 'WWA':[237, 241]}
lookup_lat = {'PacNW':[40, 45], 'Lytton':[48, 53], 'WWA':[45, 52]}



### VARIABLES
y1 = 1979
y2 = 2020
name = 'WWA'


### CONSTRAINTS
# Location constraint
loc_cons = iris.Constraint(latitude=lambda cell: lookup_lat[name][0] < cell < lookup_lat[name][1]) \
                    & iris.Constraint(longitude=lambda cell: lookup_lon[name][0] < cell < lookup_lon[name][1])
# Seasonal (JJA) constraint
season_cons = iris.Constraint(time=lambda cell: 6 <= cell.point.month <= 8)
# Years constraint
year_cons = iris.Constraint(time=lambda cell: y1 <= cell.point.year <= y2)

obs = load_ERA_tasmax(loc_cons & season_cons & year_cons)

iris.save(obs, 'tmp/WWA_ERA5_19792020.nc')
