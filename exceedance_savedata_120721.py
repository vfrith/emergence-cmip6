'''
12 Jul 2021
Loading and saving iris cubes for a specific location, for CanESM and MIROC 
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


 
def load_CanESM_tasmax(constr):
    file_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/CCCma/CanESM5/ssp585/*/day/tasmax/gn/latest/*')
    cubes = iris.load(file_list, constr)
    cubes.remove(cubes[21]) # different dimensions so removed
    for count, each in enumerate(cubes):
        print(count)
        realization_coord = iris.coords.AuxCoord(count, 'realization')
        each.add_aux_coord(realization_coord)
    equalise_attributes(cubes) # to allow merge to one cube 
    return cubes.merge_cube()

def load_MIROC_tasmax(constr):
    folder_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MIROC/MIROC6/ssp585/*')
    # each ensemble
    new_cube = iris.cube.CubeList([])
    for each in folder_list[:25]:
        print(each)
        list_files = glob.glob(each+'/day/tasmax/gn/latest/*')
        cubes = iris.load(list_files, loc_cons & season_cons)
        equalise_attributes(cubes) # to allow merge to one cube    
        new_cube.append(cubes.concatenate_cube())
    # merge ensembles
    for count, each in enumerate(new_cube):
        print(count)
        realization_coord = iris.coords.AuxCoord(count, 'realization')
        each.add_aux_coord(realization_coord)
    equalise_attributes(new_cube) # to allow merge to one cube    
    return new_cube.merge_cube()

def load_CanESM_hist_tasmax(constr):
    file_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/CMIP/CCCma/CanESM5/historical/*/day/tasmax/gn/latest/*')
    cubes = iris.load(file_list, constr)
    #cubes.remove(cubes[21]) # different dimensions so removed
    for count, each in enumerate(cubes):
        print(count)
        realization_coord = iris.coords.AuxCoord(count, 'realization')
        each.add_aux_coord(realization_coord)
    equalise_attributes(cubes) # to allow merge to one cube 
    return cubes.merge_cube()

def load_MIROC_hist_tasmax(constr):
    folder_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/CMIP/MIROC/MIROC6/historical/*')
    # each ensemble
    new_cube = iris.cube.CubeList([])
    for each in folder_list[:25]:
        print(each)
        list_files = glob.glob(each+'/day/tasmax/gn/latest/*')
        cubes = iris.load(list_files, loc_cons & season_cons)
        equalise_attributes(cubes) # to allow merge to one cube    
        new_cube.append(cubes.concatenate_cube())
    # merge ensembles
    for count, each in enumerate(new_cube):
        print(count)
        realization_coord = iris.coords.AuxCoord(count, 'realization')
        each.add_aux_coord(realization_coord)
    equalise_attributes(new_cube) # to allow merge to one cube    
    return new_cube.merge_cube()

############
' Dictionary of individual locations for assessing' 
' Dictionary of individual locations for assessing' 
lookup_lon = {'PacNW':[245, 250], 'Lytton':[236, 241], 'WWA':[237, 241]}
lookup_lat = {'PacNW':[40, 45], 'Lytton':[48, 53], 'WWA':[45, 52]}



### VARIABLES
y1 = 2020
y2 = 2040
name = 'WWA'


###### MAIN #######
# Location constraint
loc_cons = iris.Constraint(latitude=lambda cell: lookup_lat[name][0] < cell < lookup_lat[name][1]) \
                    & iris.Constraint(longitude=lambda cell: lookup_lon[name][0] < cell < lookup_lon[name][1])
# Seasonal (JJA) constraint
season_cons = iris.Constraint(time=lambda cell: 6 <= cell.point.month <= 8)

##### CanESM
cube = load_CanESM_hist_tasmax(loc_cons & season_cons)

# Ensure cube 2-d (area mean)
if len(cube.coord('latitude').points)>1:
    cube = cube.collapsed(['latitude'], iris.analysis.MEAN)
if len(cube.coord('longitude').points)>1:
    cube = cube.collapsed(['longitude'], iris.analysis.MEAN)

# add time coords (could add to identify_n_extremes function)
iris.coord_categorisation.add_year(cube, 'time') # add coord year
iris.coord_categorisation.add_month_number(cube, 'time')
iris.coord_categorisation.add_day_of_month(cube, 'time')

iris.save(cube, 'tmp/WWA_hist_CanESM.nc')

##### MIROC
cube = load_MIROC_hist_tasmax(loc_cons & season_cons)

# Ensure cube 2-d (area mean)
if len(cube.coord('latitude').points)>1:
    cube = cube.collapsed(['latitude'], iris.analysis.MEAN)
if len(cube.coord('longitude').points)>1:
    cube = cube.collapsed(['longitude'], iris.analysis.MEAN)

# add time coords (could add to identify_n_extremes function)
iris.coord_categorisation.add_year(cube, 'time') # add coord year
iris.coord_categorisation.add_month_number(cube, 'time')
iris.coord_categorisation.add_day_of_month(cube, 'time')

iris.save(cube, 'tmp/WWA_hist_MIROC.nc')



