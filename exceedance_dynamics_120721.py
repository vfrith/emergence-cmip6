'''
26 Jun 2021
Edit 06 Jul 2021
Edit 12 Jul 2021 to add multiple models
US heatwave
For CanESM model - scenario data (rcp8.5)
 or MIROC, rcp8.5
Options - 
   Plot timeseries for location
   Plot histograms for location
   Plot dynamics of top n(10) events within period (y1, y2) - (SLP and TASMAX, absolute and rel to mean)Assessing the US heatwave - chances in future climates 
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
    for each in folder_list[:20]:
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





def load_MIROC_scenario_tas_meanfield(constr):
    '''
    Extract daily field for prescribed ens member, year, month, day
    ens: ensemble member
    '''
    folder_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MIROC/MIROC6/ssp585/*')
    # individual ensemble
    new_cube = iris.cube.CubeList([])
    for each in folder_list:
        list_files = glob.glob(each+'/day/tasmax/gn/latest/*')
        cubes = iris.load(list_files, constr)
        equalise_attributes(cubes) # to allow merge to one cube    
        new_cube.append(cubes.concatenate_cube())
    # merge ensembles
    for count, each in enumerate(new_cube):
        print(count)
        realization_coord = iris.coords.AuxCoord(count, 'realization')
        each.add_aux_coord(realization_coord)
    equalise_attributes(new_cube) # to allow merge to one cube    
    full_cube = new_cube.merge_cube()
    return full_cube.collapsed(['realization', 'time'], iris.analysis.MEAN)


def load_CanESM_scenario_tas_meanfield(constr):
    '''
    Extract daily field for prescribed ens member, year, month, day
    ens: ensemble member
    '''
    folder_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/CCCma/CanESM5/ssp585/*')
    # individual ensemble
    folder_list.remove(folder_list[21])
    new_cube = iris.cube.CubeList([])
    for each in folder_list:
        list_files = glob.glob(each+'/day/tasmax/gn/latest/*')
        cubes = iris.load(list_files, constr)
        equalise_attributes(cubes) # to allow merge to one cube    
        new_cube.append(cubes.concatenate_cube())
    # merge ensembles
    for count, each in enumerate(new_cube):
        print(count)
        realization_coord = iris.coords.AuxCoord(count, 'realization')
        each.add_aux_coord(realization_coord)
    equalise_attributes(new_cube) # to allow merge to one cube    
    full_cube = new_cube.merge_cube()
    return full_cube.collapsed(['realization', 'time'], iris.analysis.MEAN)

def load_MIROC_scenario_slp_meanfield(constr):
    '''
    Extract daily field for prescribed ens member, year, month, day
    ens: ensemble member
    '''
    folder_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MIROC/MIROC6/ssp585/*')
    # individual ensemble
    new_cube = iris.cube.CubeList([])
    for each in folder_list:
        list_files = glob.glob(each+'/day/psl/gn/latest/*')
        cubes = iris.load(list_files, constr)
        equalise_attributes(cubes) # to allow merge to one cube    
        new_cube.append(cubes.concatenate_cube())
    # merge ensembles
    for count, each in enumerate(new_cube):
        print(count)
        realization_coord = iris.coords.AuxCoord(count, 'realization')
        each.add_aux_coord(realization_coord)
    equalise_attributes(new_cube) # to allow merge to one cube    
    full_cube = new_cube.merge_cube()
    return full_cube.collapsed(['realization', 'time'], iris.analysis.MEAN)

def load_CanESM_scenario_slp_meanfield(constr):
    '''
    Extract daily field for prescribed ens member, year, month, day
    ens: ensemble member
    '''
    folder_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/CCCma/CanESM5/ssp585/*')
    # individual ensemble
    folder_list.remove(folder_list[21])
    new_cube = iris.cube.CubeList([])
    for each in folder_list:
        list_files = glob.glob(each+'/day/psl/gn/latest/*')
        cubes = iris.load(list_files, constr)
        equalise_attributes(cubes) # to allow merge to one cube    
        new_cube.append(cubes.concatenate_cube())
    # merge ensembles
    for count, each in enumerate(new_cube):
        print(count)
        realization_coord = iris.coords.AuxCoord(count, 'realization')
        each.add_aux_coord(realization_coord)
    equalise_attributes(new_cube) # to allow merge to one cube    
    full_cube = new_cube.merge_cube()
    return full_cube.collapsed(['realization', 'time'], iris.analysis.MEAN)


def load_MIROC_scenario_tasmax_field(ens, constr):
    '''
    Extract daily field for prescribed ens member, year, month, day
    ens: ensemble member
    '''
    folder_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MIROC/MIROC6/ssp585/*')
    # individual ensemble
    list_files = glob.glob(folder_list[ens]+'/day/tasmax/gn/latest/*')
    return iris.load_cube(list_files, constr)

def load_MIROC_scenario_slp_field(ens, constr):
    '''
    Extract daily field for prescribed ens member, year, month, day
    ens: ensemble member
    '''
    folder_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/MIROC/MIROC6/ssp585/*')
    # individual ensemble
    list_files = glob.glob(folder_list[ens]+'/day/psl/gn/latest/*')
    return iris.load_cube(list_files, constr)

def load_CanESM_scenario_tasmax_field(ens, constr):
    '''
    Extract daily field for prescribed ens member, year, month, day
    ens: ensemble member
    '''
    folder_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/CCCma/CanESM5/ssp585/*')
    folder_list.remove(folder_list[21])
    # individual ensemble
    list_files = glob.glob(folder_list[ens]+'/day/tasmax/gn/latest/*')
    return iris.load_cube(list_files, constr)

def load_CanESM_scenario_slp_field(ens, constr):
    '''
    Extract daily field for prescribed ens member, year, month, day
    ens: ensemble member
    '''
    folder_list = glob.glob('/bp1store/geog-tropical/data/CMIP6/ScenarioMIP/CCCma/CanESM5/ssp585/*')
    folder_list.remove(folder_list[21])
    # individual ensemble
    list_files = glob.glob(folder_list[ens]+'/day/psl/gn/latest/*')
    return iris.load_cube(list_files, constr)

def identify_n_extremes(cube, year, n):
    '''
    cube: cube of tasmax
    year: array of timeslice, eg [2000, 2020]
    n: returns top n events
    '''
    mod = mine.cube_to_array(mine.time_slice(cube, year[0], year[1])) 
    mod_indices = mod.argsort()[-n:]
    ens_max = []
    yr_max = []
    mon_max = []
    day_max = []
    val_max = []
    for x in mod_indices:
        [ens_temp, pos_temp] = np.where(cube.data == mod[x])
        val_max.append(mod[x])
        yr_temp = cube[ens_temp, pos_temp].coord('year').points[0]
        mon_temp = cube[ens_temp, pos_temp].coord('month_number').points[0]
        day_temp = cube[ens_temp, pos_temp].coord('day_of_month').points[0] 
        ens_max.append(ens_temp[0])
        yr_max.append(yr_temp)
        mon_max.append(mon_temp)
        day_max.append(day_temp)
    return ens_max, yr_max, mon_max, day_max, val_max


def plot_fields(tasmax, slp, figname):
    '''
    Plot two fields - tasmax and slp, side by side
    '''
    fig, axes = plt.subplots \
                (nrows=1, ncols=1, figsize=(10, 7.), dpi=80, num=None)
    brewer_cmap = mpl_cm.get_cmap('brewer_RdBu_11')
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    im = plt.pcolormesh(tasmax.coord('longitude').points, \
                        tasmax.coord('latitude').points, tasmax.data-273.15,  \
                        transform=ccrs.PlateCarree(), cmap=brewer_cmap)
    ax.coastlines()
    cbar_ax = fig.add_axes()
    cbar = plt.colorbar(im, cax=cbar_ax, orientation='horizontal')
    plt.savefig(figname+'_tasmax.png')
    plt.close()
    fig, axes = plt.subplots \
                (nrows=1, ncols=1, figsize=(10, 7.), dpi=80, num=None)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    im = plt.pcolormesh(slp.coord('longitude').points, \
                        slp.coord('latitude').points, slp.data,  \
                        transform=ccrs.PlateCarree(), cmap=brewer_cmap)
    ax.coastlines()
    cbar_ax = fig.add_axes()
    cbar = plt.colorbar(im, cax=cbar_ax, orientation='horizontal')
    plt.savefig(figname+'_slp.png')
    plt.close()
    return



def distributions_plot(data, title, filename, record):
    '''
    Plots histograms of distribution
    Saves plot in current folder
    data1: np.array, 1d, 
    title: string, on plot and filename
    record: observed maximum
    '''
    chance = ((len(data[data > record]))/len(data)) *120
    fig = plt.figure(figsize=(10., 5.), dpi=80, num=None)
    ax1 = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
    x = 46
    y = np.min(data)-2
    nbins = np.arange(0, x, x/100)
    ax1.hist(data, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon')
    plt.xticks([10, 20, 30, 40], color = 'dimgrey')
    plt.xlim([y, x])
    plt.axvline(record, color='r', lw=2)
    #plt.ylim([0, 0.012])
    ax1.spines['top'].set_color('grey')
    ax1.spines['right'].set_color('grey')
    ax1.spines['bottom'].set_color('grey')
    ax1.spines['left'].set_color('grey') 
    ax1.tick_params(axis='x', colors='dimgrey')
    plt.yticks([])
    plt.title(title+', Chance: '+str(round(chance, 2)), color = 'dimgrey')
    plt.tight_layout()
    plt.savefig(filename+'.png')
    plt.close()
    return


############
' Dictionary of individual locations for assessing' 
' Dictionary of individual locations for assessing' 
lookup_lon = {'PacNW':[245, 250], 'Lytton':[236, 241], 'WWA':[237, 241]}
lookup_lat = {'PacNW':[40, 45], 'Lytton':[48, 53], 'WWA':[45, 52]}



### VARIABLES
y1 = 2020
y2 = 2040
name = 'WWA'
record = 39.5 # max temp in degrees celcius, taken from WWA report
model = 'CanESM' # CanESM or MIROC

# what to do, 1 = do, 0 = don't
timeseries = 0
histograms = 0
dynamics = 1


###### MAIN #######
# Location constraint
loc_cons = iris.Constraint(latitude=lambda cell: lookup_lat[name][0] < cell < lookup_lat[name][1]) \
                    & iris.Constraint(longitude=lambda cell: lookup_lon[name][0] < cell < lookup_lon[name][1])
# Seasonal (JJA) constraint
season_cons = iris.Constraint(time=lambda cell: 6 <= cell.point.month <= 8)

# Load data
if model == 'CanESM': cube = load_CanESM_tasmax(loc_cons & season_cons)
if model == 'MIROC': cube = load_MIROC_tasmax(loc_cons & season_cons)

# Ensure cube 2-d (area mean)
if len(cube.coord('latitude').points)>1:
    cube = cube.collapsed(['latitude'], iris.analysis.MEAN)
if len(cube.coord('longitude').points)>1:
    cube = cube.collapsed(['longitude'], iris.analysis.MEAN)

# add time coords (could add to identify_n_extremes function)
iris.coord_categorisation.add_year(cube, 'time') # add coord year
iris.coord_categorisation.add_month_number(cube, 'time')
iris.coord_categorisation.add_day_of_month(cube, 'time')


if timeseries == True:
    # Plot timeseries to have an eyeball
    ens = cube.shape[0]
    for i in np.arange(ens):
        plt.plot(cube.coord('year').points, cube[i,:].data-273.15, '+')
    plt.axhline(record, lw=2)
    plt.savefig(name+'_ts_'+model+'.png')
    plt.close()

if histograms == True:
    # Plot histograms
    # 2020-2040
    mod = mine.cube_to_array(mine.time_slice(cube, 2020, 2040)) - 273.15
    mine.distributions_plot(mod, name+str(record)+' 2020-2040 '+model, name+'_20202040_hist_'+model, record)
    # 2050-2070
    mod = mine.cube_to_array(mine.time_slice(cube, 2050, 2070)) - 273.15
    mine.distributions_plot(mod, name+str(record)+' 2050-2070 '+model, name+'_20502070_hist_'+model, record)
    # 2080-2100
    mod = mine.cube_to_array(mine.time_slice(cube, 2080, 2100)) - 273.15
    mine.distributions_plot(mod, name+str(record)+' 2080-2100 '+model, name+'_20802100_hist_'+model, record)







if dynamics == True:
    # Identify top 10 events to see dynamics
    ens_max, yr_max, mon_max, day_max, val_max = identify_n_extremes(cube-273.15, [y1, y2], 5)
    # Load (absolute) fields of top events
    tas_fields = iris.cube.CubeList([]) 
    slp_fields = iris.cube.CubeList([])
    for i in np.arange(5):
        ens = ens_max[i]
        year_cons = iris.Constraint(time=lambda cell: cell.point.year == yr_max[i])
        month_cons = iris.Constraint(time=lambda cell: cell.point.month == mon_max[i])
        day_cons = iris.Constraint(time=lambda cell: cell.point.day == day_max[i])
        constr = (year_cons & month_cons & day_cons)
        if model == 'MIROC':
            tas_fields.append(load_MIROC_scenario_tasmax_field(ens, constr))
            slp_fields.append(load_MIROC_scenario_slp_field(ens, constr))
        if model == 'CanESM':
            tas_fields.append(load_CanESM_scenario_tasmax_field(ens, constr))
            slp_fields.append(load_CanESM_scenario_slp_field(ens, constr))
    # Plot absolute fields
    for i in np.arange(5):
        plot_fields(tas_fields[i], slp_fields[i], model+'_'+str(i))
    # Load mean field between y1 & y2
    year_cons = iris.Constraint(time=lambda cell: y1< cell.point.year < y2)
    mon_cons = iris.Constraint(time=lambda cell: 6 <= cell.point.month <= 8)
    constr = (year_cons & mon_cons)
    if model == 'MIROC':
        SLP_meanfield = load_MIROC_scenario_slp_meanfield(constr)
        TAS_meanfield = load_MIROC_scenario_tas_meanfield(constr)
    if model == 'CanESM':
        SLP_meanfield = load_CanESM_scenario_slp_meanfield(constr)
        TAS_meanfield = load_CanESM_scenario_tas_meanfield(constr)
    # Plot anomaly fields (rel to y1-y2 mean)
    for i in np.arange(5):
        plot_fields(tas_fields[i]-TAS_meanfield, slp_fields[i]-SLP_meanfield, model+'_anom_'+str(i))
