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

### VARIABLES
y1 = 1980
y2 = 2000


obs = iris.load_cube('tmp/WWA_ERA5_19792020.nc')
obs = obs.collapsed(['longitude', 'latitude'], iris.analysis.MEAN)
iris.coord_categorisation.add_year(obs, 'time') # add coord year
obs = mine.time_slice(obs, y1, y2) - 273.15

mod1 = iris.load_cube('tmp/WWA_hist_CanESM.nc')
mod1 = mine.cube_to_array(mine.time_slice(mod1, y1, y2)) - 273.15

mod2 = iris.load_cube('tmp/WWA_hist_MIROC.nc')
mod2 = mine.cube_to_array(mine.time_slice(mod2, y1, y2)) - 273.15

fig = plt.figure(figsize=(10., 5.), dpi=80, num=None)
ax1 = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
x = 46
y = np.min(obs.data)-2
nbins = np.arange(0, x, x/100)
ax1.hist(mod1, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon')
ax1.hist(obs.data, bins=nbins, normed=True, histtype='stepfilled', \
             color='grey', alpha=0.5)
plt.xticks([10, 20, 30, 40], color = 'dimgrey')
plt.xlim([y, x])
#plt.ylim([0, 0.012])
ax1.spines['top'].set_color('grey')
ax1.spines['right'].set_color('grey')
ax1.spines['bottom'].set_color('grey')
ax1.spines['left'].set_color('grey') 
ax1.tick_params(axis='x', colors='dimgrey')
plt.yticks([])
plt.tight_layout()
plt.savefig('bias_canesm.png')
plt.close()

fig = plt.figure(figsize=(10., 5.), dpi=80, num=None)
ax1 = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
x = 46
y = np.min(obs.data)-2
nbins = np.arange(0, x, x/100)
ax1.hist(mod2, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon')
ax1.hist(obs.data, bins=nbins, normed=True, histtype='stepfilled', \
             color='grey', alpha=0.5)
plt.xticks([10, 20, 30, 40], color = 'dimgrey')
plt.xlim([y, x])
#plt.ylim([0, 0.012])
ax1.spines['top'].set_color('grey')
ax1.spines['right'].set_color('grey')
ax1.spines['bottom'].set_color('grey')
ax1.spines['left'].set_color('grey') 
ax1.tick_params(axis='x', colors='dimgrey')
plt.yticks([])
plt.tight_layout()
plt.savefig('bias_miroc.png')
plt.close()

