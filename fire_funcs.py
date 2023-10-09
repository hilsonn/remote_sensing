# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 22:14:30 2022

Functions for fire simulations and detection

@author: hills
"""

import numpy as np
from hotstuff import totalfinder

# initialise global variables
N_frames = int(input('Number of frames for simulation? '))
dims = (int(input('Length in pixels? ')),int(input('Height in pixels? ')))

temp_bg = 32.123
temp_grid = np.full((dims[0]*100,dims[1]*100),temp_bg)
bmass_grid = np.random.normal(0.9,0.4,(dims[0]*100,dims[1]*100))
# ensure bmass is between 0 and 1
bmass_grid[bmass_grid<=0] = 0
bmass_grid[bmass_grid>=1] = 1
# check_grid stores how long a fire has been burning at the coord
check_grid = np.zeros((dims[0]*100,dims[1]*100))

def temp_curve(temp_0, tpeak, biomass):
    """ Returns the temperature evolution curve for all time, needs params to control peak temp and how long it burns for """
    time = np.linspace(0,20,101)

    y = tpeak/(1 + 2*np.exp(-time)) / (1 + (1-biomass)**2/1.2*time)
    y[y<temp_bg] = temp_bg # set min temp to bg temperature
    return y

# note for future: improve the fire simulation algo

def temp_evolution(time, temp_0, tpeak, biomass):
    """ Increments timestep from time = 0 by time """
    temp_ls = temp_curve(temp_0, tpeak, biomass)
#     print(temp_ls[time])
    return temp_ls[time]

# at every timestep, check the array for temp values. make it a binomial process with chance to spread proportional to T of a location

def spread_fire(temp_self, k, coords):
    """ Given a temp_self, determines probability of spreading and spreads the fire to other grids """
    prob_spread = k*(temp_self - temp_bg)/500
    x_coord, y_coord = coords
    for i,j in [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]:
        if x_coord+i< 0 or y_coord+j< 0 or x_coord+i>= check_grid.shape[0] or y_coord+j>= check_grid.shape[1]:
                continue
        if np.random.rand() < prob_spread:
            if check_grid[x_coord+i, y_coord+j]==0:
                temp_grid[x_coord+i, y_coord+j]=500
                check_grid[x_coord+i, y_coord+j]=1
    return temp_grid

def start_sim(temp_grid = temp_grid):
  """ Starts simulation """
  temp_grid[dims[0]*100//2,dims[1]*100//2] = 500 #seed center with fire, temperature = 500 celsius

  frames = [] # make frames for the gif
  firefrac_frames = [] # list of fire frac for each frame
  for _ in range(N_frames):
      old_temp_grid = np.copy(temp_grid)
      for coord, item in np.ndenumerate(old_temp_grid):
          k = bmass_grid[coord]
          temp_grid = spread_fire(item, k, coord)
          if check_grid[coord]!=0:
              temp_grid[coord] = temp_evolution(int(check_grid[coord]-1), 500, 1100, k)
              check_grid[coord]+=1
      firefrac_frames.append(np.sum(check_grid!=0)/np.size(check_grid))
      frames.append(temp_grid.copy()+273) # convert to Kelvins
  return frames, firefrac_frames

# for the detection part
def modis_abs_detect(T4, T11):
    """
    Inputs
        T4: Brightness temp at 4 micron band
        T11: Brightness temp at 11 micron band
    Outputs
        bool What the MODIS absolute fire algo would return
    """
    if T11>360: return True
    elif T4>310 and T4-T11>10: return True
    else: return False

def modis_rel_detect(T4, T11, grid):
    """
    Inputs
        grid: M by N array of pixels containing lists of length two representing T4 and T11 at each pixel
    Outputs
        bool What the MODIS relative fire algo would return
    """
    all_detected_vals = np.reshape(grid,(-1,2))
    T4_mean = np.mean(all_detected_vals[:,0])
    T4_stddev = np.std(all_detected_vals[:,0], ddof=1)
    T4_T11_diff_median = np.median(all_detected_vals[:,0]-all_detected_vals[:,1])
    T4_T11_diff_stddev = np.std(all_detected_vals[:,0]-all_detected_vals[:,1], ddof=1)

    if modis_abs_detect(T4, T11):
        return True
    elif (T4 > T4_mean + 3*T4_stddev) and (T4 - T11 > T4_T11_diff_median + 3*T4_T11_diff_stddev):
        return True
    else:
        return False
    
def make_detected_temp(frame):
  """ Find the temperature of the different bands detected by the satellite """
  return totalfinder(frame)
