# remote_sensing
Extension to:
Simulation of satellite detection of unresolved subpixel fire for PC4262 Remote Sensing.

# files
1. `hotstuff.py`

  A helper `python` script that when imported, initializes common physical constants and imports multiple functions that compute, transfer and integrate radiances, designed for rapid computation as matrices via the `numpy` module. Written with grid simulation in mind.
  
 
2. `fire_funcs.py`

  Helper functions to do the simulations of fire and simulate detection algorithms.
  Updated: fire temperature evolution follows a S-curve with log decay
  
3. `fire_detection_plots.py`

  Run this script for a example of the simulation. Outputs: 1) MP4 of fire spread, 2) PNG of biomass grid used 3) MP4 of detection at saturation temperature, and 5) MP4 of detection using MODIS absolute and relative fire algorithms. Uses `matplotlib.animation` module to make the MP4s.
  
#

# Additional links

* https://code.earthengine.google.com/53303233dc0d27a44ed6cee64288a7ba?hideCode=true

Google Earth Engine script to obtain land surface temperature and vegetation data from MODIS.
