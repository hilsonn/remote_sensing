# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 22:16:01 2022

Running simulations and plotting

@author: hills
"""
import os
os.chdir(__file__[:-23]) # put the things in the same directory to make sure no weirdness with the importing
from fire_funcs import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter

vid_writer = FFMpegWriter(fps=10)
gif_writer = PillowWriter(fps=10)
### make simulation
print('Starting simulation...')
frames, firefrac_frames = start_sim()
print('Simulation complete\n')

## make frames for fire grid

print("Making frames for fire grid")
grid_frames = []
for frame in frames:
    pixel_dict = {}
    for i in range(dims[0]):
        for j in range(dims[1]):
            pixel_dict.update([[
                (i+1,j+1), frame[i*100:(i+1)*100, j*100:(j+1)*100]
                ]])
    grid_frames.append(pixel_dict)

### Start the fire and plot
fig,ax = plt.subplots(nrows=dims[0], ncols=dims[1],figsize=(12,12))
fig.subplots_adjust(hspace=0.01, wspace=0.01)
im = np.empty_like(ax)

for key in grid_frames[0]:
    pixel_frame = grid_frames[0][key]

    # using conditionals in order to account for pixel being 1x1, 1xN, or Nx1

    if dims[0]==1 and dims[1]==1:
        im = ax.imshow(pixel_frame, cmap = 'coolwarm', vmin = 30, vmax = 1200)
        ax.axis("off")
    elif dims[0]==1:
        im[key[1]-1] = ax[key[1]-1].imshow(pixel_frame, cmap = 'coolwarm', vmin = 30, vmax = 1200)
        [i.axis("off") for i in ax.ravel()]
    elif dims[1]==1:
        im[key[0]-1] = ax[key[0]-1].imshow(pixel_frame, cmap = 'coolwarm', vmin = 30, vmax = 1200)
        [i.axis("off") for i in ax.ravel()]
    else:
        im[key[0]-1, key[1]-1] = ax[key[0]-1, key[1]-1].imshow(pixel_frame, cmap = 'coolwarm', vmin = 30, vmax = 1200)
        [i.axis("off") for i in ax.ravel()]

def animate(i):
    fig.suptitle(f'$\\Delta t$ = {i}, mean temp = {round(np.mean(frames[i]),5)}K\nFire fraction = {firefrac_frames[i]:.5f}')
    for pixel in grid_frames[i]:
        pixel_frame = grid_frames[i][pixel]
        
        if dims[0]==1 and dims[1]==1:
            im.set_array(pixel_frame)
        elif dims[0]==1:
            im[pixel[1]-1].set_array(pixel_frame)
        elif dims[1]==1:
            im[pixel[0]-1].set_array(pixel_frame)
        else:
            im[pixel[0]-1, pixel[1]-1].set_array(pixel_frame)
    return [im]

anim = FuncAnimation(fig, animate, frames=N_frames, repeat=False)
if input("Save fire spread?(Y/N)\n").upper()=="Y":
    file_name = input('File name? ')
    print("saving...")
    anim.save(file_name+'.mp4', writer=vid_writer, dpi=100)
    print("saved\n")

if input("Save biomass image?(Y/N)\n").upper()=="Y":
    file_name = input('File name? ')
    print("saving...")
    fig = plt.figure()
    plt.title("Biomass grid")
    plt.imshow(bmass_grid, cmap = 'Greens')
    plt.axis('off')
    plt.savefig(file_name+".png", dpi=100)
    print("saved\n")

### Detections
# split up detected_list into the pixels
print("Making frames for pixel temperature detected by sensor")
detected_grid_list = []
for grid_frame in grid_frames:
    detected_pixel_dict = {}
    for key in grid_frame:
        pixel_frames = grid_frame[key]
        detected_list = make_detected_temp(pixel_frames)
        detected_pixel_dict.update([[
            key, detected_list
            ]])
    detected_grid_list.append(detected_pixel_dict) # each dict is one frame
print("done\n")

### Pure saturation
# make saturated lists
print("Making frames for 4um and 11um pixel saturation")
fire_detected_frames_4um = []
fire_detected_frames_11um = []

for frame in detected_grid_list:
    fire_pixels_frame_4um = np.full(dims, 0)
    fire_pixels_frame_11um = np.full(dims, 0)
    
    for pixel in frame:
        pixel_val = frame[pixel]
        fire_pixels_frame_4um[pixel[0]-1, pixel[1]-1] = pixel_val[0]>331
        fire_pixels_frame_11um[pixel[0]-1, pixel[1]-1] = pixel_val[1]>400

    fire_detected_frames_4um.append(fire_pixels_frame_4um)
    fire_detected_frames_11um.append(fire_pixels_frame_11um)

### initialise animation
detector_fig, detector_ax = plt.subplots(1,2,figsize=(8,6))
[i.axis("off") for i in detector_ax.ravel()]
detector_fig.suptitle(f'$\\Delta t$ = {i}, Does it saturate?\nFire fraction = {firefrac_frames[i]:.5f}')
ax_4um, ax_11um = detector_ax

im_detector_4um = ax_4um.imshow(fire_detected_frames_4um[0], cmap = 'coolwarm', vmin=-1,vmax=1.5)
ax_4um.set_title('4 microns',fontsize=16)

im_detector_11um = ax_11um.imshow(fire_detected_frames_4um[0], cmap = 'coolwarm', vmin=-1,vmax=1.5)
ax_11um.set_title('11 microns',fontsize=16)

def animate_saturator(i):
    detector_fig.suptitle(f'$\\Delta t$ = {i}, Does it saturate?\nFire fraction = {firefrac_frames[i]:.5f}')
    im_detector_4um.set_array(-1+2*fire_detected_frames_4um[i])
    im_detector_11um.set_array(-1+2*fire_detected_frames_11um[i])
    im_detector = [im_detector_4um, im_detector_11um]
    return [im_detector]

### start animation
anim_saturator = FuncAnimation(detector_fig, animate_saturator, frames=N_frames, repeat=False)
if input("Save saturation figure?(Y/N)\n").upper()=="Y":
    file_name = input('Name? ')
    print("saving...")
    anim_saturator.save(file_name+'.mp4', writer=vid_writer, dpi=100)
    print("saved\n")

### MODIS absolute fire sub-algo
print("Making frames for MODIS absolute fire sub-algorithm")
modis_abs_detected_frames = []
for frame in detected_grid_list:
    modis_abs_pixel_frame = np.full(dims,0)

    for pixel in frame:
        pixel_val = frame[pixel] # each pixel_val holds the T4 and T11 detected brightness temps
        modis_abs_pixel_frame[pixel[0]-1, pixel[1]-1] = int(modis_abs_detect(*pixel_val)) # feed the pixel_val into modis absolute fire algo

    modis_abs_detected_frames.append(modis_abs_pixel_frame)

### initialise animation
detector_fig = plt.figure()
plt.axis('off')

im_detector = plt.imshow(modis_abs_detected_frames[0], cmap = 'coolwarm',vmin=-1,vmax=1.5)
def animate_modis_abs(i):
    plt.title(f'$\\Delta t$ = {i}, Does it detect?\nFire fraction = {firefrac_frames[i]:.5f}')

    im_detector.set_array(-1+2*modis_abs_detected_frames[i])
    return [im_detector]
### start animation
anim_modis_abs = FuncAnimation(detector_fig, animate_modis_abs, frames=N_frames, repeat=False)

if input("Save MODIS absolute fire detected?(Y/N)\n").upper()=="Y":
    file_name = input('Name? ')
    print("saving...")
    anim_modis_abs.save(file_name+'.mp4', writer=vid_writer, dpi=200)
    print("saved\n")

### MODIS relative fire sub-algo
if dims != (1,1):
    print("Making frames for MODIS relative fire algorithm")
    modis_rel_detected_frames = []
    for frame in detected_grid_list:
        all_detected_vals = list(frame.values())
        modis_rel_pixel_frame = np.full(dims,0)
        for pixel in frame:
            pixel_val = frame[pixel]
            modis_rel_pixel_frame[pixel[0]-1, pixel[1]-1] = int(modis_rel_detect(*pixel_val, all_detected_vals))
        modis_rel_detected_frames.append(modis_rel_pixel_frame)

    ### initialise animation
    detector_fig = plt.figure()
    plt.axis('off')

    im_detector = plt.imshow(modis_rel_detected_frames[0], cmap = 'coolwarm',vmin=-1,vmax=1.5)
    def animate_animate_modis_rel(i):
        plt.title(f'$\\Delta t$ = {i}, Does it detect?\nFire fraction = {firefrac_frames[i]:.5f}')
        im_detector.set_array(-1+2*modis_rel_detected_frames[i])
        return [im_detector]
    ### start animation
    anim_modis_rel = FuncAnimation(detector_fig, animate_animate_modis_rel, frames=N_frames, repeat=False)

    if input("Save MODIS relative fire detected?(Y/N)\n").upper()=="Y":
        file_name = input('Name? ')
        print("saving...")
        anim_modis_rel.save(file_name+'.mp4', writer=vid_writer, dpi=200)
        print("saved\n")

input('Done, press enter to exit')
