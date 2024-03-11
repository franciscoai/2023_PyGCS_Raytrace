#!/usr/bin/env python
# coding: utf-8

# ## Librerias
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import pyGCS
from rtraytracewcs import rtraytracewcs
from get_corona_gcs_ml_diego import get_corona
import numpy as np
import datetime
import matplotlib.pyplot as plt
from astropy.io import fits
import sunpy
import sunpy.map
import pandas as pd
#from sunpy.sun.constants import radius as _RSUN
#from ext_libs.rebin import rebin
from coord_transformation import deg2px, center_rSun_pixel, pnt2arr
import scipy
from wrapper_eeggl import *
from read_GCS_variables_txt import read_GCS_variables_txt

def save_png(array, ofile=None, range=None):
    '''
    pltos array to an image in ofile without borders axes or anything else
    ofile: if not give only the image object is generated and returned
    range: defines color scale limits [vmin,vmax]
    '''    
    fig = plt.figure(figsize=(4,4), facecolor='white')
    if range is not None:
        vmin=range[0]
        vmax=range[1]
    else:
        vmin=None
        vmax=None
    plt.imshow(array, origin='lower', cmap='gray', vmin=vmin, vmax=vmax, interpolation='none')#, aspect='auto')#,extent=plotranges[sat])
    plt.axis('off')         
    if ofile is not None:
        fig.savefig(ofile, facecolor='white', bbox_inches='tight', pad_inches=0)
        plt.close(fig)
        return 1
    else:
        return fig


######Main
# CONSTANTS
#Dir where the output images are saved
OPATH = os.path.dirname(os.path.realpath(__file__)) + '/output'
#n_sat = 3 #number of satellites to  use [Cor2 A, Cor2 B, Lasco C2]

# Syntethic image options
imsize=np.array([512, 512], dtype='int32') # output image size
otype="png" # set the ouput file type: 'png' or 'fits'
otype="fits" 
im_range=2. # range of the color scale of the output final syntethyc png image in std dev around the mean
add_occ_to_btot = True # if True the occulter is added to the final btot image
level_occ=0. #mean level of the occulter relative to the background level
add_back_to_btot = False # if True the background corona is added to the final btot image
#breakpoint()
os.makedirs(OPATH, exist_ok=True)

back_corona=[]
headers=[]
#size_occ=[]
#size_occ_ext=[]
#satpos_all=[]
#plotranges_all=[]
n_sat = 3
for sat in range(n_sat):
    a,b,c,d=get_corona(sat,imsize=imsize)
    back_corona.append(a)
#    headers.append(b)
#    size_occ.append(c)
#    size_occ_ext.append(d)


frame =9
instrument = "cor2a"
fecha = "2011-02-15"
save_name = fecha+instrument+"_time_"+str(frame)+".txt"
base_images, cme_images = wrapper_20110215(frame)

for img in cme_images:
    oimg = fits.open(img)[0].data
    h0   = fits.getheader(img)
    headers.append(h0)
#Cor2A, Cor2B, C2
#fnameA0, fnameB0, fnameL0 = base_images  
#fnameA1, fnameB1, fnameL1 = cme_images 
size_occ     = [2.9, 4.0, 2.]
size_occ_ext = [16, 16, 6.]

# Get the location of sats and gcs:
satpos, plotranges = pyGCS.processHeaders(headers)
print(f'Saving images for GCS')


# GCS parameters 
# level_cme: CME intensity level relative to the mean background corona
par_names = ['CMElon', 'CMElat', 'CMEtilt', 'height', 'k','ang'] # par names
par_units = ['deg', 'deg', 'deg', 'Rsun','','deg'] # par units
#gcs_par = [45,10,50,10,0.35, 30] # min-max ranges of each parameter in par_names
gcs_par = read_GCS_variables_txt(save_name)
#breakpoint()
n_sat = len(base_images)
#breakpoint()
for sat in range(n_sat):
    #defining ranges and radius of the occulter
    x = np.linspace(plotranges[sat][0], plotranges[sat][1], num=imsize[0])
    y = np.linspace(plotranges[sat][2], plotranges[sat][3], num=imsize[1])
    xx, yy = np.meshgrid(x, y)
    x_cS, y_cS = center_rSun_pixel(headers, plotranges, sat)  
    r = np.sqrt((xx - x_cS)**2 + (yy - y_cS)**2)
    
    #Total intensity (Btot) figure from raytrace:               
    btot = rtraytracewcs(headers[sat], gcs_par[0], gcs_par[1], gcs_par[2], gcs_par[3], gcs_par[4], gcs_par[5],
                        imsize=imsize, occrad=size_occ[sat], in_sig=0.5, out_sig=0.25, nel=1e5)
    
    #mask for occulter
    arr = np.zeros(xx.shape)
    arr[r <= size_occ[sat]] = 1    
    arr[r >= size_occ_ext[sat]] = 1 

    #background corona
    back = back_corona[sat]
    level_back = np.mean(back)               

    #adds background
    if add_back_to_btot:
        btot+=back

    #adds occulter  
    if add_occ_to_btot:   
        noise =0        
        noise_ext=0 
        btot[r <= size_occ[sat]]     = level_occ*level_back + noise
        btot[r >= size_occ_ext[sat]] = level_occ*level_back + noise_ext
    #cme = fits.PrimaryHDU(btot)
    #cme.writeto(OPATH +'/sat{}_btot.fits'.format(sat+1), overwrite=True)
    
    
    if otype=="fits":
        cme = fits.PrimaryHDU(data=btot,header=headers[sat])
        cme.writeto(OPATH +'/{:05.2f}_{:05.2f}_{:05.1f}_{:05.2f}_{:05.2f}_{:05.2f}_sat{}_date{}_frame{}_btot.fits'.format(
            gcs_par[0], gcs_par[1], gcs_par[2], gcs_par[3], gcs_par[4], gcs_par[5], sat+1,fecha,frame), overwrite=True)
    elif otype =="png":
    #full image
        m = np.mean(btot)
        sd = np.std(btot)
        vmin=m-im_range*sd
        vmax=m+im_range*sd
        ofile = OPATH +'/{:08.3f}_{:08.3f}_{:08.3f}_{:08.3f}_{:08.3f}_{:08.3f}_sat{}_btot.png'.format(
            gcs_par[0], gcs_par[1], gcs_par[2], gcs_par[3], gcs_par[4], gcs_par[5], sat+1)
        fig=save_png(btot,ofile=ofile, range=[vmin, vmax])

###check this.
#                clouds = pyGCS.getGCS(df['CMElon'][row], df['CMElat'][row], df['CMEtilt'][row], df['height'][row], df['k'][row], df['ang'][row], satpos)                
#            x = clouds[sat, :, 1]
#            y = clouds[0, :, 2]         
#            arr_cloud=pnt2arr(x,y,plotranges,imsize, sat)
#            ofile = folder +'/{:08.3f}_{:08.3f}_{:08.3f}_{:08.3f}_{:08.3f}_{:08.3f}_sat{}_mesh.png'.format(
#                df['CMElon'][row], df['CMElat'][row], df['CMEtilt'][row], df['height'][row], df['k'][row], df['ang'][row], sat+1)
#            fig=save_png(arr_cloud,ofile=ofile, range=[0, 1])  



    # overplot  GCS mesh to cme figure
    #clouds = pyGCS.getGCS(gcs_par[0], gcs_par[1], gcs_par[2], gcs_par[3], gcs_par[4], gcs_par[5], satpos)                
    #x = clouds[sat, :, 1]
    #y = clouds[0, :, 2]    
    #arr_cloud=pnt2arr(x,y,plotranges,imsize, sat)
    #ofile = OPATH +'/{:08.3f}_{:08.3f}_{:08.3f}_{:08.3f}_{:08.3f}_{:08.3f}_sat{}_mesh.png'.format(
    #    gcs_par[0], gcs_par[1], gcs_par[2], gcs_par[3], gcs_par[4], gcs_par[5], sat+1)
    #fig=save_png(arr_cloud,ofile=ofile, range=[0, 1]) 
