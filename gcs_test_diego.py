import sunpy
from sunpy.net import Fido, attrs as a
import astropy.units as u
from glob import glob
from constants import LASCO_PATH, AIA193_PATH, EUVI195_A_PATH, EUVI195_B_PATH, COR2A_PATH, COR2B_PATH, SAVE_PATH, \
    GCS_DELTA_T
from sunpy.time import parse_time
from pyGCSgui import runGCSgui
import pyGCS
import get_event_gcs
from astropy.io import fits
import numpy as np
from wrapper_eeggl import *


if __name__ == "__main__":
    frame = 8
    instrument = ""
    dir = "/data_local/GCS/2023_PyGCS_Raytrace/output/2011-02-15/"
    fecha = "2011-02-15"
    save_name = fecha+instrument+"_time_"+str(frame)+".txt"
    base_images, cme_images = wrapper_20110215(frame)

    #Cor2A, Cor2B, C2
    fnameA0, fnameB0, fnameL0 = base_images  #cor2a_0, cor2b_0, lasco_0
    fnameA1, fnameB1, fnameL1 = cme_images   #cor2a_1, cor2b_1, lasco_1
    
    myfitsA0 = fits.open(fnameA0)
    ima0 = myfitsA0[0].data
    hdra0 = myfitsA0[0].header
    myfitsA1 = fits.open(fnameA1)
    ima1 = myfitsA1[0].data
    hdra1 = myfitsA1[0].header

    myfitsB0 = fits.open(fnameB0)
    imb0 = myfitsB0[0].data
    hdrb0 = myfitsB0[0].header
    myfitsB1 = fits.open(fnameB1)
    imb1 = myfitsB1[0].data
    hdrb1 = myfitsB1[0].header

    myfitsL0 = fits.open(fnameL0)
    imL0 = myfitsL0[0].data
    hdrL0 = myfitsL0[0].header
    myfitsL1 = fits.open(fnameL1)
    imL1 = myfitsL1[0].data
    hdrL1 = myfitsL1[0].header

    headers = [hdra1,hdrb1, hdrL1]
    ims = [np.transpose(ima1 - ima0), np.transpose(imb1 - imb0), np.transpose((imL1 - imL0))]
    
    #Lasco-C2 needs to be rotated in heliographic latitude
    do_rotate_lat=[]
    for instr in range(len(headers)):
        if headers[instr]['TELESCOP'] == 'SOHO': 
            do_rotate_lat.append(True) 
        else:
            do_rotate_lat.append(False)

    satpos, plotranges = pyGCS.processHeaders(headers)
    sats = [["STEREO-A", "COR2"], ["STEREO-B", "COR2"], ["LASCO", "C2"]]
    date_obs = [hdra1['DATE-OBS'], hdrb1['DATE-OBS'], hdrL1['DATE-OBS']]
    filenames = [hdra1['filename'], hdrb1['filename'], hdrL1['filename']]
    #breakpoint()
    runGCSgui(ims, satpos, plotranges, sats,date_obs, save_name,do_rotate_lat, [5, 20, 30])

    # This is a test for running several events at once
    # start_time = parse_time("2011-06-05 07:00:00")
    # end_time = parse_time("2011-06-05 14:00:00")
    # get_event_gcs.run_event_gcs(start_time=start_time, end_time=end_time, event_id=52, cadence=8, do_low_corona=False)
