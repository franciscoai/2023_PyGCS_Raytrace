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

if __name__ == "__main__":

    lasco_0 = '/media/gehme/gehme/data/soho/lasco/level_1/c2/20110605/25374840.fts'
    cor2a_0 = '/media/gehme/gehme/data/stereo/secchi/L1/a/img/cor2/20110605/20110605_073900_14c2A.fts'
    cor2b_0 = '/media/gehme/gehme/data/stereo/secchi/L1/b/img/cor2/20110605/20110605_073900_14c2B.fts'
    lasco_1 = '/media/gehme/gehme/data/soho/lasco/level_1/c2/20110605/25374842.fts'
    cor2a_1 = '/media/gehme/gehme/data/stereo/secchi/L1/a/img/cor2/20110605/20110605_075400_14c2A.fts'
    cor2b_1 = '/media/gehme/gehme/data/stereo/secchi/L1/b/img/cor2/20110605/20110605_075400_14c2B.fts'

    fnameA0, fnameB0, fnameL0 = lasco_0, cor2a_0, cor2b_0
    fnameA1, fnameB1, fnameL1 = lasco_1, cor2a_1, cor2b_1

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

    headers = [hdrb1, hdrL1, hdra1]
    ims = [np.transpose(imb1 - imb0), np.transpose((imL1 - imL0)), np.transpose(ima1 - ima0)]
    
    
    satpos, plotranges = pyGCS.processHeaders(headers)
    sats = [["STEREO-A SECCHI", "COR2"], ["STEREO-B SECCHI", "COR2"], ["SOHO-LASCO", "C2"]]
    
    runGCSgui(ims, satpos, plotranges, sats, [5, 20, 30])

    # This is a test for running several events at once
    # start_time = parse_time("2011-06-05 07:00:00")
    # end_time = parse_time("2011-06-05 14:00:00")
    # get_event_gcs.run_event_gcs(start_time=start_time, end_time=end_time, event_id=52, cadence=8, do_low_corona=False)
