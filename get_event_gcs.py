from glob import glob
from typing import List, Tuple

import astropy.units as apu
import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta
from sunpy.time import parse_time

import pyGCS
from constants import LASCO_PATH, AIA193_PATH, EUVI195_A_PATH, EUVI195_B_PATH, COR2A_PATH, COR2B_PATH, SAVE_PATH, \
    GCS_DELTA_T


def run_event_gcs(start_time: Time, end_time: Time, event_id: int, cadence: int = 4, do_low_corona: bool = False
                  ) -> List[Tuple]:

    """
        run_event_gcs Esta función toma la lista de datos para un evento, los ordena y
             llama a do_gcs para computar cgs con cierta cadencia.

        start_time: Event start time (at photosphere?)
        end_time: Event end time (in Lasco c2?)
        event_id: event id to give name to the final output
        cadence: Amount of gcs per hour
        do_low_corona (default value = False):
                          - False: Does gcs with lasco_c2, cor2_a & cor2_b
                          - True: Does gcs with aia, euvi_a & euvi_b
    """

    # Let's assume for now all events happen in the same day and don´t last till the other
    # or start the previous day, to simplify the code.
    # Later we'll fix this checking start and end time

    event_date = start_time.to_value(format='iso', subfmt='date')
    event_date = event_date.replace('-', '')

    event_save_path = SAVE_PATH + event_date + '/'
    event_lasco_path = LASCO_PATH + event_date + '/c2/'
    event_aia_path = AIA193_PATH + event_date + '/preped/'
    event_euvia_path = EUVI195_A_PATH + event_date + '/preped/'
    event_euvib_path = EUVI195_B_PATH + event_date + '/preped/'
    event_cor2a_path = COR2A_PATH + event_date + '/preped/'
    event_cor2b_path = COR2B_PATH + event_date + '/preped/'

    # Select files needed from files in folder
    lasco_files = event_file_select(files=(event_lasco_path + '*.fts'), start_time=start_time, end_time=end_time)
    aia_files = event_file_select(files=glob(event_aia_path + '*.fits'), start_time=start_time, end_time=end_time)
    euvia_files = event_file_select(files=glob(event_euvia_path + '*fts'), start_time=start_time, end_time=end_time)
    euvib_files = event_file_select(files=glob(event_euvib_path + '*fts'), start_time=start_time, end_time=end_time)
    cor2a_files = event_file_select(files=glob(event_cor2a_path + '*fts'), start_time=start_time, end_time=end_time)
    cor2b_files = event_file_select(files=glob(event_cor2b_path + '*fts'), start_time=start_time, end_time=end_time)

    # Check which instruments to run gcs with
    if do_low_corona:
        # Run gcs with aia & euvi files
        gcs_data = data_sort(aia_files, euvia_files, euvib_files, start_time, end_time, cadence)
    else:
        # Run gcs with lasco & cor files
        gcs_data = data_sort(lasco_files, cor2a_files, cor2b_files, start_time, end_time, cadence)

    event_gcs_out = []
    # Run Gcs for each set of instruments date in gcs_data
    for i in range(0, len(gcs_data) - 2):
        event_snapshot_0 = gcs_data[i]
        event_snapshot_1 = gcs_data[i + 1]
        event_gcs_out.append(do_gcs(input_data=(event_snapshot_0, event_snapshot_1),
                                    savepath=event_save_path,
                                    event_id=event_id))
    return event_gcs_out


def do_gcs(input_data: Tuple[List[str], List[str]], savepath: str, event_id: int) -> Tuple:

    """
        do_cgs: This function actually runs gcs for a given set of lasco, aia & secchi files
        savepath: Final path to save the ouput
        event_id: Event id to name the output
    """
    # check if we can run IDL gcs with pidly

    save_file = savepath + f'event_{event_id}.sav'  # check this on python i.e. how to save this data.

    snapshot_0, snapshot_1 = input_data

    # Read in fits files
    fnameA0, fnameB0, fnameL0 = snapshot_0
    fnameA1, fnameB1, fnameL1 = snapshot_1

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
    clouds = pyGCS.getGCS(0, 30., 50., 12., 0.3, 30, satpos, nleg=5, ncirc=20, ncross=40)

    # Missing How to save output into some kind of file
    # or what else to do here.

    return headers, ims, satpos, plotranges, clouds


def data_sort(view1_files: List[str], view2_files: List[str], view3_files: List[str], start_time: Time, end_time: Time,
              cadence: int) -> List:

    """
        data_gcs_sort:
        This function sort the file list into a chronologically ordered sets of lasco, cor, aia, euvi files
        and returns this list of ordered sets
        Inputs:
            lasco_files: list of lasco files in Map objects
            aia_files: list of aia files in Map objects
            secchia_files: list of euvia or cor2a files in Map objects
            secchib_files: list of euvib or cor2b files in Map objects
        Outputs:
            List of ordered sets
            example output:
                [ [aia_data1, euvia_data1, euvib_data1],
                  [aia_data2, euvia_data2, euvib_data2],
                  [aia_data3, euvia_data3, euvib_data3],
                  [...................................],
                  [aia_dataN, euvia_dataN, euvib_dataN] ]
    """

    event_duration = end_time - start_time

    time_tolerance = TimeDelta(GCS_DELTA_T)

    delta_t = TimeDelta(1 * apu.h) / cadence

    # Define list of time instants to match images times
    instants_number = int(event_duration.sec / delta_t.sec)
    time_instants = start_time + delta_t * np.arange(0, instants_number)

    # Match time with images
    view1_matched_fits = [match_time_with_image(time, view1_files, time_tolerance) for time in time_instants]
    view2_matched_fits = [match_time_with_image(time, view2_files, time_tolerance) for time in time_instants]
    view3_matched_fits = [match_time_with_image(time, view3_files, time_tolerance) for time in time_instants]

    # return zipped view1, view2, view3 files
    return list(zip(view1_matched_fits, view2_matched_fits, view3_matched_fits))


def match_time_with_image(time_instant: Time, images_files: List, time_tolerance: TimeDelta) -> str:

    """
        match_time_with_image:
            This function takes a Time object and a list of open fits files, and looks for the closest in time fits,
            taking into account max time tolerance
    """

    matched_files = []
    time_distances = []

    for file in images_files:
        file_time = get_time_from_file(file)
        time_instant2 = time_instant
        time_distance = abs(file_time - time_instant2)

        if time_distance.sec <= time_tolerance.sec:
            matched_files.append(file)
            time_distances.append(time_distance)

    # Return the closest one
    try:
        closest_one = matched_files[time_distances.index(min(time_distances))]
    except ValueError:
        closest_one = ""

    return closest_one


def get_time_from_file(filename: str) -> Time:

    hdr = fits.getheader(filename)

    if "lasco" in filename:
        return parse_time(hdr["DATE-OBS"] + " " + hdr['TIME-OBS'])
    else:
        return parse_time(hdr["DATE-OBS"])


def event_file_select(files: List[str], start_time: Time, end_time: Time) -> List[str]:

    """
        event_file_select: This function checks for events hours and return only
                            the needed list of files for the needed event
        files: list of files to select from.
        start_time:
        end_time:
    """

    fits_files = [fits.open(file) for file in files]

    # This may work on newer fits files but on some lasco files it will have to be DATE-OBS + TIME-OBS
    # haven´t tested it yet but should be working
    fits_date_obs = [parse_time(fits_file[0].header['DATE-OBS']) for fits_file in fits_files]

    start_index = 0
    end_index = 0

    # look for indexes
    for i in range(0, len(fits_date_obs)):
        if fits_date_obs[i] < start_time:
            start_index += 1
        if fits_date_obs[i] < end_time:
            end_index += 1

    try:
        selected_files = files[start_index:end_index+1]
    except IndexError:
        selected_files = None
        print("Error while figuring out start/end indexes")

    return selected_files


if __name__ == '__main__':

    # Ask for start && end times to look for and event
    event_start_time = parse_time(input("\nEnter start-time (YYYY-MM-DD hh:mm:ss format):\n"))
    event_end_time = parse_time(input("\nEnter start-time (YYYY-MM-DD hh:mm:ss format):\n"))
    event_number = int(input("\nEnter event number\n"))  # ask this

    run_event_gcs(start_time=event_start_time, end_time=event_end_time, cadence=2, do_low_corona=True, event_id=event_number)


