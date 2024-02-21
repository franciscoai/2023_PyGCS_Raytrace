import astropy.units as apu

DATA_PATH = '/media/gehme/gehme/data'
SECCHI_PATH = DATA_PATH + '/stereo/secchi/L1'

LASCO_PATH = DATA_PATH + '/soho/lasco/level_1/c2/'
AIA193_PATH = DATA_PATH + '/sdo/aia/l1/193/'
EUVI195_B_PATH = SECCHI_PATH + '/b/img/euvi/'
EUVI195_A_PATH = SECCHI_PATH + '/a/img/euvi/'
COR2A_PATH = SECCHI_PATH + '/a/img/cor2/'
COR2B_PATH = SECCHI_PATH + '/b/img/cor2/'

SAVE_PATH = '~/2021_cme_expansions_sources/GCS/'

GCS_DELTA_T = 10 * apu.min
