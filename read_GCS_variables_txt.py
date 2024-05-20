
def read_GCS_variables_txt(filename):
    dir='/data_local/GCS/2023_PyGCS_Raytrace/output/'
    with open(dir+filename, 'r') as file:
        # Read all lines
        lines = file.readlines()

    # Define variables to store values
    lon = None
    lat = None
    tilt = None
    height = None
    half_ang = None
    ratio = None
    scaling = None
    sat1min = None
    sat1max = None
    sat2min = None
    sat2max = None
    sat3min = None
    sat3max = None

    # Iterate through each line and extract values
    for line in lines:
        
        # Split the line by colon and strip whitespace
        key, value = map(str.strip, line.split(':'))
        
        # Assign value to the corresponding variable
        if key == 'Lon':
            lon = float(value)
        elif key == 'Lat':
            lat = float(value)
        elif key == 'Tilt':
            tilt = float(value)
        elif key == 'Height':
            height = float(value)
        elif key == 'HalfAng':
            half_ang = float(value)
        elif key == 'Ratio':
            ratio = float(value)
        elif key == 'Scaling':
            scaling = int(value)
        elif key == 'Sat1min':
            sat1min = int(value)
        elif key == 'Sat1max':
            sat1max = int(value)
        elif key == 'Sat2min':
            sat2min = int(value)
        elif key == 'Sat2max':
            sat2max = int(value)
        elif key == 'Sat3min':
            sat3min = int(value)
        elif key == 'Sat3max':
            sat3max = int(value)
    #order of gcs_parameters: lon, lat, tilt, height, ratio, half_ang
    return [lon, lat, tilt, height, ratio, half_ang]