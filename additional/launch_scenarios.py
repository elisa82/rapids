def write_input_earthquake(lon, lat, depth, Mo, strike, dip, rake):
    input_rt = {
        "lon": lon,
        "lat": lat,
        "depth": depth,
        "Mo": Mo,
        "strike": strike,
        "dip": dip,
        "rake": rake
    }
    return input_rt

def read_moment_tensor_file(file_input):
    import numpy as np
    input = {}
    with open(file_input) as fp:
        line = fp.readline()
        while line:
            if line.strip().find('=') >= 0:
                key, value = line.strip().split('=', 1)
                input[key.strip()] = value.strip()
            line = fp.readline()
    strike_double_planes = [x.strip() for x in input['Strike'].strip('').split(';')]
    strike_arr = np.array(strike_double_planes, dtype=float)
    rake_double_planes = [x.strip() for x in input['Rake'].strip('').split(';')]
    rake_arr = np.array(rake_double_planes, dtype=float)
    dip_double_planes = [x.strip() for x in input['Dip'].strip('').split(';')]
    dip_arr = np.array(dip_double_planes, dtype=float)
    Mo = float(input['Mo'])
    line_depth = input['isoflag']
    depth = line_depth.split("=")[1]
    return strike_arr, dip_arr, rake_arr, Mo, depth

#manca along strike e along dip
def write_fileini_rapids(fileini, input_rt, folder):
    f = open(fileini, "w")
    f.write('{}{}\n'.format('output_folder = ', folder))
    f.write('{}\n'.format(''))
    f.write('{}\n'.format('#velocity model'))
    f.write('{}\n'.format('vel_model = NAC_1D_Friuli'))
    f.write('{}\n'.format(''))
    f.write('{}\n'.format('topography = 1'))
    f.write('{}\n'.format('fault_geolocation = from_hypo'))
    f.write('{}\n'.format('#fault'))

    Mw = 2 / 3 * math.log10(input_rt['Mo']) - 9.05  # in Nm Hanks & Kanamori (1979)

    print(Mw)
    threshold_ext = 4.5
    if Mw >= threshold_ext:
        fault_type = 'extended'
        if Mw > 5.5:
            fault_slip_mode = 'Archuleta'
        else:
            fault_slip_mode = 'constant'
    else:
        fault_type = 'point'

    f.write('{} {}\n'.format('fault_type =', fault_type))
    if Mw >= threshold_ext:
        f.write('{} {}\n'.format('slip_mode =', fault_slip_mode))
    f.write('{} {}\n'.format('strike =', input_rt['strike']))
    f.write('{} {}\n'.format('dip =', input_rt['dip']))
    f.write('{} {}\n'.format('rake =', input_rt['rake']))
    f.write('{} {}\n'.format('lon_hypo =', input_rt['lon']))
    f.write('{} {}\n'.format('lat_hypo =', input_rt['lat']))
    f.write('{} {}\n'.format('depth_hypo =', input_rt['depth']))
    f.write('{} {}\n'.format('Mo =', input_rt['Mo']))
    f.write('{}\n'.format(''))
    f.write('{}\n'.format('# Computational params SPEED'))
    f.write('{}\n'.format('spectral_degree = 4'))
    f.write('{}\n'.format('fval_quality_factor = 3'))
    f.write('{}\n'.format('damping = 2'))
    f.write('{}\n'.format('dt_speed = 0.001'))
    f.write('{}\n'.format('tmonitor = 10'))
    f.write('{}\n'.format('stoptime = 50'))
    f.write('{}\n'.format('mlst = [-100, 0]'))
    f.write('{}\n'.format('setuponl = F'))
    f.write('{}\n'.format('ndiv = 16'))
    f.write('{}\n'.format('npplambda = 4'))
    f.write('{}\n'.format('fmax_speed = 1.5'))
    f.write('{}\n'.format('distance_absorbing_boundaries = 10'))
    f.write('{}\n'.format(''))
    f.write('{}\n'.format('# Computational params UCSB'))
    f.write('{}\n'.format('kappa = 0.03'))
    f.write('{}\n'.format('fmin_ucsb = 1.0'))
    f.write('{}\n'.format('Tp_Tr = 0.2'))
    f.write('{}\n'.format('npts_ucsb = 2048'))
    f.write('{}\n'.format('dt_ucsb = 0.01'))
    f.write('{}\n'.format('fmax_ucsb = 10'))
    f.write('{}\n'.format('realizations = 1'))
    f.write('{}\n'.format('seed = [247, -559, -123]'))
    f.write('{}\n'.format(''))
    f.write('{}\n'.format('# receivers'))
    f.write('{}\n'.format('site_grid_step_km = 2'))
    f.write('{}\n'.format('site_maximum_dist_km_N = 20'))
    f.write('{}\n'.format('site_maximum_dist_km_S = 20'))
    f.write('{}\n'.format('site_maximum_dist_km_E = 20'))
    f.write('{}\n'.format('site_maximum_dist_km_W = 20'))
    f.write('{}\n'.format(''))
    f.write('{}\n'.format('#plot'))
    f.write('{}\n'.format('fmax_filter = 1.5'))
    f.write('{}\n'.format('time_max_plot = 50'))
    f.write('{}\n'.format(''))
    f.write('{}\n'.format('#cineca'))
    f.write('{}\n'.format('nnodes = 1'))
    f.write('{}\n'.format('nproc = 1'))
    f.write('{}\n'.format('memory = 5'))
    f.write('{}\n'.format('partition = g100_usr_dbg'))
    f.write('{}\n'.format('duration = 1:00:00'))
    f.write('{}\n'.format('ntask_per_node = 1'))
    f.write('{}\n'.format('account = IscrC_URASS-2'))
    f.write('{}\n'.format('job_name = Friuli'))
#cineca
#nnodes = 20
#nproc = nnodes * 48
#memory = 300
#partition = 'g100_usr_prod'
#duration = '2:00:00'
    f.close()
    return

import sys
import math

folder = sys.argv[1]

fileini = 'rts.ini'
mag = 5.5
lat = 46.00
lon = 13.00
depth = 10.00

moment_tensor_file = 'mt_inv.out'
strike_arr, dip_arr, rake_arr, Mo, depth = read_moment_tensor_file(moment_tensor_file)
strike = strike_arr[0]
dip = dip_arr[0]
rake = rake_arr[0]

input_rt = write_input_earthquake(lon, lat, depth, Mo, strike, dip, rake)

write_fileini_rapids(fileini, input_rt, folder)
