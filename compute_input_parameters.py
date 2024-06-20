def compute_vs_average_on_fault(fault, layers):
    import numpy as np
    nsuby = fault['number_subfaults_dip']
    nsubx = fault['number_subfaults_strike']
    vs = layers['vs']
    thk = layers['thk']
    y_hypc = fault['hypo_down_dip']
    dy = fault['subfault_width']
    ncyp = int(y_hypc / dy) + 1
    cyp = (ncyp - 0.5) * dy
    y_hypc = cyp
    rdip = fault['dip'] * np.pi / 180
    depth_hypc = fault['hypo']['Z']

    nsubxy = nsubx*nsuby
    beta = np.zeros((nsubxy))
    for j in range(nsuby):
        yij = ((j+1) - 0.5) * dy - cyp
        hk = yij * np.sin(rdip) + depth_hypc
        thk[len(thk)-1] = hk + 0.5
        lyr = -1
        sumh = 0
        while hk > sumh:
            lyr = lyr + 1
            sumh = sumh + thk[lyr]
        sv = vs[lyr]
        for i in range(nsubx):
            k = ((i+1) - 1) * nsuby + j
            beta[k] = sv
    vs_average = 0.
    for k in range(nsubxy):
        vs_average = vs_average + beta[k]
    vs_average = vs_average/nsubxy
    return vs_average

def WC1994(rake, mag):
    #Wells and Coppersmith 1994
    if -45 <= rake <= 45 or rake >= 135 or rake <= -135:
        # strike slip
        area = 10.0 ** (-3.42 + 0.90 * mag)
        length = 10.0 ** (-2.57 + 0.62 * mag)
    elif rake > 0:
        # thrust/reverse
        area = 10.0 ** (-3.99 + 0.98 * mag)
        length = 10.0 ** (-2.42 + 0.58 * mag)
    else:
        # normal
        area = 10.0 ** (-2.87 + 0.82 * mag)
        length = 10.0 ** (-1.88 + 0.50 * mag)
    return area, length


def weibull_distribution(scale,shape):
    import numpy as np
    x = np.random.weibull(shape)
    return scale*x


def hypo_strike_Causse_2008():
    import numpy as np
    position_along_length = np.random.normal(loc=0.5, scale=0.23, size=None) #Il riferimento indicato nel 3Ptool è (Cotton 2008 BSSA)
    #Xnuc/L normal distribution mu=0.5, sigma=0.23
    if position_along_length < 1./11.:
        position_along_length = 1./11.
    if position_along_length > 10./11.:
        position_along_length = 10./11.
    return position_along_length


def hypo_dip_Causse_2008(rake):
    # Il riferimento indicato nel 3Ptool è (Cotton 2008 BSSA) Da controllare perchè nel 3Ptool usano una sola relazione
    if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
        # strike slip
        # Mai et al. (2005) Hypocenter Locations in Finite-Source Rupture Models strike-slip
        # Weibull distribution a=0.626 (scale parameter) b=3.921 Ynuc/L (shape parameter)
        position_along_width = weibull_distribution(0.626,3.921)
    elif rake > 0:
        # thrust/reverse
        # Mai et al. (2005) Hypocenter Locations in Finite-Source Rupture Models crustal dip-slip
        # Weibull distribution a=0.692   b=3.394 Ynuc/L
        position_along_width = weibull_distribution(0.692, 3.394)
    else:
        # normal
        # Mai et al. (2005) Hypocenter Locations in Finite-Source Rupture Models crustal dip-slip
        # Weibull distribution a=0.692   b=3.394 Ynuc/L
        position_along_width = weibull_distribution(0.692, 3.394)

    # set minimum and maximum possible values according to Causse et al. 2008
    if position_along_width < 0.5:
        position_along_width = 0.5
    if position_along_width > 0.7:
        position_along_width = 0.7
    return position_along_width


def define_slip_extended_source(slip_file, fault, computational_param):
    import os
    import numpy as np
    split_tup = os.path.splitext(slip_file)
    file_extension = split_tup[1]
    if file_extension == '.srcmod':
        #slip_txt = np.loadtxt(slip_file)
        nseg = 1
        segL = fault['length']
        segW = fault['width']
        segStk = fault['strike']
        segDip = fault['dip']
        Hypo_strike = fault['hypo_along_strike']
        Hypo_dip = fault['hypo_down_dip']
        nStk = fault['number_subfaults_strike']
        nDip = fault['number_subfaults_dip']
        ndiv = computational_param['ndiv']
        slip_orig = np.loadtxt(slip_file, comments='seg')

        dx_orig = segL / nStk
        dy_orig = segW / nDip

    return slip_orig


def find_origin_layer(layers, fault):
    thk_sum = 0
    if layers['thk'][len(layers['thk'])-1] == 0:
        layers['thk'][len(layers['thk']) - 1] = 1000
    for i in range(len(layers['thk'])):
        thk_sum = thk_sum + layers['thk'][i]
        if fault['hypo']['Z'] <= thk_sum:
            nlayer_hypo = i
            return nlayer_hypo
    return


def compute_slip_point_source(layers, fault):
    nlayer_hypo = find_origin_layer(layers, fault)
    mu = layers['rho'][nlayer_hypo]*1000 * layers['vs'][nlayer_hypo]*1000 * layers['vs'][nlayer_hypo]*1000
    area = fault['length']*1000 * fault['width']*1000
    slip = fault['Mo'] / (mu * area)
    return [slip]


def determine_origin_point(lon, lat, fault):
    import numpy as np
    from shapely.geometry import Point

    strike = fault['strike']
    dip = fault['dip']
    length = fault['length']
    width = fault['width']
    hypo_along_strike = fault['hypo_along_strike']
    hypo_down_dip = fault['hypo_down_dip']
    hypo_depth = fault['hypo']['Z']

    #
    #     .....     the dotted line is the hdist
    #     \      |
    #      \     |  this dashed vertical line is the height
    #       \    |
    #        \   |
    # rupture \  |
    #

    # if ztor is not None:
    #	depth = ztor + height / 2

    height = width * np.sin(np.radians(dip))
    hdist = width * np.cos(np.radians(dip))

    # Move hor. hypo_down_dip*hdist in direction -90
    hypo_top = point_at(lon, lat, strike - 90, hypo_down_dip * hdist)  # dist in km
    # Move hor. (1-hypo_down_dip)* hdist in direction +90
    hypo_bot = point_at(lon, lat, strike + 90, (1 - hypo_down_dip) * hdist)

    # compute corner points at the surface
    top_right = point_at(hypo_top[0], hypo_top[1], strike, hypo_along_strike * length)
    top_left = point_at(hypo_top[0], hypo_top[1], strike + 180, (1 - hypo_along_strike) * length)
    bot_right = point_at(hypo_bot[0], hypo_bot[1], strike, hypo_along_strike * length)
    bot_left = point_at(hypo_bot[0], hypo_bot[1], strike + 180, (1 - hypo_along_strike) * length)

    # compute corner points in 3D
    pbl = Point(bot_left[0], bot_left[1], hypo_depth + (1 - hypo_down_dip) * height)
    pbr = Point(bot_right[0], bot_right[1], hypo_depth + (1 - hypo_down_dip) * height)
    ptl = Point(top_left[0], top_left[1], hypo_depth - hypo_down_dip * height)
    ptr = Point(top_right[0], top_right[1], hypo_depth - hypo_down_dip * height)

    return pbl, pbr, ptl, ptr


def determine_utm_coord(lon, lat):
    import math
    import sys
    from pyproj import Proj
    if lon > 180 or lon < -180:
        sys.exit('Error with longitudes')
    utm_zone_num = math.ceil((lon + 180) / 6)
    print('UTM zone:', utm_zone_num)
    if lat > 0:
        hemisphere = 'north'
    else:
        hemisphere = 'south'
    myProj = Proj("+proj=utm +zone=" + str(utm_zone_num) + " +" + hemisphere + " +datum=WGS84 +units=m +no_defs ")
    UTMx, UTMy = myProj(lon, lat, inverse=False)
    # x Easting
    # y Northing
    return UTMx, UTMy


def point_at(lon, lat, azimuth, distance):
    """
	Perform a forward geodetic transformation: find a point lying at a given
	distance from a given one on a great circle arc defined by azimuth.

	:param float lon, lat:
	    Coordinates of a reference point, in decimal degrees.
	:param azimuth:
	    An azimuth of a great circle arc of interest measured in a reference
	    point in decimal degrees.
	:param distance:
	    Distance to target point in km.
	:returns:
	    Tuple of two float numbers: longitude and latitude of a target point
	    in decimal degrees respectively.

	Implements the same approach as :func:`npoints_towards`.
	"""
    import numpy as np

    EARTH_RADIUS = 6371.0
    # this is a simplified version of points_towards().
    # code duplication is justified by performance reasons.
    lon, lat = np.radians(lon), np.radians(lat)
    tc = np.radians(360 - azimuth)
    sin_dists = np.sin(distance / EARTH_RADIUS)
    cos_dists = np.cos(distance / EARTH_RADIUS)
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)

    sin_lats = sin_lat * cos_dists + cos_lat * sin_dists * np.cos(tc)
    lats = np.degrees(np.arcsin(sin_lats))

    dlon = np.arctan2(np.sin(tc) * sin_dists * cos_lat, cos_dists - sin_lat * sin_lats)
    lons = np.mod(lon - dlon + np.pi, 2 * np.pi) - np.pi
    lons = np.degrees(lons)

    return lons, lats


def create_computational_param():
    dt = None
    npts = None
    fmin = None
    fmax = None
    kappa = None
    spectral_degree = None
    fval_quality_factor = None
    damping = None
    monfile = None
    mpifile = None
    optiout = None
    timestep = None
    tmonitor = None
    stoptime = None
    mlst = None
    setuponl = None
    ndiv = None
    computational_param = {
        "dt": dt,
        "npts": npts,
        "fmax": fmax,
        "fmin": fmin,
        "kappa": kappa,
        "spectral_degree" : spectral_degree,
        "fval_quality_factor" : fval_quality_factor,
        "damping" : damping,
        "monfile" : monfile,
        "mpifile" : mpifile,
        "optiout" : optiout,
        "timestep" : timestep,
        "tmonitor" : tmonitor,
        "stoptime" : stoptime,
        "mlst" : mlst,
        "setuponl" : setuponl,
        "ndiv" : ndiv

    }
    return computational_param


def create_point_utm(point_x, point_y, point_z):
    point = {
        "X": point_x,
        "Y": point_y,
        "Z": point_z
    }
    return point


def create_point(point_lon, point_lat, point_z):
    point = {
        "lon": point_lon,
        "lat": point_lat,
        "Z": point_z
    }
    return point


def create_fault():
    length = 0.100  # in km
    width = 0.100  # in km
    strike = None
    dip = None
    rake = None
    rake_distribution = None
    rise_time = None
    rupture_velocity = None
    number_subfaults_strike = 1
    number_subfaults_dip = 1
    subfault_length = length
    subfault_width = width
    hypo = None
    hypo_utm = None
    origin = None
    slip = None
    Mo = None
    hypo_along_strike = 0.5  # The reference is in the top-left corner of the rupture plane when viewed at an angle of 90° from strike
    hypo_down_dip = 0.5
    STF = None
    Tp_Tr = None
    fault_type = None
    Ztor = 0
    depth_max_fault = 999
    fault = {
        "length": length,
        "width": width,
        "strike": strike,
        "dip": dip,
        "rake": rake,
        "rise_time": rise_time,
        "rupture_velocity": rupture_velocity,
        "number_subfaults_strike": number_subfaults_strike,
        "number_subfaults_dip": number_subfaults_dip,
        "hypo": hypo,
        "hypo_utm": hypo_utm,
        "origin": origin,
        "slip": slip,
        "Mo": Mo,
        "hypo_along_strike": hypo_along_strike,
        "hypo_down_dip": hypo_down_dip,
        "IDx": STF,
        "fault_type": fault_type,
        "Tp_Tr": Tp_Tr,
        "rake_distribution": rake_distribution,
        "subfault_length": subfault_length,
        "subfault_width": subfault_width,
        "Ztor": Ztor,
        "depth_max_fault" : depth_max_fault
    }
    return fault


def create_material_properties(vp, vs, qp, qs, rho, thk, fqp, fqs):
    layers = {
        "vp": vp,
        "vs": vs,
        "qp": qp,
        "qs": qs,
        "fqp": fqp,
        "fqs": fqs,
        "rho": rho,
        "thk": thk
    }
    return layers


def create_receivers(receivers_X, receivers_Y, receivers_Z, receivers_ID):
    sites = {
        "X": receivers_X,
        "Y": receivers_Y,
        "Z": receivers_Z,
        "ID": receivers_ID
    }
    return sites


def create_plot_param():
    fmin_filter = None
    fmax_filter = None
    time_max_plot = None
    plot_param = {
        "fmin_filter": fmin_filter,
        "fmax_filter": fmax_filter,
        "time_max_plot": time_max_plot
    }
    return plot_param


def read_input_data(fileini, code, calculation_mode):
    import sys
    import numpy as np
    import math

    input = {}
    with open(fileini) as fp:
        line = fp.readline()
        while line:
            if line.strip().find('=') >= 0:
                key, value = line.strip().split('=', 1)
                input[key.strip()] = value.strip()
            line = fp.readline()

    vp = [x.strip() for x in input['vp'].strip('[]').split(',')]
    vp = np.array(vp, dtype=float)
    vs = [x.strip() for x in input['vs'].strip('[]').split(',')]
    vs = np.array(vs, dtype=float)

    thk = [x.strip() for x in input['thk'].strip('[]').split(',')]
    thk = np.array(thk, dtype=float)
    try:
        rho = [x.strip() for x in input['rho'].strip('[]').split(',')]
        rho = np.array(rho, dtype=float)
    except KeyError:
        rho = np.zeros((len(vp)))
        for i in range(len(vp)):
            #“Nafe-Drake curve” (Ludwig et al., 1970).
            #Ludwig, W. J., J. E. Nafe, and C. L. Drake (1970).
            # Seismic refraction, in The Sea, A. E. Maxwell, (Editor) Vol. 4, Wiley-Interscience, New York, 53–84.
            if vp[i] >= 1.5 and vp[i] <= 8.5:
                rho[i] = 1.6612*vp[i]-0.4721*vp[i]**2+0.0671*vp[i]**3-0.0043*vp[i]**4+0.000106*vp[i]**5
            else:
                sys.exit('Error: TO DO - need to find a regression for rho')
        print('Velocity model: rho computed through Nafe-Drake relationships')

    try:
        qs = [x.strip() for x in input['qs'].strip('[]').split(',')]
        qs = np.array(qs, dtype=float)
    except KeyError:
        qs = np.zeros((len(vs)))
        for i in range(len(vs)):
            if vs[i] < 0.5:
                qs[i] = 10
                # (based on shallow borehole studies by Tullos and Reid (1969), Hamilton (1972), Gibbs et al.(1994), Liu et al.(1994),
                # and Kudo and Shima (1970)).San Francisco Bay Area specific studies include Gibbs et al.(1994) and Liu et al.(1994).
            elif vs[i] >= 0.5 and vs[i] < 1.5:
                qs[i] = 20*vs[i] # vs in km/s (Olsen et al., 2003)
            else:
                qs[i] = 100*vs[i] #Olsen et al., 2003
        print('Velocity model: vs computed from qs')

    try:
        qp = [x.strip() for x in input['qp'].strip('[]').split(',')]
        qp = np.array(qp, dtype=float)
    except KeyError:
        qp = np.zeros((len(qs)))
        for i in range(len(vs)):
            qp[i] = 9/4*qs[i] #Lay, T., and Wallace, T. (1995). Modern Global Seimology, Vol. 58.
        print('Velocity model: qp computed from qs')


    fqp = None
    fqs = None
    if 'hisada' in code:
        fqp = [x.strip() for x in input['fqp'].strip('[]').split(',')]
        fqp = np.array(fqp, dtype=float)
        fqs = [x.strip() for x in input['fqs'].strip('[]').split(',')]
        fqs = np.array(fqs, dtype=float)

    layers = create_material_properties(vp, vs, qp, qs, rho, thk, fqp, fqs)

    fault = create_fault()
    lon_hypo = float(input['lon_hypo'])
    lat_hypo = float(input['lat_hypo'])
    depth_hypo = float(input['depth_hypo'])

    fault['fault_type'] = input['fault_type']

    fault['strike'] = float(input['strike'])
    fault['dip'] = float(input['dip'])
    fault['rake'] = float(input['rake'])
    if fault['rake'] > 180:
        fault['rake'] = fault['rake'] - 360

    if fault['fault_type'] == 'extended':

        try:
            fault['length'] = float(input['length'])
        except KeyError:
            areaWC94, lengthWC94 = WC1994(fault['rake'], fault['Mw'])
            try:
                fault['width'] = float(input['width'])
                fault['length'] = areaWC94 / fault['width']
                print('Fault: length computed trough width and Wells and Coppersmith relationships: ', fault['length'])
            except KeyError:
                fault['length'] = lengthWC94
                print('Fault: length computed trough Wells and Coppersmith relationships: ', fault['length'])

        try:
            fault['width'] = float(input['width'])
        except KeyError:
            areaWC94, lengthWC94 = WC1994(fault['rake'], fault['Mw'])
            fault['width'] = areaWC94 / fault['length']
            print('Fault: width computed trough Wells and Coppersmith relationships: ', fault['width'])

        try:
            fault['number_subfaults_strike'] = int(input['number_subfaults_strike'])
            fault['subfault_length'] = fault['length']/input['number_subfaults_strike']
        except KeyError:
            try:
                fault['subfault_length'] = float(input['subfault_length'])
                fault['number_subfaults_strike'] = math.ceil(fault['length']/fault['subfault_length'])
            except KeyError:
                sys.exit('Error: number_subfaults_strike/subfault_length not defined')

        try:
            fault['number_subfaults_dip'] = int(input['number_subfaults_dip'])
            fault['subfault_width'] = fault['width'] / input['number_subfaults_dip']
        except KeyError:
            try:
                fault['subfault_width'] = float(input['subfault_width'])
                fault['number_subfaults_dip'] = math.ceil(fault['width'] / fault['subfault_width'])
            except KeyError:
                sys.exit('Error: number_subfaults_dip/subfault_width not defined')

        if 'ucsb' in code or fault['slip_mode'] == 'Archuleta':
            try:
                fault['subfault_length'] = float(input['subfault_length'])
            except KeyError:
                sys.exit('Error: subfault_length not defined')
            try:
                fault['subfault_width'] = float(input['subfault_width'])
            except KeyError:
                sys.exit('Error: subfault_width not defined')

        try:
            fault['hypo_along_strike'] = float(input['hypo_along_strike'])
        except KeyError:
            hypo_strike = hypo_strike_Causse_2008()
            fault['hypo_along_strike'] = hypo_strike

        try:
            fault['hypo_down_dip'] = float(input['hypo_down_dip'])
        except KeyError:
            fault['hypo_down_dip'] = hypo_dip_Causse_2008(fault['rake'])

    hypo = create_point(lon_hypo, lat_hypo, depth_hypo)  # in degrees and depth in km
    fault['hypo'] = hypo

    fault['Ztor'] = fault['hypo']['Z'] - fault['hypo_down_dip'] * fault['width'] * np.sin(np.radians(fault['dip']))
    try:
        min_allowed_fault_depth = float(input['min_fault_depth'])
    except KeyError:
        min_allowed_fault_depth = 0.
    try:
        max_allowed_fault_depth = float(input['max_fault_depth'])
    except KeyError:
        max_allowed_fault_depth = 999.

    if fault['Ztor'] < min_allowed_fault_depth:
        fault['Ztor'] = min_allowed_fault_depth
        fault['hypo']['Z'] = fault['Ztor'] + fault['hypo_down_dip'] * fault['width'] * np.sin(np.radians(fault['dip']))

    fault['depth_max_fault'] = fault['Ztor'] + fault['width'] * np.sin(np.radians(fault['dip']))
    if fault['depth_max_fault'] > max_allowed_fault_depth:
        sys.exit('Error: deepest fault point exceeds the maximum fault depth')

    hypo_x, hypo_y = determine_utm_coord(lon_hypo, lat_hypo)
    hypo_utm = create_point_utm(hypo_y, hypo_x, depth_hypo * 1000)  # in m
    fault['hypo_utm'] = hypo_utm

    fault['IDx'] = input['STF']
    try:
        fault['rupture_velocity'] = float(input['rupture_velocity'])
    except KeyError:
        vs_average = compute_vs_average_on_fault(fault, layers)
        try:
            fault['percentage_rupture_velocity'] = float(input['percentage_rupture_velocity'])
            fault['rupture_velocity'] = fault['percentage_rupture_velocity'] * vs_average
        except KeyError:
            # span the range of rupture velocities expected for most earthquakes (e.g. Heaton 1990)
            # Heaton T. Evidence for and implications of self-healing pulses of slip in earthquake rupture,
            # Phys. Earth planet. Inter., 1990, vol. 64 (pg. 1-20)
            fault['rupture_velocity'] = np.random.uniform(0.65, 0.85) * vs_average
            # upture velocity values reported in source studies mainly range between 0.65Vs and 0.85Vs [e.g., Heaton, 1990].
        print('Fault: computed rupture velocity = ', fault['rupture_velocity'])

    try:
        fault['rise_time'] = float(input['rise_time'])
    except KeyError:
        if fault['fault_type'] == 'extended':
            TD = fault['length'] / fault['rupture_velocity']
            fault['rise_time'] = 0.1 * TD  # Haskell 1964; Aki 1967; Heaton 1990
        else:
            sys.exit('The rise time must be specified for a point source')

    try:
        fault['Mo'] = float(input['Mo'])
    except KeyError:
        try:
            mag = float(input['Mw'])
            fault['Mo'] = 10. ** (1.5 * mag + 9.05)  # in Nm Hanks & Kanamori (1979)
            print('Fault: Mo computed from Mw: ', fault['Mo'])
        except KeyError:
            sys.exit('Error: Mo/Mw not found')

    try:
        fault['mag'] = float(input['Mw'])
    except KeyError:
        try:
            Mo = float(input['Mo'])
            fault['Mw'] = 2/3 * math.log10(Mo) - 9.05
            print('Fault: Mw computed from Mo: ', fault['Mw'])
        except KeyError:
            sys.exit('Error: Mo/Mw not found')

    if fault['fault_type'] == 'point':
        slip = compute_slip_point_source(layers, fault)
        fault['slip'] = slip
        fault['rake_distribution'] = [fault['rake']]

    if 'hisada' in code:
        pbl, pbr, ptl, ptr = determine_origin_point(lon_hypo, lat_hypo, fault)
        origin_point = pbl
        origin_utmx, origin_utmy = determine_utm_coord(origin_point.x, origin_point.y)
        origin = create_point_utm(origin_utmy, origin_utmx, origin_point.z*1000)
        fault['origin'] = origin

    if 'ucsb' in code:
        fault['Tp_Tr'] = float(input['Tp_Tr'])
        origin_point = determine_origin_point(lon_hypo, lat_hypo, fault)


    computational_param = create_computational_param()
    if 'speed' in code:
        computational_param['spectral_degree'] = int(input['spectral_degree'])
        computational_param['fval_quality_factor'] = float(input['fval_quality_factor'])
        computational_param['damping'] = int(input['damping'])
        computational_param['monfile'] = input['monfile']
        computational_param['mpifile'] = input['mpifile']
        optiout = [x.strip() for x in input['optiout'].strip('[]').split(',')]
        computational_param['optiout'] = np.array(optiout, dtype=int)
        computational_param['timestep'] = float(input['timestep'])
        computational_param['tmonitor'] = int(input['tmonitor'])
        computational_param['stoptime'] = float(input['stoptime'])
        mlst = [x.strip() for x in input['mlst'].strip('[]').split(',')]
        computational_param['mlst'] = np.array(mlst)
        computational_param['setuponl'] = input['setuponl']
        try:
            computational_param['ndiv'] = input['ndiv']
            # Should be 4 ^ n, where n is a real integer = 0, 1, 2, 3....
            # Each rectangular subfault in Each Fault segment will be divided
            # into ndiv rectangles.This is to increase No.of Point Sources
            # representing Fault Plane
            # Suggested values in ThreePtool 1, 4, 16, 64, 256, 1024
        except KeyError:
            computational_param['ndiv'] = 1
    if 'hisada' in code:
        computational_param['dt'] = float(input['dt'])
        computational_param['npts'] = int(input['npts'])
        computational_param['fmax'] = float(input['fmax'])

    if 'ucsb' in code:
        computational_param['fmin'] = float(input['fmin'])
        computational_param['fmax'] = int(input['fmax'])
        computational_param['kappa'] = float(input['kappa'])
        computational_param['dt'] = float(input['dt'])
        computational_param['npts'] = int(input['npts'])

    if fault['fault_type'] == 'extended' and 'ucsb' not in code:
        try:
            slip_file = input['slip_file']
        except KeyError:
            try:
                slip_mode = input['slip_mode']
                if slip_mode == 'Archuleta':
                    fault['slip'] = 'Archuleta'
            except KeyError:
                sys.exit('Error: the slip file must be defined for an extended seismic source')
        #slip = define_slip_extended_source(slip_file, fault, computational_param)

    receivers_lat = [x.strip() for x in input['receivers_lat'].strip('[]').split(',')]
    receivers_lat = np.array(receivers_lat, dtype=float)
    receivers_lon = [x.strip() for x in input['receivers_lon'].strip('[]').split(',')]
    receivers_lon = np.array(receivers_lon, dtype=float)
    try:
        receivers_depth = [x.strip() for x in input['receivers_depth'].strip('[]').split(',')]
        receivers_depth = np.array(receivers_depth, dtype=float)
    except KeyError:
        receivers_depth = np.zeros((len(receivers_lat)))
    try:
        receivers_ID = [x.strip() for x in input['receivers_ID'].strip('[]').split(',')]
        receivers_ID = np.array(receivers_ID, dtype=str)
    except KeyError:
        receivers_ID = []
        for l in range(len(receivers_depth)):
            receivers_ID.append('site'+str(l+1).zfill(3))
        receivers_ID = np.asarray(receivers_ID)

    path_code_speed = None
    path_code_hisada = None
    path_code_ucsb = None
    folder_speed = None
    folder_hisada = None
    folder_ucsb = None
    if 'hisada' in code:
        folder_hisada = input['folder_hisada']
        path_code_hisada = input['path_code_hisada']

    if 'speed' in code:
        folder_speed = input['folder_speed']
        path_code_speed = input['path_code_speed']

    if 'ucsb' in code:
        folder_ucsb = input['folder_ucsb']
        path_code_ucsb = input['path_code_ucsb']

    receivers_X = []
    receivers_Y = []
    for i in range(len(receivers_depth)):
        receivers_utmx, receivers_utmy = determine_utm_coord(receivers_lon[i], receivers_lat[i])
        receivers_X.append(receivers_utmx)
        receivers_Y.append(receivers_utmy)
    sites = create_receivers(receivers_Y, receivers_X, receivers_depth*1000, receivers_ID)

    plot_param = create_plot_param()
    folder_plot = None
    if calculation_mode == '--plot':
        try:
            plot_param['fmin_filter'] = float(input['fmin_filter'])
        except KeyError:
            plot_param['fmin_filter'] = None
        try:
            plot_param['fmax_filter'] = float(input['fmax_filter'])
        except KeyError:
            plot_param['fmax_filter'] = None
        plot_param['time_max_plot'] = float(input['time_max_plot'])
        folder_plot = input['folder_plot']

    return layers, fault, computational_param, sites, path_code_speed, path_code_hisada, path_code_ucsb, \
        folder_speed, folder_hisada, folder_ucsb, folder_plot, plot_param
