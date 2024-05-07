def compute_maximum_distance(mag):
    import rapids.bragato_slejko_2005 as bragato_slejko_2005
    import numpy as np

    #to define threshold and dist according to fragility functions of target structures
    threshold = 0.01

    dist = 0
    while dist < 300:
        mean, stdv = bragato_slejko_2005.BragatoSlejko2005().ground_motion('PGA', mag, dist)
        if np.exp(mean-stdv) < threshold:
            break
        dist += 1
    return dist


def grid_sites(fault, receiver_grid_step_km, receiver_maximum_dist):
    import math
    import numpy as np
    from rapids.conversions import determine_utm_coord, utm_to_lon_lat

    nstep_grid = math.floor(receiver_maximum_dist / receiver_grid_step_km) + 1
    receivers_lat = []
    receivers_lon = []
    receiver_grid_step = receiver_grid_step_km * 1000
    hypo_x, hypo_y, zone, letter = determine_utm_coord(fault['hypo']['lon'], fault['hypo']['lat'])
    x0 = hypo_x - nstep_grid * receiver_grid_step
    y0 = hypo_y - nstep_grid * receiver_grid_step
    for i in range(2 * nstep_grid + 1):
        x_grid = x0 + i * receiver_grid_step
        for j in range(2 * nstep_grid + 1):
            y_grid = y0 + j * receiver_grid_step
            lon_grid, lat_grid = utm_to_lon_lat(x_grid, y_grid, zone)
            receivers_lon.append(lon_grid)
            receivers_lat.append(lat_grid)
    receivers_lon = np.asarray(receivers_lon)
    receivers_lat = np.asarray(receivers_lat)
    receivers_depth = np.zeros((len(receivers_lat)))
    print(receivers_lon, receivers_lat)
    receivers_ID = []
    for l in range(len(receivers_depth)):
        receivers_ID.append('site' + str(l + 1).zfill(3))
    receivers_ID = np.asarray(receivers_ID)
    return receivers_lon, receivers_lat, receivers_depth, receivers_ID


def grid_sites_4lengths(fault, receiver_grid_step_km, receiver_maximum_dist_N, receiver_maximum_dist_S,
                       receiver_maximum_dist_E, receiver_maximum_dist_W):
    import math
    import numpy as np
    from rapids.conversions import determine_utm_coord, utm_to_lon_lat

    nstep_grid_N = math.floor(receiver_maximum_dist_N / receiver_grid_step_km) + 1
    nstep_grid_S = math.floor(receiver_maximum_dist_S / receiver_grid_step_km) + 1
    nstep_grid_E = math.floor(receiver_maximum_dist_E / receiver_grid_step_km) + 1
    nstep_grid_W = math.floor(receiver_maximum_dist_W / receiver_grid_step_km) + 1
    receivers_lat = []
    receivers_lon = []
    receiver_grid_step = receiver_grid_step_km * 1000
    hypo_x, hypo_y, zone, letter = determine_utm_coord(fault['hypo']['lon'], fault['hypo']['lat'])
    x0 = hypo_x - nstep_grid_W * receiver_grid_step
    y0 = hypo_y - nstep_grid_S * receiver_grid_step
    for i in range(nstep_grid_W + nstep_grid_E + 1):
        x_grid = x0 + i * receiver_grid_step
        for j in range(nstep_grid_S + nstep_grid_N + 1):
            y_grid = y0 + j * receiver_grid_step
            lon_grid, lat_grid = utm_to_lon_lat(x_grid, y_grid, zone)
            receivers_lon.append(lon_grid)
            receivers_lat.append(lat_grid)
    receivers_lon = np.asarray(receivers_lon)
    receivers_lat = np.asarray(receivers_lat)
    receivers_depth = np.zeros((len(receivers_lat)))
    receivers_ID = []
    for l in range(len(receivers_depth)):
        receivers_ID.append('site' + str(l + 1).zfill(3))
    receivers_ID = np.asarray(receivers_ID)
    return receivers_lon, receivers_lat, receivers_depth, receivers_ID


def create_computational_param():
    import numpy as np
    dt_ucsb = None
    dt_speed = None
    dt_hisada = None
    npts_ucsb = None
    npts_hisada = None
    fmax_speed = None
    fmax_ucsb = None
    fmin_ucsb = None
    fmax_hisada = None
    kappa = None
    qzero = None
    alpha = None
    spectral_degree = None
    fval_quality_factor = None
    damping = None
    optiout = np.zeros((1, 6))[0]
    tmonitor = None
    stoptime = None
    mlst = None
    setuponl = None
    ndiv = None
    seed = None
    output_type = None
    abso_bound_dist_km = None
    freq_join = None
    realizations = 1
    computational_param = {
        "dt_ucsb": dt_ucsb,
        "dt_speed": dt_speed,
        "dt_hisada": dt_speed,
        "npts_ucsb": npts_ucsb,
        "npts_hisada": npts_hisada,
        "fmax_ucsb": fmax_ucsb,
        "fmax_speed": fmax_speed,
        "fmax_hisada": fmax_hisada,
        "fmin_ucsb": fmin_ucsb,
        "kappa": kappa,
        "qzero": qzero,
        "alpha": alpha,
        "spectral_degree": spectral_degree,
        "fval_quality_factor": fval_quality_factor,
        "damping": damping,
        "optiout": optiout,
        "output_type": output_type,
        "tmonitor": tmonitor,
        "stoptime": stoptime,
        "mlst": mlst,
        "setuponl": setuponl,
        "ndiv": ndiv,
        "seed": seed,
        "abso_bound_dist_km": abso_bound_dist_km,
        "freq_join": freq_join,
        "realizations": realizations
    }
    return computational_param


def create_fault():
    length = None  # in km
    width = None  # in km
    strike = None
    dip = None
    rake = None
    rise_time = None
    rupture_velocity = None
    number_subfaults_strike = None
    number_subfaults_dip = None
    subfault_length = None
    subfault_width = None
    hypo = None
    hypo_utm = None
    origin = None
    slip = None
    Mw = None
    Mo = None
    hypo_along_strike = None  # The reference is in the top-left corner of the rupture plane when viewed at an angle of 90° from strike
    hypo_down_dip = None
    STF = None
    Tp_Tr = 0.2
    fault_type = None
    min_fault_depth = 0
    max_fault_depth = 999
    percentage_rupture_velocity = None
    # 	#Allmann, B. P. and P. M. Shearer (2009). Global variations of stress drop for moderate to large
    # 	#earthquakes, J. Geophys. Res. 114, B01310, doi:10.1029/2008JB005821, 22 pp.
    stress_drop = None
    slip_mode = None
    slip_file = None
    fc = None
    vertex = None
    vertex_utm = None
    fault_geolocation = None
    rupture_duration = None
    nsubcells_length_xta = None
    nsubcells_width_xta = None
    file_xta = ''
    Ztor = None
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
        "Mw": Mw,
        "hypo_along_strike": hypo_along_strike,
        "hypo_down_dip": hypo_down_dip,
        "IDx": STF,
        "fault_type": fault_type,
        "Tp_Tr": Tp_Tr,
        "subfault_length": subfault_length,
        "subfault_width": subfault_width,
        "min_fault_depth": min_fault_depth,
        "max_fault_depth": max_fault_depth,
        "stress_drop": stress_drop,
        "percentage_rupture_velocity": percentage_rupture_velocity,
        "slip_mode": slip_mode,
        "slip_file": slip_file,
        "fc": fc,
        "vertex": vertex,
        "vertex_utm": vertex_utm,
        "fault_geolocation": fault_geolocation,
        "rupture_duration": rupture_duration,
        "file_xta": file_xta,
        "nsubcells_length_xta": nsubcells_length_xta,
        "nsubcells_width_xta": nsubcells_width_xta,
        "Ztor": Ztor
    }
    return fault


def create_cineca_slurm():
    nnodes = None
    memory = None  
    partition = None
    duration = None
    ntask_per_node = None
    account = None
    job_name = None

    cineca = {
        "nnodes": nnodes,
        "memory": memory,
        "partition": partition,
        "duration": duration,
        "ntask_per_node": ntask_per_node,
        "account": account,
        "job_name": job_name
    }
    return cineca


def create_material_properties():
    vp = None
    vs = None
    qp = None
    qs = None
    rho = None
    thk = None
    fqp = None
    fqs = None
    vel_model = None
    depth_top_layer = 0
    layers = {
        "vp": vp,
        "vs": vs,
        "qp": qp,
        "qs": qs,
        "fqp": fqp,
        "fqs": fqs,
        "rho": rho,
        "thk": thk,
        "vel_model": vel_model,
        "depth_top_layer": depth_top_layer
    }
    return layers


def create_receivers(receivers_X, receivers_Y, receivers_Z, receivers_ID, receivers_lon, receivers_lat):
    sites = {
        "X": receivers_X,
        "Y": receivers_Y,
        "Z": receivers_Z,
        "ID": receivers_ID,
        "lon": receivers_lon,
        "lat": receivers_lat
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
    from shapely.geometry import Point
    from rapids.conversions import determine_utm_coord, create_point, create_vertex
    import math

    layers = create_material_properties()
    fault = create_fault()

    input = {}
    with open(fileini) as fp:
        line = fp.readline()
        while line:
            if line.strip().find('=') >= 0:
                key, value = line.strip().split('=', 1)
                input[key.strip()] = value.strip()
            line = fp.readline()

    try:
        layers['vel_model'] = input['vel_model']
    except KeyError:
        vp = [x.strip() for x in input['vp'].strip('[]').split(',')]
        layers['vp'] = np.array(vp, dtype=float)
        vs = [x.strip() for x in input['vs'].strip('[]').split(',')]
        layers['vs'] = np.array(vs, dtype=float)

        thk = [x.strip() for x in input['thk'].strip('[]').split(',')]
        layers['thk'] = np.array(thk, dtype=float)
        try:
            rho = [x.strip() for x in input['rho'].strip('[]').split(',')]
            layers['rho'] = np.array(rho, dtype=float)
        except KeyError:
            pass

        try:
            qs = [x.strip() for x in input['qs'].strip('[]').split(',')]
            layers['qs'] = np.array(qs, dtype=float)
        except KeyError:
            pass

        try:
            qp = [x.strip() for x in input['qp'].strip('[]').split(',')]
            layers['qp'] = np.array(qp, dtype=float)
        except KeyError:
            pass

        if 'hisada' in code:
            fqp = [x.strip() for x in input['fqp'].strip('[]').split(',')]
            layers['fqp'] = np.array(fqp, dtype=float)
            fqs = [x.strip() for x in input['fqs'].strip('[]').split(',')]
            layers['fqs'] = np.array(fqs, dtype=float)

    fault['fault_geolocation'] = input['fault_geolocation']
    fault['fault_type'] = input['fault_type']

    if fault['fault_geolocation'] == 'from_hypo' or fault['fault_type'] == 'point' or \
            fault['fault_geolocation'] == 'from_hypo_geometry':
        lon_hypo = float(input['lon_hypo'])
        lat_hypo = float(input['lat_hypo'])
        depth_hypo = float(input['depth_hypo'])
        hypo = create_point(lon_hypo, lat_hypo, depth_hypo)  # in degrees and depth in km
        fault['hypo'] = hypo
        try:
            fault['Ztor'] = float(input['Ztor'])
        except KeyError:
            pass

    if fault['fault_geolocation'] == 'from_geometry' or \
            fault['fault_geolocation'] == 'from_hypo_geometry':
        ul = [x.strip() for x in input['upper-left_vertex'].strip('[]').split(',')]
        ur = [x.strip() for x in input['upper-right_vertex'].strip('[]').split(',')]
        lr = [x.strip() for x in input['lower-right_vertex'].strip('[]').split(',')]
        ll = [x.strip() for x in input['lower-left_vertex'].strip('[]').split(',')]
        fault['min_fault_depth'] = float(input['min_fault_depth'])
        fault['max_fault_depth'] = float(input['max_fault_depth'])
        pbl = Point(ll[0], ll[1], fault['max_fault_depth'])
        pbr = Point(lr[0], lr[1], fault['max_fault_depth'])
        ptl = Point(ul[0], ul[1], fault['min_fault_depth'])
        ptr = Point(ur[0], ur[1], fault['min_fault_depth'])
        vertex = create_vertex(pbl, pbr, ptl, ptr)
        fault['vertex'] = vertex
        fault['Ztor'] = fault['min_fault_depth']

    fault['strike'] = float(input['strike'])
    fault['dip'] = float(input['dip'])
    fault['rake'] = float(input['rake'])
    if fault['rake'] > 180:
        fault['rake'] = fault['rake'] - 360

    if fault['fault_type'] == 'extended':
        try:
            fault['length'] = float(input['length'])
        except KeyError:
            pass

        try:
            fault['width'] = float(input['width'])
        except KeyError:
            pass

        try:
            fault['number_subfaults_strike'] = int(input['number_subfaults_strike'])
        except KeyError:
            pass

        try:
            fault['subfault_length'] = float(input['subfault_length'])
        except KeyError:
            pass

        try:
            fault['number_subfaults_dip'] = int(input['number_subfaults_dip'])
        except KeyError:
            pass

        try:
            fault['subfault_width'] = float(input['subfault_width'])
        except KeyError:
            pass

        #Forse questa parte va aggiustata
        try:
            fault['slip_file'] = input['slip_file']
        except KeyError:
            try:
                fault['slip_mode'] = input['slip_mode']
                if fault['slip_mode'] == 'file_xta':
                    fault['file_xta'] = input['file_xta']
                    fault['nsubcells_length_xta'] = int(input['nsubcells_length_xta'])
                    fault['nsubcells_width_xta'] = int(input['nsubcells_width_xta'])
            except KeyError:
                sys.exit('Error: the slip file must be defined for an extended seismic source')

        #Lo zero della faglia è sempre in alto a sinistra, ma attenzione allo strike!
        try:
            fault['hypo_along_strike'] = float(input['hypo_along_strike'])
        except KeyError:
            pass

        try:
            fault['hypo_down_dip'] = float(input['hypo_down_dip'])
        except KeyError:
            pass

    try:
        fault['IDx'] = input['STF']
    except KeyError:
        pass

    try:
        fault['rupture_velocity'] = float(input['rupture_velocity'])
    except KeyError:
        pass

    try:
        fault['percentage_rupture_velocity'] = float(input['percentage_rupture_velocity'])
    except KeyError:
        pass

    try:
        fault['rise_time'] = float(input['rise_time'])
    except KeyError:
        pass

    try:
        fault['fc'] = float(input['fc'])
    except KeyError:
        pass

    try:
        fault['stress_drop'] = float(input['stress_drop'])
    except KeyError:
        pass
    #stress drop in Mpa

    try:
        fault['Mo'] = float(input['Mo'])
    except KeyError:
        fault['Mw'] = float(input['Mw'])

    if 'ucsb' in code:
        try:
            fault['Tp_Tr'] = float(input['Tp_Tr'])
        except KeyError:
            pass

    computational_param = create_computational_param()
    #computational_param['output_type'] = input['type_output']

    if 'speed' in code:
        computational_param['spectral_degree'] = int(input['spectral_degree'])  # ordine spettrale 4 significa 5 integration
        # points along each direction. Cambia il numero di integration points nella mesh, dentro ogni elemento. 4 divisioni
        computational_param['fval_quality_factor'] = float(input['fval_quality_factor'])
        computational_param['damping'] = int(input['damping'])
        computational_param['tmonitor'] = int(input['tmonitor'])
        computational_param['stoptime'] = float(input['stoptime'])
        mlst = [x.strip() for x in input['mlst'].strip('[]').split(',')]
        computational_param['mlst'] = np.array(mlst, dtype=int)
        computational_param['setuponl'] = input['setuponl']
        computational_param['dt_speed'] = float(input['dt_speed'])

        try:
            computational_param['ndiv'] = int(input['ndiv'])
            # Should be 4 ^ n, where n is a real integer = 0, 1, 2, 3....
            # Each rectangular subfault in Each Fault segment will be divided
            # into ndiv rectangles.This is to increase No.of Point Sources
            # representing Fault Plane
            # Suggested values in ThreePtool 1, 4, 16, 64, 256, 1024
        except KeyError:
            computational_param['ndiv'] = 1

        computational_param['abso_bound_dist_km'] = float(input['distance_absorbing_boundaries'])
        computational_param['fmax_speed'] = float(input['fmax_speed'])  # needed for mesh generation
        computational_param['npplambda'] = float(input['npplambda'])  # pt per lunghezza d'onda (5-10)


    if 'hisada' in code:
        computational_param['npts_hisada'] = int(input['npts_hisada'])
        computational_param['fmax_hisada'] = float(input['fmax_hisada'])
        computational_param['dt_hisada'] = float(input['dt_hisada'])

    if 'ucsb' in code:
        computational_param['kappa'] = float(input['kappa'])
        computational_param['dt_ucsb'] = float(input['dt_ucsb'])
        computational_param['npts_ucsb'] = int(input['npts_ucsb'])
        computational_param['qzero'] = float(input['qzero'])
        computational_param['alpha'] = float(input['alpha'])

    if 'ucsb' in code or fault['slip_mode'] == 'Archuleta':
        computational_param['fmax_ucsb'] = float(input['fmax_ucsb'])
        computational_param['fmin_ucsb'] = float(input['fmin_ucsb'])
        computational_param['realizations'] = int(input['realizations'])
        try:
            seed = [x.strip() for x in input['seed'].strip('[]').split(',')]
            computational_param['seed'] = np.array(seed, dtype=int)
        except KeyError:
            pass

    try:
        receiver_grid_step_km = float(input['site_grid_step_km'])
        try:
            receiver_maximum_dist = input['site_maximum_dist_km']
            if receiver_maximum_dist == 'default':
                distance_max = compute_maximum_distance(fault['Mw'])
            else:
                distance_max = float(receiver_maximum_dist)
            receivers_lon, receivers_lat, receivers_depth, receivers_ID = \
                grid_sites(fault, receiver_grid_step_km, distance_max)
        except KeyError:
            receiver_maximum_dist_N = float(input['site_maximum_dist_km_N'])
            receiver_maximum_dist_S = float(input['site_maximum_dist_km_S'])
            receiver_maximum_dist_E = float(input['site_maximum_dist_km_E'])
            receiver_maximum_dist_W = float(input['site_maximum_dist_km_W'])
            receivers_lon, receivers_lat, receivers_depth, receivers_ID = \
                grid_sites_4lengths(fault, receiver_grid_step_km, receiver_maximum_dist_N, receiver_maximum_dist_S,
                           receiver_maximum_dist_E, receiver_maximum_dist_W)

    except KeyError:
        try:
            receivers_file = input['file_receivers']
            rf = np.loadtxt(receivers_file, skiprows=1)
            receivers_lon = rf[:, 0]
            receivers_lat = rf[:, 1]
            receivers_depth = np.zeros((len(receivers_lat)))
            receivers_ID = []
            for l in range(len(receivers_depth)):
                receivers_ID.append('site' + str(l + 1).zfill(3))
            receivers_ID = np.asarray(receivers_ID)

        except KeyError:
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
                    receivers_ID.append('site' + str(l + 1).zfill(3))
                receivers_ID = np.asarray(receivers_ID)


    output_folder = input['output_folder']

    receivers_X = []
    receivers_Y = []
    for i in range(len(receivers_depth)):
        receivers_utmx, receivers_utmy, zone, letter = determine_utm_coord(receivers_lon[i], receivers_lat[i])
        receivers_X.append(receivers_utmx)
        receivers_Y.append(receivers_utmy)
    sites = create_receivers(receivers_Y, receivers_X, receivers_depth * 1000, receivers_ID,
                             receivers_lon, receivers_lat)

    plot_param = create_plot_param()
    if calculation_mode == '--post' or calculation_mode == '--run' or calculation_mode == '--stitch' or calculation_mode == '--run-nogreen':
        try:
            plot_param['fmin_filter'] = float(input['fmin_filter'])
        except KeyError:
            plot_param['fmin_filter'] = None
        try:
            plot_param['fmax_filter'] = float(input['fmax_filter'])
        except KeyError:
            plot_param['fmax_filter'] = None
        plot_param['time_max_plot'] = float(input['time_max_plot'])

    topo = int(input['topography'])

    #STITCHED
    if calculation_mode == '--stitch':
        computational_param['freq_join'] = float(input['freq_join'])

    cineca = create_cineca_slurm()

    if 'speed' in code:
        cineca['nnodes'] = int(input['nnodes'])
        cineca['memory'] = int(input['memory'])
        cineca['partition'] = input['partition']
        cineca['duration'] = input['duration']
        cineca['ntask_per_node'] = int(input['ntask_per_node'])
        cineca['account'] = input['account']
        cineca['job_name'] = input['job_name']
   
    cineca = create_cineca_slurm()

    return layers, fault, computational_param, sites, plot_param, topo, output_folder, cineca
