def compute_vs_average_on_fault(fault, layers):
    import numpy as np
    #nsuby = fault['number_subfaults_dip']
    #nsubx = fault['number_subfaults_strike']
    nsuby = 100
    nsubx = 200
    vs = layers['vs']
    thickness = list(layers['thk'])
    y_hypc = fault['hypo_down_dip']
    dy = fault['width'] / nsuby
    ncyp = int(y_hypc / dy) + 1
    cyp = (ncyp - 0.5) * dy
    rdip = fault['dip'] * np.pi / 180
    depth_hypc = fault['hypo']['Z']

    nsubxy = nsubx * nsuby
    beta = np.zeros(nsubxy)
    for j in range(nsuby):
        yij = ((j + 1) - 0.5) * dy - cyp
        hk = yij * np.sin(rdip) + depth_hypc
        thickness[len(thickness) - 1] = hk + 0.5
        lyr = -1
        sumh = 0
        while hk > sumh:
            lyr = lyr + 1
            sumh = sumh + thickness[lyr]
        sv = vs[lyr]
        for i in range(nsubx):
            k = ((i + 1) - 1) * nsuby + j
            beta[k] = sv
    vs_average = 0.
    for k in range(nsubxy):
        vs_average = vs_average + beta[k]
    vs_average = vs_average / float(nsubxy)
    return vs_average


def WC1994_dimensions(rake, mag):
    # Wells and Coppersmith 1994
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


def weibull_distribution(scale, shape):
    import numpy as np
    x = np.random.weibull(shape)
    return scale * x


def hypo_strike_Causse_2008():
    import numpy as np
    position_along_length = np.random.normal(loc=0.5, scale=0.23,
                                             size=None)  # Il riferimento indicato nel 3Ptool è (Cotton 2008 BSSA)
    # Xnuc/L normal distribution mu=0.5, sigma=0.23
    if position_along_length < 1. / 10.:
        position_along_length = 1. / 10.
    if position_along_length > 9. / 10.:
        position_along_length = 9. / 10.
    return position_along_length


def hypo_dip_Causse_2008(rake):
    # Il riferimento indicato nel 3Ptool è (Cotton 2008 BSSA) Da controllare perchè nel 3Ptool usano una sola relazione
    if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
        # strike slip
        # Mai et al. (2005) Hypocenter Locations in Finite-Source Rupture Models strike-slip
        # Weibull distribution a=0.626 (scale parameter) b=3.921 Ynuc/L (shape parameter)
        position_along_width = weibull_distribution(0.626, 3.921)
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
        # slip_txt = np.loadtxt(slip_file)
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
    if layers['thk'][len(layers['thk']) - 1] == 0:
        layers['thk'][len(layers['thk']) - 1] = 1000
    for i in range(len(layers['thk'])):
        thk_sum = thk_sum + layers['thk'][i]
        if fault['hypo']['Z'] <= thk_sum:
            nlayer_hypo = i
            return nlayer_hypo
    return


def compute_slip_point_source(layers, fault):
    nlayer_hypo = find_origin_layer(layers, fault)
    mu = layers['rho'][nlayer_hypo] * 1000 * layers['vs'][nlayer_hypo] * 1000 * layers['vs'][nlayer_hypo] * 1000
    area = fault['length'] * 1000 * fault['width'] * 1000
    slip = fault['Mo'] / (mu * area)
    return [slip]


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

def define_missing_parameters(code, layers, fault, computational_param, path_data, topo):
    import sys
    import numpy as np
    import math
    import random
    from rapids.conversions import find_rotated_coords, determine_utm_coord, create_point, create_point_utm
    from rapids.conversions import determine_fault_coordinates_from_hypocentre, create_vertex, create_vertex_utm
    from rapids.conversions import utm_to_lon_lat

    if layers['vel_model'] == 'GNDT_14':
        vel_model_file = path_data + '/VelModel/GNDT/m550014.stp'
        profile = np.loadtxt(vel_model_file, skiprows=1)
        thk = profile[:, 0]
        rho = profile[:, 1]
        vp = profile[:, 2]
        vs = profile[:, 3]
        qp = profile[:,4]
        qs = profile[:,5]
        depth_gndt = 0
        for i in range(len(thk)):
            depth_gndt = depth_gndt + thk[i]
            if depth_gndt > 80: #Deve essere limitato in lunghezza altrimenti da' segmentation fault nel calcolo GF. Ho deciso di prendere il modello fino a circa 80 km di depth (il layer subito dopo 80 km, ma èuna mia scelta)
                layer_number = i
                break

        layers['thk'] = thk[0:layer_number+1]
        layers['rho'] = rho[0:layer_number+1]
        layers['vp'] = vp[0:layer_number+1]
        layers['vs'] = vs[0:layer_number+1]
        layers['qp'] = qp[0:layer_number+1]
        layers['qs'] = qs[0:layer_number+1]

    elif layers['vel_model'] == 'NAC_1D_Friuli':
        vel_model_file = path_data + '/VelModel/friuli_1D.xyz'
        profile = np.loadtxt(vel_model_file, skiprows=3)
        z = profile[:, 0]
        vp = profile[:, 1]
        vs = profile[:, 3]
        rho = profile[:, 5]
        vp_grouped = []
        vs_grouped = []
        rho_grouped = []
        thk_grouped = []
        next_line = 0
        if topo == 1:
            for i in range(len(rho)):
                if not np.isnan(vp[i]):
                    layers['depth_top_layer'] = z[i]
                    break
        else:
            layers['depth_top_layer'] = 0
        for i in range(len(rho)):
            if not np.isnan(vp[i]):
                if i == next_line:
                    vp_grouped.append(vp[i])
                    vs_grouped.append(vs[i])
                    rho_grouped.append(rho[i])
                    thk = 0
                    for ii in range(i, len(rho)):
                        if vp[ii] == vp[i] and vs[ii] == vs[i] and rho[ii] == rho[i]:
                            if ii < len(rho) - 1:
                                thk += z[ii + 1] - z[ii]
                            else:
                                thk += 0
                            istop = ii
                    thk_grouped.append(thk)
                    next_line = istop + 1
            else:
                next_line = i + 1

        layers['vp'] = np.asarray(vp_grouped)
        layers['vs'] = np.asarray(vs_grouped)
        layers['rho'] = np.asarray(rho_grouped)
        layers['thk'] = np.asarray(thk_grouped)

    else:
        if layers['rho'] is None:
            layers['rho'] = np.zeros((len(layers['vp'])))
            for i in range(len(layers['vp'])):
                # “Nafe-Drake curve” (Ludwig et al., 1970).
                # Ludwig, W. J., J. E. Nafe, and C. L. Drake (1970).
                # Seismic refraction, in The Sea, A. E. Maxwell, (Editor) Vol. 4, Wiley-Interscience, New York, 53–84.
                if 1.5 <= layers['vp'][i] <= 8.5:
                    layers['rho'][i] = 1.6612 * layers['vp'][i] - 0.4721 * layers['vp'][i] ** 2 + \
                                   0.0671 * layers['vp'][i] ** 3 - 0.0043 * layers['vp'][i] ** 4 + \
                                   0.000106 * layers['vp'][i] ** 5
                else:
                    sys.exit('Error: TO DO - need to find a regression for rho')

    if layers['qs'] is None:
        layers['qs'] = np.zeros((len(layers['vs'])))
        thk_sum = 0
        for i in range(len(layers['vs'])):
            #if layers['vs'][i] < 0.5:
            #    layers['qs'][i] = 10
                # (based on shallow borehole studies by Tullos and Reid (1969), Hamilton (1972), Gibbs et al.(1994), Liu et al.(1994),
                # and Kudo and Shima (1970)).San Francisco Bay Area specific studies include Gibbs et al.(1994) and Liu et al.(1994).
            #elif 0.5 <= layers['vs'][i] < 1.5:
            #    layers['qs'][i] = 20 * layers['vs'][i]  # vs in km/s (Olsen et al., 2003)
            #else:
            #    layers['qs'][i] = 100 * layers['vs'][i]  # Olsen et al., 2003
            #thk_sum += layers['thk'][i]
            #if thk_sum <= 2:
            #    layers['qs'][i] = 0.05 * layers['vs'][i]*1000
            #else:
            #    layers['qs'][i] = 0.1 * layers['vs'][i]*1000

            if layers['vs'][i] <= 0.3:
                layers['qs'][i] = 13
            elif 0.3 < layers['vs'][i] < 5:
                layers['qs'][i] = -16 + 104.13 * layers['vs'][i] - 25.225 * layers['vs'][i]**2 + 8.2184 * layers['vs'][i]**3
            else:
                sys.exit('Error:vs larger than 5 km')

    if layers['qp'] is None:
        layers['qp'] = np.zeros((len(layers['qs'])))
        for i in range(len(layers['vs'])):
            #layers['qp'][i] = 9 / 4 * layers['qs'][i]  # Lay, T., and Wallace, T. (1995). Modern Global Seimology, Vol. 58.
            layers['qp'][i] = 2 * layers['qs'][i]  # (Day and Bradley, 2001; Graves and Day, 2003, Brocher et al., 2007).
        print('Velocity model: qp computed from qs')


    if fault['Mo'] is None and fault['Mw'] is not None:
        fault['Mo'] = 10. ** (1.5 * (fault['Mw']) + 9.05)  # in Nm Hanks & Kanamori (1979)

    elif fault['Mo'] is not None and fault['Mw'] is None:
        fault['Mw'] = 2 / 3 * math.log10(fault['Mo']*10**7) - 10.7 # Hanamori
    else:
        sys.exit('Error: Mo/Mw not found')

    #if code == 'ucsb' or fault['slip_mode'] == 'Archuleta':
    #    if fault['Mw'] < 6.0:
    #        sys.exit('Error: UCSB method works well with earthquakes with M >6:0. '
    #                 'When the source is small, for example, M 5.5 or less, '
    #                 'there are not enough subfaults in our method to establish the correlation structure '
    #                 'between source parameters. Moreover, the initial distribution for each of the '
    #                 'kinematic parameters was based on much larger simulated earthquakes. '
    #                 'It may be that small events have a narrower distribution of source parameters '
    #                 'for each subfault.')

    if fault['fault_type'] == 'point':
        fault['number_subfaults_strike'] = 1
        fault['number_subfaults_dip'] = 1
        fault['length'] = 0.100  # in km
        fault['width'] = 0.100  # in km
        fault['subfault_length'] = fault['length']
        fault['subfault_width'] = fault['width']
        fault['hypo_along_strike'] = 0.5  # The reference is in the top-left corner of the rupture plane when viewed at an angle of 90° from strike
        fault['hypo_down_dip'] = 0.5
        slip = compute_slip_point_source(layers, fault)
        fault['slip'] = slip


    if fault['fault_type'] == 'point':

        if fault['Ztor'] is None:
            fault['Ztor'] = fault['hypo']['Z'] - fault['hypo_down_dip'] * fault['width'] * np.sin(np.radians(fault['dip']))
            if fault['Ztor'] < 0:
                fault['Ztor'] = 0.
            if fault['Ztor'] < fault['min_fault_depth']:
                fault['Ztor'] = fault['min_fault_depth']
        else:
            fault['hypo_down_dip'] = (fault['hypo']['Z'] - fault['Ztor'])/(fault['width'] * np.sin(np.radians(fault['dip'])))

        fault['max_fault_depth'] = fault['Ztor'] + fault['width'] * np.sin(np.radians(fault['dip']))
        hypo_x, hypo_y, zone, letter = determine_utm_coord(fault['hypo']['lon'], fault['hypo']['lat'])
        hypo_utm = create_point_utm(hypo_y, hypo_x, fault['hypo']['Z'] * 1000)  # in m
        fault['hypo_utm'] = hypo_utm

        pbl, pbr, ptl, ptr = determine_fault_coordinates_from_hypocentre(fault)
        vertex = create_vertex(pbl, pbr, ptl, ptr)
        fault['vertex'] = vertex

    print(fault['fault_type'])

    if fault['fault_type'] == 'extended':

        if fault['fault_geolocation'] == 'from_hypo':
    
            if fault['length'] is None:
                areaWC94, lengthWC94 = WC1994_dimensions(fault['rake'], fault['Mw'])
                if fault['width'] is not None:
                    fault['length'] = areaWC94 / fault['width']
                else:
                    fault['length'] = lengthWC94
            if fault['width'] is None:
                areaWC94, lengthWC94 = WC1994_dimensions(fault['rake'], fault['Mw'])
                fault['width'] = areaWC94 / fault['length']

            if fault['hypo_along_strike'] is None:
                fault['hypo_along_strike'] = hypo_strike_Causse_2008()

            if fault['Ztor'] is None:
                if fault['hypo_down_dip'] is None:
                    fault['hypo_down_dip'] = hypo_dip_Causse_2008(fault['rake'])
                fault['Ztor'] = fault['hypo']['Z'] - fault['hypo_down_dip'] * fault['width'] * np.sin(np.radians(fault['dip']))
                if fault['Ztor'] < 0:
                    fault['Ztor'] = 0.
                if fault['Ztor'] < fault['min_fault_depth']:
                    fault['Ztor'] = fault['min_fault_depth']
            else:
                fault['hypo_down_dip'] = (fault['hypo']['Z'] - fault['Ztor'])/(fault['width'] * np.sin(np.radians(fault['dip'])))

            fault['max_fault_depth'] = fault['Ztor'] + fault['width'] * np.sin(np.radians(fault['dip']))
            hypo_x, hypo_y, zone, letter = determine_utm_coord(fault['hypo']['lon'], fault['hypo']['lat'])
            hypo_utm = create_point_utm(hypo_y, hypo_x, fault['hypo']['Z'] * 1000)  # in m
            fault['hypo_utm'] = hypo_utm

            pbl, pbr, ptl, ptr = determine_fault_coordinates_from_hypocentre(fault)
            vertex = create_vertex(pbl, pbr, ptl, ptr)
            fault['vertex'] = vertex

            pbl_utm_X, pbl_utm_Y, zone, letter = determine_utm_coord(fault['vertex']['pbl']['lon'], fault['vertex']['pbl']['lat'])
            pbr_utm_X, pbr_utm_Y, zone, letter = determine_utm_coord(fault['vertex']['pbr']['lon'], fault['vertex']['pbr']['lat'])
            ptl_utm_X, ptl_utm_Y, zone, letter = determine_utm_coord(fault['vertex']['ptl']['lon'], fault['vertex']['ptl']['lat'])
            ptr_utm_X, ptr_utm_Y, zone, letter = determine_utm_coord(fault['vertex']['ptr']['lon'], fault['vertex']['ptr']['lat'])
            pbl_utm = create_point_utm(pbl_utm_X, pbl_utm_Y, fault['vertex']['pbl']['Z'] * 1000)
            pbr_utm = create_point_utm(pbr_utm_X, pbr_utm_Y, fault['vertex']['pbr']['Z'] * 1000)
            ptl_utm = create_point_utm(ptl_utm_X, ptl_utm_Y, fault['vertex']['ptl']['Z'] * 1000)
            ptr_utm = create_point_utm(ptr_utm_X, ptr_utm_Y, fault['vertex']['ptr']['Z'] * 1000)
            vertex_utm = create_vertex_utm(pbl_utm, pbr_utm, ptl_utm, ptr_utm)
            fault['vertex_utm'] = vertex_utm


        if fault['fault_geolocation'] == 'from_hypo_geometry':
            hypo_x, hypo_y, zone, letter = determine_utm_coord(fault['hypo']['lon'], fault['hypo']['lat'])
            hypo_utm = create_point_utm(hypo_y, hypo_x, fault['hypo']['Z'] * 1000)  # in m
            fault['hypo_utm'] = hypo_utm

            fault['hypo_down_dip'] = (fault['hypo']['Z'] - fault['Ztor'])/(fault['width'] * np.sin(np.radians(fault['dip'])))

            a = np.sqrt((hypo_x-ptl_utm_X)**2+(hypo_y-ptl_utm_Y)**2)
            c = np.sqrt((ptr_utm_X-ptl_utm_X)**2+(ptr_utm_Y-ptl_utm_Y)**2)
            b = np.sqrt((ptr_utm_X-hypo_x)**2+(ptr_utm_Y-hypo_y)**2)
            cos_angle_B = (a**2+c**2-b**2)/(2*a*c)

            fault['hypo_along_strike'] = a*cos_angle_B / c

        if fault['fault_geolocation'] == 'from_geometry':
            if fault['hypo_down_dip'] is None:
                fault['hypo_down_dip'] = hypo_dip_Causse_2008(fault['rake'])
            depth_hypo = fault['Ztor'] + \
                         fault['hypo_down_dip'] * fault['width'] * np.sin(np.radians(fault['dip']))
            if fault['hypo_along_strike'] is None:
                fault['hypo_along_strike'] = hypo_strike_Causse_2008()
            Hypo_strike = fault['hypo_along_strike'] * fault['length'] * 1000
            Hypo_dip = fault['hypo_down_dip'] * fault['width'] * 1000
            # i - sono per la convenzione di SPEED
            HypLoc = [Hypo_strike, -Hypo_dip * np.cos(np.deg2rad(fault['dip'])),
                      -Hypo_dip * np.sin(np.deg2rad(fault['dip']))]
            HypGlo = np.zeros(3)
            HypGlo[0], HypGlo[1] = find_rotated_coords(HypLoc[0], HypLoc[1], -(fault['strike'] - 90))
            HypGlo[0] = HypGlo[0] + ptl_utm_X
            HypGlo[1] = HypGlo[1] + ptl_utm_Y
            HypGlo[2] = -HypLoc[2] + fault['Ztor']
            fault['hypo_utm'] = create_point_utm(HypGlo[1], HypGlo[0], HypGlo[2])
            lon_hypo, lat_hypo = utm_to_lon_lat(fault['hypo_utm']['Y'], fault['hypo_utm']['X'], zone)
            hypo = create_point(lon_hypo, lat_hypo, depth_hypo)  # in degrees and depth in km
            fault['hypo'] = hypo


        if fault['rupture_velocity'] is None:
            vs_average = compute_vs_average_on_fault(fault, layers)
            if fault['percentage_rupture_velocity'] is not None:
                fault['rupture_velocity'] = fault['percentage_rupture_velocity'] * vs_average
            else:
                # span the range of rupture velocities expected for most earthquakes (e.g. Heaton 1990)
                # Heaton T. Evidence for and implications of self-healing pulses of slip in earthquake rupture,
                # Phys. Earth planet. Inter., 1990, vol. 64 (pg. 1-20)
                fault['rupture_velocity'] = np.random.uniform(0.65, 0.85) * vs_average
                # rupture velocity values reported in source studies mainly range between 0.65Vs and 0.85Vs [e.g.,
                # Heaton, 1990].

        if fault['slip_mode'] == 'uniform':
            # TODO con Wells and Coppersmith e vedere anche cosa fare con dimensioni sottosorgenti e SPEED
            slip = compute_slip_point_source(layers, fault)
            print('slip')

        if fault['slip_mode'] == 'Archuleta' or fault['slip_mode'] == 'uniform':

            if fault['subfault_length'] is None and fault['number_subfaults_strike'] is None:
                fault['subfault_length'] = fault['rupture_velocity']/computational_param['fmax_ucsb']
            if fault['subfault_width'] is None and fault['number_subfaults_dip'] is None:
                fault['subfault_width'] = fault['rupture_velocity']/computational_param['fmax_ucsb']

            if fault['subfault_length'] is None and fault['number_subfaults_strike'] is not None:
                fault['subfault_length'] = fault['length'] / (fault['number_subfaults_strike'] - 1) #il -1 aggiunto dal file con il plot
            elif fault['subfault_length'] is not None and fault['number_subfaults_strike'] is None:
                fault['number_subfaults_strike'] = math.ceil(fault['length'] / fault['subfault_length'] + 1)
            else:
                sys.exit('Error: number_subfaults_strike/subfault_length not defined')

            if fault['subfault_width'] is None and fault['number_subfaults_dip'] is not None:
                fault['subfault_width'] = fault['width'] / (fault['number_subfaults_dip'] - 1) #il -1 aggiunto dal file con plot
            elif fault['subfault_width'] is not None and fault['number_subfaults_dip'] is None:
                fault['number_subfaults_dip'] = math.ceil(fault['width'] / fault['subfault_width'] + 1)
            else:
                sys.exit('Error: number_subfaults_dip/subfault_width not defined')

            fault['number_subfaults_strike'] = \
                int(2 ** math.ceil(math.log(fault['length'] / fault['subfault_length'] + 1, 2)))
            fault['number_subfaults_dip'] = \
                int(2 ** math.ceil(math.log(fault['width'] / fault['subfault_width'] + 1, 2)))

            fault['subfault_length'] = fault['length'] / (fault['number_subfaults_strike'] - 1)
            fault['subfault_width'] = fault['width'] / (fault['number_subfaults_dip'] - 1)


    if fault['IDx'] is None:
        fault['IDx'] = 'Archuleta'

    #Compute rupture_duration
    ydlnuc = fault['hypo_down_dip'] - 0.5
    xdlnuc = fault['hypo_along_strike'] - 0.5
    # ovvero la distanza massima
    Distrup = np.sqrt(((0.5 + np.abs(xdlnuc)) * fault['length']) ** 2 +
                              ((0.5 + np.abs(ydlnuc)) * fault['width']) ** 2)
    if fault['fault_type'] == 'extended':
        fault['rupture_duration'] = Distrup/fault['rupture_velocity']

    if fault['rise_time'] is None:
        if fault['fault_type'] == 'point':
            fault['rise_time'] = 10 ** (0.5 * fault['Mw'] - 3.34)  # Somerville et al. (1999)
        # Paul Somerville, Kojiro Irikura, Robert Graves, Sumio Sawada, David Wald,
        # Norman Abrahamson, Yoshinori Iwasaki, Takao Kagawa, Nancy Smith, Akira Kowada;
        # Characterizing Crustal Earthquake Slip Models for the Prediction of Strong Ground Motion.
        # Seismological Research Letters 1999;; 70 (1): 59–80. doi: https://doi.org/10.1785/gssrl.70.1.59
        else:
            #Uso la notazione e formule del pulsyn
            C_heaton = 0.125 #è un parametro di ingresso (la larghezza della striscia che sta nucleando in termini di frazione di Distrup)
            fault['rise_time'] = C_heaton * fault['rupture_duration']    #local slip / rise / HF radiation time:  we believe that ideally m1 - rise = Trise / 2

    if fault['fc'] is None and fault['stress_drop'] is not None:
        # fc da Allmann e Shearer 2009
        # Allmann, B. P. and P. M. Shearer (2009). Global variations of stress drop for moderate to large
        # earthquakes, J. Geophys. Res. 114, B01310, doi:10.1029/2008JB005821, 22 pp.
        vs_average = compute_vs_average_on_fault(fault, layers)
        #controllare unita di misura
        fault['fc'] = (0.42 * vs_average * 10 ** 3) * (fault['stress_drop'] * 10 ** 6 / fault['Mo']) ** (1. / 3.)
    elif fault['fc'] is None and fault['stress_drop'] is None:
        average_stress_drop = 4  # Allman and Shearer
        vs_average = compute_vs_average_on_fault(fault, layers)
        fault['fc'] = (0.42 * vs_average * 10 ** 3) * (average_stress_drop * 10 ** 6 / fault['Mo']) ** (1. / 3.)
        fault['Mo'] = 10. ** (1.5 * fault['Mw'] + 9.05)  # in Nm Hanks & Kanamori (1979)
    else:
        pass

    if 'hisada' in code:
        origin_point = fault['vertex']['pbl']
        origin_utmx, origin_utmy, zone, letter = determine_utm_coord(origin_point['lon'], origin_point['lat'])
        origin = create_point_utm(origin_utmy, origin_utmx, origin_point['Z'] * 1000)
        fault['origin'] = origin

    if 'ucsb' in code or fault['slip_mode'] == 'Archuleta':
        if computational_param['seed'] is None:
            seed1 = (-1) ** random.randint(1, 2) * random.randint(1, 1000)
            seed2 = (-1) ** random.randint(1, 2) * random.randint(1, 1000)
            seed3 = (-1) ** random.randint(1, 2) * random.randint(1, 1000)
            computational_param['seed'] = [seed1, seed2, seed3]

    if 'speed' in code:
        # computational_param['optiout'][0] = 0
        computational_param['optiout'][0] = 1 #Ho messo di default SPEED in spostamento
        computational_param['optiout'][1] = 0
        computational_param['optiout'][2] = 0
        computational_param['optiout'][3] = 0
        computational_param['optiout'][4] = 0
        computational_param['optiout'][5] = 0
        # if computational_param['output_type'] == 'dis':
        #     computational_param['optiout'][0] = 1
        # if computational_param['output_type'] == 'vel':
        #     computational_param['optiout'][1] = 1
        # if computational_param['output_type'] == 'acc':
        #     computational_param['optiout'][2] = 1

    return fault, layers
