from scipy.interpolate import RectBivariateSpline, bisplrep

def nearest_power_of_3(x):
    import numpy as np
    """Arrotonda x alla potenza di 3 più vicina e restituisce l'esponente."""
    exp = round(np.log(x) / np.log(3))  # Trova l'esponente più vicino
    return 3 ** exp, exp  # Restituisce la potenza di 3 e l'esponente


def redefine_layers(layers_orig, fault, computational_param, topo):
    import numpy as np
    from rapids.read_input_data import create_material_properties

    layers_new = create_material_properties()

    depth_max_layers = fault['Ztor'] + fault['width'] * np.sin(np.radians(fault['dip'])) + \
                       computational_param['abso_bound_dist_km']


    vp_newlayers = []
    vs_newlayers = []
    thk_newlayers = []
    rho_newlayers = []
    qp_newlayers = []
    qs_newlayers = []
    h = np.cumsum(layers_orig['thk'])
    for ipt in range(len(layers_orig['thk'])):
        if h[ipt] <= depth_max_layers:
            index = ipt
            vp_newlayers.append(layers_orig['vp'][ipt])
            vs_newlayers.append(layers_orig['vs'][ipt])
            rho_newlayers.append(layers_orig['rho'][ipt])
            qp_newlayers.append(layers_orig['qp'][ipt])
            qs_newlayers.append(layers_orig['qs'][ipt])
            thk_newlayers.append(layers_orig['thk'][ipt])
    if depth_max_layers < h[index + 1]:
        if depth_max_layers - h[index] > 1000: #deve essere maggiore di 1 km altrimenti non lo considera
            vp_newlayers.append(layers_orig['vp'][index + 1])
            vs_newlayers.append(layers_orig['vs'][index + 1])
            rho_newlayers.append(layers_orig['rho'][index + 1])
            qp_newlayers.append(layers_orig['qp'][index + 1])
            qs_newlayers.append(layers_orig['qs'][index + 1])
            thk_newlayers.append(depth_max_layers - h[index])

    layers_new['vp'] = np.asarray(vp_newlayers)
    layers_new['vs'] = np.asarray(vs_newlayers)
    layers_new['thk'] = np.asarray(thk_newlayers)
    layers_new['qp'] = np.asarray(qp_newlayers)
    layers_new['qs'] = np.asarray(qs_newlayers)
    layers_new['rho'] = np.asarray(rho_newlayers)
    layers_new['depth_top_layer'] = layers_orig['depth_top_layer']
    layers_new['thk_topo'] = layers_orig['thk_topo']
    if topo == 'yes' and layers_new['depth_top_layer']<0:
        layers_new['thk'][0] = layers_new['thk_topo']

    return layers_new


def slip_function_Archuleta(time, trupt, trise, dt, fault):
    import sys
    import numpy as np
    if time <= trupt:
        slip = 0
    elif time >= trupt + trise:
        slip = 1
    else:
        Tr = trise
        Te = (1 - fault['Tp_Tr']) * Tr
        Tp = fault['Tp_Tr'] * Tr

        np0 = int(Tp / dt + 1)
        np1 = int(Te / dt + 1)
        np2 = int(Tr / dt + 1)
        np_current = int((time - trupt) / dt + 1)

        if np_current > np2:
            sys.exit('Error in Archuleta function')

        psv = np.sqrt(1 + 100 / (np0 * dt))
        sum1 = 0
        svf = np.zeros(np2)

        for i in range(1, np0 + 1):
            ts = (i - 1) * dt
            svi = ts * psv / Tp * np.sin(0.5 * np.pi / Tp * ts)
            svf[i - 1] = svi
            sum1 += svi

        for i in range(np0 + 1, np1 + 1):
            ts = (i - 1) * dt
            svi = np.sqrt(1. + 100. / ts)
            svf[i - 1] = svi
            sum1 += svi

        for i in range(np1 + 1, np2 + 1):
            ts = (i - 1) * dt
            svi = np.sqrt(1. + 100. / ts) * np.sin((np2 - i) * dt * np.pi * 0.5 / (Tr - Te))
            svf[i - 1] = svi
            sum1 += svi

        sum1 *= dt
        svf /= sum1

        int_svf = np.zeros(np2)
        integ_svf = np.zeros(np2)

        for i in range(1, np2):
            int_svf[i] = 0.5 * (svf[i] + svf[i - 1]) * dt
            integ_svf[i] = np.sum(int_svf[1:i + 1])

        slip = integ_svf[np_current - 1]

        if slip >= 1.05:
            sys.exit('Error in Archuleta slip function')

    return slip


def slip_function_YoffeDCF(time, trupt, trise, dt, fault):
    slip = 0
    return slip


def compute_STF(folder, fault):
    import numpy as np
    dt = 0.001
    trise = 1
    tend = 1
    trupt = 0
    time = np.arange(0, tend + dt, dt)
    npts = len(time)
    slip = []
    time = []
    for i in range(npts):
        time.append(i * dt)
        if fault['IDx'] == 'Archuleta':
            slip.append(slip_function_Archuleta(i * dt, trupt, trise, dt, fault))
        if fault['IDx'] == 'Yoffe-DCF':
            slip.append(slip_function_YoffeDCF(i * dt, trupt, trise, dt, fault))
    fileSTF = 'SrcTimeFunc_STD.txt'
    fid = open(folder + '/' + fileSTF, 'w')
    fid.write('{}\n'.format(npts))
    for i in range(npts):
        fid.write('{:13.6f} {:13.6f}\n'.format(time[i], slip[i]))
    fid.close()
    return npts, fileSTF


def create_griddata3ptool(x_ref, y_ref, x_orig, y_orig, x_glo, y_glo, z_glo):
    griddata3ptool = {
        "x_ref": x_ref,
        "y_ref": y_ref,
        "x_orig": x_orig,
        "y_orig": y_orig,
        "x_glo": x_glo,
        "y_glo": y_glo,
        "z_glo": z_glo
    }
    return griddata3ptool


def create_segdata(x_glo, y_glo, z_glo, Rho, Vs, Vp, Mo, slip_ref, rake_ref, tau_ref, trup_ref, vr_ref):
    segdata = {
        "x": x_glo.flatten('F'),
        "y": y_glo.flatten('F'),
        "z": z_glo.flatten('F'),
        "Rho": Rho,
        "Vs": Vs,
        "Vp": Vp,
        "Mo": Mo,
        "slip": slip_ref.flatten('F'),
        "slip_ref": slip_ref,
        "rake": rake_ref,
        "tau": tau_ref,
        "trup": trup_ref,
        "vr": vr_ref
    }
    return segdata


def create_vertdata(xvert, yvert, zvert):
    vertdata = {
        "xvert": xvert,
        "yvert": yvert,
        "zvert": zvert
    }
    return vertdata


def nextpow2(x):
    import numpy as np
    from numpy import ceil, log2
    """returns the smallest power of two that is greater than or equal to the
    absolute value of x.

    This function is useful for optimizing FFT operations, which are
    most efficient when sequence length is an exact power of two.

    """
    res = ceil(log2(x))
    return res.astype('int')  # we want integer values only but ceil gives float


def read_srcmodfile_multiseg(file_name):
    import numpy as np
    import itertools
    # Getting Dimensions of File
    nseg = 0
    iline = -1
    startend = np.zeros((10, 2))  # TODO: now it is initialized at 10 segments. This can be adjusted
    seg_tag = []
    with open(file_name, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            iline = iline + 1
            tline = line.strip()
            if tline[0:3] == "SEG":
                nseg = nseg + 1
                startend[nseg - 1, 0] = iline
                key, value = tline.split('SEG', 1)
                seg_tag.append(float(value))
                # hdr = f.readline().strip().split()
                # ncol = len(hdr)
            elif tline == '' and nseg > 0:
                startend[nseg - 1, 1] = iline

    if startend[nseg, 1] == 0:
        startend[nseg, 1] = iline + 1

    srcmoddata = []
    for iseg in range(nseg):
        file = open(file_name)
        content = file.readlines()
        nlines = int(startend[iseg, 1] - startend[iseg, 0] - 1)
        lines = itertools.islice(open(file_name), int(startend[iseg, 0] + 1), nlines + int(startend[iseg, 0]))
        srcmoddata.append(np.loadtxt(lines))
        # srcmoddata.append(np.loadtxt(file_name, max_rows=nlines, skiprows=int(startend[iseg, 0] + 1)))
    return srcmoddata


def define_param_sism_commands(ndiv, folder, fault):
    # This script generates a slip distribution for a given srcmod file with multiple
    # Fault segments.In SRCMOD each segment of fault plane is divided into Rectangular
    # subfaults, and slip value at each rectangular subfault is given.Here each rectangular
    # subfault is subdivided into more partitions, so that we can have more no.of point
    # sources representing fault plane j
    #
    # INPUTS:
    # ndiv = Should be 4^n, where n is a real integer = 0, 1, 2, 3....
    # Each rectangular subfault in Each Fault segment will be divided into ndiv rectangles.This is to
    # increase No.of Point Sources representing Fault Plane
    import numpy as np
    import sys
    from scipy.interpolate import griddata, interp2d
    from rapids.conversions import find_rotated_coords
    import scipy

    # NB dictionaries inizializzati con max 10 segmenti!! da sistemare

    fault_dat = np.loadtxt(folder + '/fault.dat', skiprows=2)
    if fault_dat.ndim == 1:
        fault_dat_test = np.zeros((1, np.size(fault_dat)))
        fault_dat_test[0, :] = fault_dat
        fault_dat = fault_dat_test

    nseg = np.size(fault_dat, 0)  # Fault NUM and also IDROW. IDCOL is fault_dat[1], which is the id number
    segL = fault_dat[:, 2]  # in m
    segW = fault_dat[:, 3]  # in m
    segStk = fault_dat[:, 4]
    segDip = fault_dat[:, 5]
    faultF1 = np.zeros((nseg, 3))
    faultF2 = np.zeros((nseg, 3))
    faultF3 = np.zeros((nseg, 3))
    faultF4 = np.zeros((nseg, 3))
    for i in range(nseg):
        faultF1[i, :] = [fault_dat[i, 6], fault_dat[i, 7], fault_dat[i, 8]]
        faultF2[i, :] = [fault_dat[i, 9], fault_dat[i, 10], fault_dat[i, 11]]
        faultF3[i, :] = [fault_dat[i, 12], fault_dat[i, 13], fault_dat[i, 14]]
        faultF4[i, :] = [fault_dat[i, 15], fault_dat[i, 16], fault_dat[i, 17]]

    # Hypocenter
    Hypo_seg = 1  # maybe to be defined by the user. Questo significa che l'ipocentro è nel primo ed unico segmento!!!
    Hypo_strike = fault['hypo_along_strike'] * segL  # Al momento questo funziona solo se c'è un solo segmento!!!!
    # Da verificare e sistemare quando ci sono più di un segmento!!!!!!!
    Hypo_dip = fault['hypo_down_dip'] * segW

    mech_prop = np.loadtxt(folder + '/mech_prop.dat', skiprows=1)
    if mech_prop.ndim == 1:
        mech_prop = np.reshape(mech_prop, (1, mech_prop.size))
    mech_prop[:, 0:4] = mech_prop[:, 0:4] * 1000
    lay_bot_surf_z = -np.cumsum(mech_prop[:, 3])

    # Slip
    slip_sism = read_srcmodfile_multiseg(folder + '/SLIP.srcmod')

    if len(slip_sism) != nseg:
        sys.exit('No.of Segments mismatch in slip SRCMOD file')

    n = 10  # Example size of the vector
    griddata3ptool = {}
    segdata = {}
    vertdata = {}
    for i in range(n):
        key = i
        griddata3ptool[key] = None
        segdata[key] = None
        vertdata[key] = None

    # Finding Grid Points -> coordinates of centroids of Rectangular Subfaults
    # in Each Fault Segment

    dx_ref = np.zeros((nseg, 1))
    dy_ref = np.zeros((nseg, 1))
    for iseg in range(nseg):
        nDip = np.size(slip_sism[iseg], 0)
        nStk = np.size(slip_sism[iseg], 1)
        dx_orig = segL[iseg] / nStk
        dy_orig = segW[iseg] / nDip
        x_orig, y_orig = np.meshgrid(np.arange(dx_orig / 2, segL[iseg], dx_orig),
                                     np.arange(dy_orig / 2, segW[iseg], dy_orig))
        slip_orig = slip_sism[iseg]

        # Making Finer grid and Interpolating Slip
        # ndiv Should be 4^n, where n is a real integer = 0,1,2,3....
        # Each rectangular subfault in Each Fault segment will be divided
        # into ndiv rectangles. This is to increase No. of Point Sources
        # representing Fault Plane
        if np.remainder(nextpow2(ndiv), 2) == 0:
            if ndiv == 1:
                npart = 1
            else:
                npart = int(np.sqrt(ndiv))
        else:
            sys.exit('Give correct value for ndiv')

        dx_ref[iseg] = dx_orig / npart
        dy_ref[iseg] = dy_orig / npart
        x_ref_linear = np.arange(dx_ref[iseg] / 2, segL[iseg], dx_ref[iseg])
        y_ref_linear = np.arange(dy_ref[iseg] / 2, segW[iseg], dy_ref[iseg])
        x_ref, y_ref = np.meshgrid(x_ref_linear, y_ref_linear)
        arr1 = np.arange(dx_orig / 2, segL[iseg], dx_orig)
        arr1 = np.concatenate([[0], arr1, [segL[iseg]]])
        arr2 = np.arange(dy_orig / 2, segW[iseg], dy_orig)
        arr2 = np.concatenate([[0], arr2, [segW[iseg]]])
        x_orig1, y_orig1 = np.meshgrid(arr1, arr2)
        slip_orig1 = np.zeros((np.size(x_orig1, 0), np.size(x_orig1, 1)))
        slip_orig1[1:-1, 1:-1] = slip_orig
        slip_orig1[:, -1] = slip_orig1[:, -2]
        slip_orig1[:, 0] = slip_orig1[:, 1]
        slip_orig1[0, :] = slip_orig1[1, :]
        slip_orig1[-1, :] = slip_orig1[-2, :]

        if npart == 1:
            slip_ref = slip_orig
        else:
            f_int = interp2d(arr1, arr2, slip_orig1, kind='linear')
            slip_ref = f_int(x_ref_linear, y_ref_linear)

            # import matplotlib.pyplot as plt
            # fig, ax = plt.subplots(nrows=1, ncols=2)
            # ax[0].pcolormesh(x_orig1, y_orig1, slip_orig1)
            # ax[1].pcolormesh(x_ref, y_ref, slip_ref)
            # plt.show()

            if iseg + 1 == 7:
                for j in range(np.size(x_orig, 1)):
                    for i in range(np.size(x_orig, 0)):
                        indr1 = i * npart + 1
                        indc1 = j * npart + 1
                        slip_ref[indr1: (npart * (i + 1)), indc1: (npart * (j + 1))] = \
                            slip_orig[i, j]  # to check if index are ok in this formulation

        # Local Coordinate system is a Right hand Cartesian Coordinate System
        # +x_loc = Strike Direction
        # +z_loc = Up direction
        # +y_loc = following Right hand coordinate rule
        # Origin = Top Left Vertex of Fault Segment
        x_loc = x_ref
        y_loc = -y_ref * np.cos(np.deg2rad(segDip[iseg]))
        z_loc = -y_ref * np.sin(np.deg2rad(segDip[iseg]))

        # Global Coordinate system is a Right hand Cartesian Coorfinate System
        # +x_glo = East Direction
        # +z_glo = Up direction
        # +y_glo = North Direction
        # Origin = As per User specified when defining Fault segment Vertices in fault.dat
        # Rotating Local coordinate to Global coordinate system
        x_glo, y_glo = find_rotated_coords(x_loc, y_loc, -(segStk[iseg] - 90))

        # Translating the Origin to Global reference Coordinate System
        x_glo = x_glo + faultF1[iseg, 0]
        y_glo = y_glo + faultF1[iseg, 1]
        z_glo = z_loc + faultF1[iseg, 2]

        # hypo
        if iseg + 1 == Hypo_seg:
            HypGlo = np.zeros((3))
            HypLoc = [Hypo_strike, -Hypo_dip * np.cos(np.deg2rad(segDip[iseg])),
                      -Hypo_dip * np.sin(np.deg2rad(segDip[iseg]))]
            HypGlo[0], HypGlo[1] = find_rotated_coords(HypLoc[0], HypLoc[1], -(segStk[iseg] - 90))
            HypGlo[0] = HypGlo[0] + faultF1[iseg, 0]
            HypGlo[1] = HypGlo[1] + faultF1[iseg, 1]
            HypGlo[2] = HypLoc[2] + faultF1[iseg, 2]

        # Material Property at each point
        lay_num = np.zeros((len(slip_ref.flatten('F'))))
        for ipt in range(len(slip_ref.flatten('F'))):
            indx1 = np.where(lay_bot_surf_z <= z_glo.flatten('F')[ipt])
            if not indx1:
                lay_num[ipt] = len(lay_bot_surf_z) - 1
            else:
                lay_num[ipt] = indx1[0][0]
        lay_num = np.asarray(lay_num, dtype=int)
        mu = mech_prop[:, 2] * mech_prop[:, 1] * mech_prop[:, 1]
        Rho = np.zeros((len(lay_num)))
        Vs = np.zeros((len(lay_num)))
        Vp = np.zeros((len(lay_num)))
        Mo = np.zeros((len(lay_num)))
        for ipt in range(len(lay_num)):
            Rho[ipt] = mech_prop[lay_num[ipt], 2]
            Vs[ipt] = mech_prop[lay_num[ipt], 1]
            Vp[ipt] = mech_prop[lay_num[ipt], 0]
            Mo[ipt] = mu[lay_num[ipt]] * slip_ref.flatten('F')[ipt] * (dx_ref[iseg] * dy_ref[iseg])

        # Saving data
        griddata3ptool[iseg] = create_griddata3ptool(x_ref, y_ref, x_orig, y_orig, x_glo, y_glo, z_glo)
        rake_ref = None
        tau_ref = None
        trup_ref = None
        vr_ref = None
        segdata[iseg] = create_segdata(x_glo, y_glo, z_glo, Rho, Vs, Vp, Mo, slip_ref,
                                       rake_ref, tau_ref, trup_ref, vr_ref)

        # 4 vertices of Each subfault
        eps = 0.00001
        xv_loc, yv_loc = np.meshgrid(np.arange(0, segL[iseg] + eps, dx_ref[iseg]),
                                     np.arange(0, segW[iseg] + eps, dy_ref[iseg]))

        xv_loc = xv_loc
        zv_loc = -yv_loc * np.sin(np.deg2rad(segDip[iseg]))
        yv_loc = -yv_loc * np.cos(np.deg2rad(segDip[iseg]))
        xv_glo, yv_glo = find_rotated_coords(xv_loc, yv_loc, -(segStk[iseg] - 90))

        # Translating the Origin to Global reference Coordinate System
        xv_glo = xv_glo + faultF1[iseg, 0]
        yv_glo = yv_glo + faultF1[iseg, 1]
        zv_glo = zv_loc + faultF1[iseg, 2]

        ipt = -1
        xvert = np.zeros((x_ref.shape[1] * x_ref.shape[0], 4))
        yvert = np.zeros((x_ref.shape[1] * x_ref.shape[0], 4))
        zvert = np.zeros((x_ref.shape[1] * x_ref.shape[0], 4))
        for j in range(x_ref.shape[1]):
            for i in range(x_ref.shape[0]):
                ipt += 1
                vrow = [i, i, i + 1, i + 1]
                vcol = [j, j + 1, j + 1, j]
                for i1 in range(4):
                    xvert[ipt, i1] = xv_glo[vrow[i1], vcol[i1]]
                    yvert[ipt, i1] = yv_glo[vrow[i1], vcol[i1]]
                    zvert[ipt, i1] = zv_glo[vrow[i1], vcol[i1]]

        # Cross Verify
        tol = 5  # Tolerance is 5 m
        dist_diff = np.abs(xv_glo[-1, -1] - faultF3[iseg, 0]) ** 2
        dist_diff = dist_diff + np.abs(yv_glo[-1, -1] - faultF3[iseg, 1]) ** 2
        dist_diff = dist_diff + np.abs(zv_glo[-1, -1] - faultF3[iseg, 2]) ** 2
        dist_diff = np.sqrt(dist_diff)
        if dist_diff > tol:
            print('dist_diff = ', dist_diff)
            sys.exit('Tolerance is high. Check Coordinates of Fault Vertices in fault.dat file')

        vertdata[iseg] = create_vertdata(xvert, yvert, zvert)

    # Seismic Moment
    Mo_calc = 0
    for iseg in range(nseg):
        Mo_calc = Mo_calc + np.sum(segdata[iseg]['Mo'])
    Mw_calc = 2 / 3 * np.log10(Mo_calc) - 6.033
    print('Target Mw = ', fault['Mw'])
    print('Calculated Mw =', Mw_calc)

    # Reading other *.SRCMOD Files

    # Trup
    trup = read_srcmodfile_multiseg(folder + '/TRUP.srcmod')
    if len(trup) != nseg:
        sys.exit('No.of Segments mismatch in trup SRCMOD file')

    # Rise
    rise = read_srcmodfile_multiseg(folder + '/RISE.srcmod')
    if len(rise) != nseg:
        sys.exit('No.of Segments mismatch in rise SRCMOD file')

    # Rake
    rake = read_srcmodfile_multiseg(folder + '/RAKE.srcmod')
    if len(rake) != nseg:
        sys.exit('No.of Segments mismatch in rake SRCMOD file')

    for iseg in range(nseg):
        x_orig = griddata3ptool[iseg]['x_orig']
        y_orig = griddata3ptool[iseg]['y_orig']
        x_ref = griddata3ptool[iseg]['x_ref']
        y_ref = griddata3ptool[iseg]['y_ref']

        trup_orig = trup[iseg]
        f_int = RectBivariateSpline(y_orig[:, 0], x_orig[0, :], trup_orig)
        trup_ref = f_int(y_ref[:, 0], x_ref[0, :])
        ind2 = np.isnan(trup_ref)
        ind3 = np.where(ind2 == True)
        trup_ref[ind3] = 0

        dist = np.sqrt((griddata3ptool[iseg]['x_glo'] - HypGlo[0]) ** 2 +
                       (griddata3ptool[iseg]['y_glo'] - HypGlo[1]) ** 2 +
                       (griddata3ptool[iseg]['z_glo'] - HypGlo[2]) ** 2)
        vr_ref = dist / trup_ref

        rake_orig = rake[iseg]
        tau_orig = rise[iseg]
        rake_ref = np.zeros((np.size(x_ref, 0), np.size(x_ref, 1)))
        tau_ref = np.zeros((np.size(x_ref, 0), np.size(x_ref, 1)))

        for j in range(x_orig.shape[1]):
            for i in range(x_orig.shape[0]):
                indr1 = i * npart
                indc1 = j * npart
                rake_ref[indr1:(npart * (i + 1)) + 1, indc1:(npart * (j + 1)) + 1] = rake_orig[i, j]
                tau_ref[indr1:(npart * (i + 1)) + 1, indc1:(npart * (j + 1)) + 1] = tau_orig[i, j]

        segdata[iseg]['rake'] = rake_ref.flatten('F')
        segdata[iseg]['tau'] = tau_ref.flatten('F')
        segdata[iseg]['trup'] = trup_ref.flatten('F')
        segdata[iseg]['vrup'] = vr_ref.flatten('F')

    return segdata, vertdata, nseg, segStk, segDip, HypGlo


def create_cartesian_vertex_fault(fault):
    # the program outputs the source coordinates(Fault according to the
    # reference system used in Hisada program(Hisada & Bielak, 2003) provided
    # the geometry of the EXTENDED fault
    import numpy as np

    azim = fault['strike']
    azim = azim * np.pi / 180 - np.pi / 2
    dip = fault['dip']
    dip = dip * (-np.pi) / 180
    rake = fault['rake']
    height = fault['width'] * 1000
    width = fault['length'] * 1000
    depth = fault['Ztor'] * 1000

    # 1) defines corners
    # 1 -------------2
    # |              |
    # |              |
    # 4 -------------3
    x = np.zeros(4)
    y = np.zeros(4)
    z = np.zeros(4)

    x[0] = 0
    y[0] = 0
    z[0] = 0
    x[1] = width
    y[1] = 0
    z[1] = 0
    x[2] = width
    y[2] = 0
    z[2] = height
    x[3] = 0
    y[3] = 0
    z[3] = height

    # 2) rotation  according to azimuth(vertical axis through points 1, 4)
    x2 = x[1]
    y2 = y[1]
    x[1] = np.cos(azim) * x2 + np.sin(azim) * y2
    y[1] = -np.sin(azim) * x2 + np.cos(azim) * y2
    x3 = x[2]
    y3 = y[2]
    x[2] = np.cos(azim) * x3 + np.sin(azim) * y3
    y[2] = np.sin(azim) * x3 - np.cos(azim) * y3

    # 3) rotation according to dip(horizontal axis through points 1, 2)
    vx = 0
    vy = 0
    m = np.abs(z[2] - z[1])
    vz = (z[2] - z[1]) / m

    ix = x[1] - x[0]
    iy = y[1] - y[0]
    iz = 0
    n = np.sqrt(ix * ix + iy * iy)
    ix = ix / n
    iy = iy / n

    wx = iy * vz - iz * vy
    wy = iz * vx - ix * vz
    wz = ix * vy - iy * vx
    ux = (np.cos(dip) * wx + np.sin(dip) * vx) * m
    uy = (np.cos(dip) * wy + np.sin(dip) * vy) * m
    uz = (np.cos(dip) * wz + np.sin(dip) * vz) * m
    x[2] = x[1] + ux
    y[2] = y[1] + uy
    z[2] = z[1] + uz
    x[3] = x[0] + ux
    y[3] = y[0] + uy
    z[3] = z[0] + uz

    # i = np.zeros(3)
    # j = np.zeros(3)
    # i[0] = x[1] - x[0]
    # i[1] = y[1] - y[0]
    # i[2] = z[1] - z[0]
    # i = i / np.linalg.norm(i, 2)
    # j[0] = x[3] - x[0]
    # j[1] = y[3] - y[0]
    # j[2] = z[3] - z[0]
    # j = j / np.linalg.norm(j, 2)
    # xh = i * xh + j * zh

    # TRANSLATION of vector
    x0 = fault['vertex_utm']['ptl']['X']
    y0 = fault['vertex_utm']['ptl']['Y']
    x = x + x0
    y = y + y0
    z = z - depth

    return x, y, z


def create_input_extended(folder, fault, layers, folder_ucsb):
    import numpy as np
    from rapids.conversions import convert_ucsb_moment2slip, read_slip_xta

    fid = open(folder + '/mech_prop.dat', 'w')
    fid.write('{}    {}\n'.format(len(layers['thk']), 1.0))
    for i in range(len(layers['thk'])):
        fid.write('{} {} {} {} {} {}\n'.format(layers['vp'][i], layers['vs'][i], layers['rho'][i],
                                               layers['thk'][i], layers['qs'][i], layers['qp'][i]))
    fid.close()

    # N.B. le coordinate UTM che venivano fuori dall'altra procedura divergevano di 50 m e 200 m da quelle calcolate
    # direttamente in coordinate cartesiane e mi davano errore con il 3ptool, per questo ho implementato la subroutine
    # create_cartesian_vertex_fault. Sarebbe allora più corretto definire direttamente le UTM da qua e da qui quelle geogradiche
    # Quindi la routine presa da GEM andrebbe usata solo per definire ptl, che andrebbe poi dato in pasto a create_cartesian_vertex_fault

    x, y, z = create_cartesian_vertex_fault(fault)

    idrow = 1
    idcol = 1
    fid = open(folder + '/fault.dat', 'w')
    fid.write('{}\n'.format('% SOURCE PARAMETERS'))
    fid.write('{}\n'.format('%ID     LF(m)   WF(m) st(°) dp(°) F1_x(m)   F1_y(m)     F1_z(m)     F2_x(m)    F2_y(m)'
                            '      F2_z(m)    F3_x(m)   F3_y(m)      F3_z(m)     F4_x(m)    F4_y(m)     F4_z(m)'))
    fid.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(idrow, idcol, fault['length'] * 1000,
                                                                               fault['width'] * 1000, fault['strike'],
                                                                               fault['dip'],
                                                                               x[0], y[0], z[0], x[1], y[1], z[1],
                                                                               x[2], y[2], z[2], x[3], y[3], z[3]))
    fid.close()

    if fault['slip_mode'] == 'Archuleta':
        if fault['IDx'] == 'Yoffe-DCF':
            file_slip = folder_ucsb + '/Source.bst'
        else:
            file_slip = folder_ucsb + '/Source.001'

        matrix_ucsb = np.loadtxt(file_slip, skiprows=1)
        # 1-3.  are the coordinates in meters (North, East, Down) of each point source,
        # 4-7.  are the moment (N-m), rupture time (s), rise time(s), and another source
        # parameters (any value if not used)
        # 8-10. are the strike, dip, and rake in degree
        x_ps = matrix_ucsb[:, 0]
        y_ps = matrix_ucsb[:, 1]
        z_ps = matrix_ucsb[:, 2]

        moment_ps = matrix_ucsb[:, 3]
        slip_ps = convert_ucsb_moment2slip(moment_ps, z_ps, layers, fault)
        rupture_time_ps = matrix_ucsb[:, 4]
        rise_time_ps = matrix_ucsb[:, 5]
        strike_ps = matrix_ucsb[:, 7]
        dip_ps = matrix_ucsb[:, 8]
        rake_ps = matrix_ucsb[:, 9]

        slip_srcmod = np.zeros((fault['number_subfaults_dip'], fault['number_subfaults_strike']))
        rupt_srcmod = np.zeros((fault['number_subfaults_dip'], fault['number_subfaults_strike']))
        rise_srcmod = np.zeros((fault['number_subfaults_dip'], fault['number_subfaults_strike']))
        rake_srcmod = np.zeros((fault['number_subfaults_dip'], fault['number_subfaults_strike']))

        inc = -1
        for j in range(fault['number_subfaults_dip']):
            m = inc + 1
            for k in range(fault['number_subfaults_strike']):
                slip_srcmod[j, k] = slip_ps[m]
                rupt_srcmod[j, k] = rupture_time_ps[m]
                rise_srcmod[j, k] = rise_time_ps[m]
                rake_srcmod[j, k] = rake_ps[m]
                m = m + fault['number_subfaults_dip']
            inc = inc + 1

    if fault['slip_mode'] == 'Aoudia_et_2000':
        moment_matrix = np.loadtxt('/Users/elisa/Documents/Progetti/SimulazioniFriuli/Frane/angela.amp')
        nx = np.size(moment_matrix, 0)
        ny = np.size(moment_matrix, 1)
        fault['number_subfaults_strike'] = nx
        fault['number_subfaults_dip'] = ny

        moment_matrix = np.transpose(moment_matrix)
        moment_matrix = moment_matrix * 10 ** 20  # dyne*cm
        moment_matrix = moment_matrix * 10 ** (-7)  # Nm

        fault['subfault_length'] = fault['length'] / (fault['number_subfaults_strike'])
        fault['subfault_width'] = fault['width'] / (fault['number_subfaults_dip'])

        moment_ps = []
        z_ps = []
        rake_ps = []
        rupture_time_ps = []
        rise_time_ps = []
        for j in range(fault['number_subfaults_dip']):
            for k in range(fault['number_subfaults_strike']):
                moment_ps.append(moment_matrix[j, k])
                depth_subsource = j * fault['subfault_width'] + fault['subfault_width'] / 2
                z_ps.append(depth_subsource * 1000)  # in m
                dist_x_subsource = k * fault['subfault_length'] + fault['subfault_length'] / 2  # in km
                rake_ps.append(fault['rake'])
                rise_time_ps.append(fault['rise_time'])
                dist_subsource = np.sqrt((fault['hypo_down_dip'] * fault['width'] - depth_subsource) ** 2 + (
                        fault['hypo_along_strike'] * fault['length'] - dist_x_subsource) ** 2)
                trup_subsource = dist_subsource / fault['rupture_velocity']
                rupture_time_ps.append(trup_subsource)
        moment_ps = np.asarray(moment_ps)
        z_ps = np.asarray(z_ps)
        rake_ps = np.asarray(rake_ps)
        rise_time_ps = np.asarray(rise_time_ps)
        rupture_time_ps = np.asarray(rupture_time_ps)
        slip_ps = convert_ucsb_moment2slip(moment_ps, z_ps, layers, fault)
        # plot_slip_from_file(folder, fault, slip_ps, rupture_time_ps)

    if fault['slip_mode'] == 'file_xta':
        slip_ps, rake_ps, rise_time_ps, rupture_time_ps = read_slip_xta(folder, fault, layers)

        slip_srcmod = slip_ps.reshape(fault['number_subfaults_dip'], fault['number_subfaults_strike'])
        rupt_srcmod = rupture_time_ps.reshape(fault['number_subfaults_dip'], fault['number_subfaults_strike'])
        rise_srcmod = rise_time_ps.reshape(fault['number_subfaults_dip'], fault['number_subfaults_strike'])
        rake_srcmod = rake_ps.reshape(fault['number_subfaults_dip'], fault['number_subfaults_strike'])

    #if fault['slip_mode'] == 'constant':



    mat = np.matrix(slip_srcmod)
    with open('SLIP.srcmod', 'w') as f:
        f.write('\n')
        f.write('\n')
        f.write('SEG 1\n')
        for line in mat:
            np.savetxt(f, line)
        f.write('\n')
        f.write('\n')

    mat = np.matrix(rupt_srcmod)
    with open('TRUP.srcmod', 'w') as f:
        f.write('\n')
        f.write('\n')
        f.write('SEG 1\n')
        for line in mat:
            np.savetxt(f, line)
        f.write('\n')
        f.write('\n')

    mat = np.matrix(rake_srcmod)
    with open('RAKE.srcmod', 'w') as f:
        f.write('\n')
        f.write('\n')
        f.write('SEG 1\n')
        for line in mat:
            np.savetxt(f, line)
        f.write('\n')
        f.write('\n')

    mat = np.matrix(rise_srcmod)
    with open('RISE.srcmod', 'w') as f:
        f.write('\n')
        f.write('\n')
        f.write('SEG 1\n')
        for line in mat:
            np.savetxt(f, line)
        f.write('\n')
        f.write('\n')

    return


def compute_slip_vec(strike, rake, dip):
    # calculates the slip vector for a given fault geometry
    import numpy as np
    slip_vec = {
        "X": None,
        "Y": None,
        "Z": None
    }

    # IMPORTANT: the vector is given with respect to the UTM geographic coordinate system(X, Y, Z):
    # X = EAST, Y = NORTH, Z = UP

    strike = strike * np.pi / 180  # measured clockwise from north
    dip = dip * np.pi / 180  # measured down from the horizontal
    rake = rake * np.pi / 180  # measured counterclockwise from horizon. strike direct., as in Aki&Richard

    slip_vec['X'] = +np.cos(rake) * np.sin(strike) - np.sin(rake) * np.cos(dip) * np.cos(strike)
    slip_vec['Y'] = +np.cos(rake) * np.cos(strike) + np.sin(rake) * np.cos(dip) * np.sin(strike)
    slip_vec['Z'] = np.sin(rake) * np.sin(dip)

    return slip_vec


def compute_norm_vec(strike, rake, dip):
    # calculates the  normal fault vector for a given fault geometry
    import numpy as np
    norm_vec = {
        "X": None,
        "Y": None,
        "Z": None
    }

    # IMPORTANT: the vector is given with respect to the UTM geographic coordinate system(X, Y, Z):
    # X = EAST, Y = NORTH, Z = UP

    strike = strike * np.pi / 180  # measured clockwise from north
    dip = dip * np.pi / 180  # measured down from the horizontal

    norm_vec['X'] = +np.sin(dip) * np.cos(strike)
    norm_vec['Y'] = -np.sin(dip) * np.sin(strike)
    norm_vec['Z'] = np.cos(dip)

    if np.abs(norm_vec['X']) <= 1.e-6:
        norm_vec['X'] = 0
    if np.abs(norm_vec['Y']) <= 1.e-6:
        norm_vec['Y'] = 0
    if np.abs(norm_vec['Z']) <= 1.e-6:
        norm_vec['Z'] = 0
    return norm_vec


def listToString(s):
    # initialize an empty string
    str1 = ""
    # traverse in the string
    for ele in s:
        str1 += str(ele) + ' '
    # return string
    return str1


def define_area(fault, sites, computational_param):
    abso_bound_dist = computational_param['abso_bound_dist_km'] * 1000
    for i in range(2):
        if i == 0:
            str_coord = 'X'
            min_val = min(sites['Y'])
            max_val = max(sites['Y'])
            min_val = min(min_val, fault['hypo_utm']['Y'])
            max_val = max(max_val, fault['hypo_utm']['Y'])
        else:
            str_coord = 'Y'
            min_val = min(sites['X'])
            max_val = max(sites['X'])
            min_val = min(min_val, fault['hypo_utm']['X'])
            max_val = max(max_val, fault['hypo_utm']['X'])

        if fault['fault_type'] == 'extended':
            if fault['vertex_utm']['pbl'][str_coord] > max_val:
                max_val = fault['vertex_utm']['pbl'][str_coord]
            if fault['vertex_utm']['pbr'][str_coord] > max_val:
                max_val = fault['vertex_utm']['pbr'][str_coord]
            if fault['vertex_utm']['ptr'][str_coord] > max_val:
                max_val = fault['vertex_utm']['ptr'][str_coord]
            if fault['vertex_utm']['ptl'][str_coord] > max_val:
                max_val = fault['vertex_utm']['ptl'][str_coord]

        if str_coord == 'X':
            maxX = max_val
        else:
            maxY = max_val

        if fault['fault_type'] == 'extended':
            if fault['vertex_utm']['pbl'][str_coord] < min_val:
                min_val = fault['vertex_utm']['pbl'][str_coord]
            if fault['vertex_utm']['pbr'][str_coord] < min_val:
                min_val = fault['vertex_utm']['pbr'][str_coord]
            if fault['vertex_utm']['ptr'][str_coord] < min_val:
                min_val = fault['vertex_utm']['ptr'][str_coord]
            if fault['vertex_utm']['ptl'][str_coord] < min_val:
                min_val = fault['vertex_utm']['ptl'][str_coord]

        if str_coord == 'X':
            minX = min_val
        else:
            minY = min_val

    coord1 = [minX - abso_bound_dist, minY - abso_bound_dist]
    coord2 = [minX - abso_bound_dist, maxY + abso_bound_dist]
    coord3 = [maxX + abso_bound_dist, maxY + abso_bound_dist]
    coord4 = [maxX + abso_bound_dist, minY - abso_bound_dist]

    return coord1, coord2, coord3, coord4


def define_mesh_size(layers, computational_param):
    lunghezza_onda_list = []
    dz_list = []
    lato_spettrale_list = []
    for i in range(len(layers['thk'])):
        lunghezza_onda = layers['vs'][i] * 1000 / computational_param['fmax_speed']
        lunghezza_onda_list.append(lunghezza_onda)
        dz = lunghezza_onda / computational_param['npplambda']
        dz_list.append(dz)
        lato_spettrale = dz * computational_param['spectral_degree']
        lato_spettrale_list.append(lato_spettrale)
        # Il lato spettrale calcolato è il max. Se ne può usare uno più piccolo
    mesh_size = min(lato_spettrale_list)
    return mesh_size




def find_repeat_interval(filename):
    repeating_lines = []
    repeating_x = None
    with open(filename, 'r') as xyz_file:
        for line_number, line in enumerate(xyz_file, start=1):
            columns = line.strip().split()
            x_coordinate = float(columns[0])
            if line_number == 1:
                x_coordinate_test = x_coordinate
            else:
                if x_coordinate == x_coordinate_test:
                    repeating_lines.append(line_number)
                    repeating_x = x_coordinate
    if repeating_x is not None and len(repeating_lines) > 1:
        return repeating_lines[0] - 1
    else:
        return None


def create_read_topo_py(python_script_cubit, topo_utm_file, interval, folder_mesh):
    # file_stl = '/topo.stl'
    file_cub = folder_mesh + '/topo.cub'
    lines = [
        "#!/usr/bin/env python\n",
        "from __future__ import print_function\n",
        "\n",
        "import os\n",
        "import sys\n",
        "import fileinput\n",
        "import string\n",
        "import math\n",
        "\n",
        "import cubit\n",
        "\n",
        "print(sys.path)\n",
        "\n",
        "#############################################################\n",
        "# USER PARAMETERS\n",
        "\n",
        "# topography file, data points per lon-increment\n",
        "inputFile = '" + topo_utm_file + "'\n",
        "\n",
        "# X coordinate in topography file repeats after line\n",
        "nstep = " + str(interval) + "\n",
        "\n",
        "#############################################################\n",
        "\n",
        "# converts xyz to utm\n"
        "#os.system('./convert_lonlat2utm.pl ptopo.mean.xyz 33 > ptopo.mean.utm ')\n",
        "\n"
        "print('#reading from file: ', inputFile)\n",
        "\n",
        "cubit.cmd('reset')\n",
        "cubit.cmd('echo off')\n",
        "cubit.cmd('Graphics Pause')\n"
        "\n"
        "\n"
        "# creates point vertices\n",
        "print('#creating points...')\n",
        "count = 0\n",
        "for line in fileinput.input(inputFile):\n",
        "  count = count + 1\n",
        "  lineitems = line.split()\n",
        "  x = lineitems[0]\n",
        "  y = lineitems[1]\n",
        "  z = lineitems[2]\n",
        "  xyz = str(x) + ' ' + str(y) + ' ' + str(z)\n",
        "  cubit.cmd('create vertex ' + xyz)\n",
        "fileinput.close()\n",
        "print('#done points: ' + str(count))\n",
        "print('')\n",
        #"cubit.cmd('Display')\n",
        "\n",
        "# creates smooth spline curves for surface\n",
        "print('#creating curves...')\n",
        "countcurves = 0\n",
        "for i in range(1, count+1):\n",
        "  if i > 1:\n",
        "    if i % nstep == 0:\n",
        "      countcurves = countcurves + 1\n",
        "      cubit.cmd('create curve spline vertex ' + str(i-nstep+1) + ' to ' + str(i) + ' delete')\n",
        "print('#done curves: '+str(countcurves))\n",
        "print('')\n",
        #"cubit.cmd('Display')\n",
        #"cubit.cmd('pause')\n"
        "print('')\n",
        "# creates surface\n",
        "print('#creating skin surface...')\n",
        "cubit.cmd('create surface skin curve all')\n",
        "\n",
        #"cubit.cmd('Display')\n",
        #"cubit.cmd('pause')\n",
        "\n",
        "# cleans up\n",
        "cubit.cmd('merge all ')\n",
        "cubit.cmd('delete vertex all')\n",
        "cubit.cmd('delete curve all')\n",
        "\n",
        "print('#done cleaning up')\n",
        #"cubit.cmd('Display')\n",
        #"cubit.cmd('echo on')\n",
        "\n",
        "# saves and exports surface\n",
        "# cubit file (uses ACIS file format)\n",
        "cubit.cmd('save as \"" + file_cub + "\" overwrite')\n",
        #       "# export surface as STL file\n",
        #       "cubit.cmd('export stl ascii \""+file_stl+"\" overwrite')\n",
        "\n",
        "print('#exporting done')\n",
        "print('#finished')\n",
        "cubit.cmd('exit')"
    ]

    with open(python_script_cubit, "w") as file:
        file.writelines(lines)
    return


def rename_vertices(val):
    if val == 1:
        num = 4
    elif val == 2:
        num = 2
    elif val == 3:
        num = 1
    elif val == 4:
        num = 3
    else:
        num = val
    return num


def rename_curves(val):
    if val == 1:
        num = 4
    elif val == 2:
        num = 1
    elif val == 3:
        num = 2
    elif val == 4:
        num = 3
    else:
        num = val
    return num


def create_mesh(folder, computational_param, layers, fault, sites, topo, path_cubit, path_data):
    import math
    import os
    import sys
    from rapids.conversions import determine_utm_coord, utm_to_lon_lat

    folder_mesh = folder + '/MESH'
    command_make_folder_MESH = 'mkdir -p ' + folder_mesh
    os.system(command_make_folder_MESH)
    #I need to remove this since I have included it in the computation of LS
    coord1, coord2, coord3, coord4 = define_area(fault, sites, computational_param)

    if topo == 'yes':
        path_topo = path_data+'/topo'
        temp1, temp2, zone, letter = determine_utm_coord(fault['hypo']['lon'], fault['hypo']['lat'])
#        #coord1_lon, coord1_lat = utm_to_lon_lat(coord1[0], coord1[1], zone)
#        #coord2_lon, coord2_lat = utm_to_lon_lat(coord2[0], coord2[1], zone)
#        #coord3_lon, coord3_lat = utm_to_lon_lat(coord3[0], coord3[1], zone)
#        #coord4_lon, coord4_lat = utm_to_lon_lat(coord4[0], coord4[1], zone)
#        #buffer = 0.2
#        #minlon = min(coord1_lon, coord2_lon, coord3_lon, coord4_lon) - buffer
#        #maxlon = max(coord1_lon, coord2_lon, coord3_lon, coord4_lon) + buffer
#        #minlat = min(coord1_lat, coord2_lat, coord3_lat, coord4_lat) - buffer
#        #maxlat = max(coord1_lat, coord2_lat, coord3_lat, coord4_lat) + buffer
#        #topo_xyz_file = folder_mesh + '/ptopo.mean.xyz'
#        #command_extract = 'gmt blockmean ' + path_topo + '/bedrock.xyz -R' + str(minlon) + '/' + str(
#        #    maxlon) + '/' + str(minlat) \
#        #                  + '/' + str(maxlat) + ' -I15s+e/15s+e > ' + topo_xyz_file
#        #os.system(command_extract)
#        topo_xyz_file = path_topo + '/bedrock.xyz' #in questo modo considera tutta l'area, altrimenti con gmt blockmean vado a selezionare solo una sottoarea
#        topo_utm_file = folder_mesh + '/ptopo.mean.utm'
#        command_cp = 'cp ' + '$HOME/rapids/topo_tools/convert_lonlat2utm.pl' + ' ' + folder_mesh
#        os.system(command_cp)
#        command_convert2utm = folder_mesh + '/convert_lonlat2utm.pl ' + topo_xyz_file + ' ' + str(zone) + ' > ' + \
#                              topo_utm_file
#        os.system(command_convert2utm)
#        result_find_repeat_line = find_repeat_interval(topo_xyz_file)
#        if result_find_repeat_line is not None:
#            interval = result_find_repeat_line
#        else:
#            sys.exit("Error: No repeating x coordinates found in " + topo_xyz_file)
#        python_script_cubit = folder_mesh + '/read_topo.py'
#        create_read_topo_py(python_script_cubit, topo_utm_file, interval, folder_mesh)
#        command_cubit = path_cubit + ' -nographics python3 ' + python_script_cubit
#        os.system(command_cubit)


    cubit_journal = folder_mesh + '/crustal_model.jou'
    fid = open(cubit_journal, 'w')
    fid.write('{}\n'.format('reset # clear the CUBIT database of the current geometry and mesh model'))
    fid.write('{}\n'.format('undo off # hopefully this should improve speed'))
    fid.write('{}\n'.format('graphics off # hopefully this should improve speed'))
    fid.write('{}\n'.format(''))
    if topo == 'yes':
        fid.write('{}\n'.format('# Reading topography CUB file, and cropping out the unwanted portion'))
        #file_cub = folder_mesh + '/topo.cub'
        file_cub = path_topo + '/Friuli.cub' #non lo prende più dal file creato in real time, ma dalla cartella DATA
        #file_cub = path_topo + '/topo.cub' #non lo prende più dal file creato in real time, ma dalla cartella DATA
        # file_stl = folder_mesh + '/topo.stl'
        # command = 'import stl \"'+file_stl+'\" nofreesurfaces heal attributes_on  separate_bodies'
        command = 'open \"' + file_cub + '\"'
        fid.write('{}\n'.format(command))
        # healing the surface...
        fid.write('{}\n'.format('auto_clean volume 1 small_surfaces small_curve_size 10'))
        fid.write('{}\n'.format('regularize volume 1'))
        fid.write('{}\n'.format('compress ids'))
        fid.write('{}\n'.format(''))

    if topo == 'no':
        fid.write('{} {} {} {}\n'.format('create vertex', coord1[0], coord1[1], 0))
        fid.write('{} {} {} {}\n'.format('create vertex', coord2[0], coord2[1], 0))
        fid.write('{} {} {} {}\n'.format('create vertex', coord3[0], coord3[1], 0))
        fid.write('{}\n'.format(''))
        fid.write('{}\n'.format('create surface parallelogram vertex 1 2 3'))
        fid.write('{}\n'.format('volume all scale 0.001'))
    else:
        fid.write('{} {} {} {}\n'.format('create vertex', coord1[0], coord1[1], -(layers['depth_top_layer']+layers['thk'][0]) * 1000))
        fid.write('{} {} {} {}\n'.format('create vertex', coord2[0], coord2[1], -(layers['depth_top_layer']+layers['thk'][0]) * 1000))
        fid.write('{} {} {} {}\n'.format('create vertex', coord3[0], coord3[1], -(layers['depth_top_layer']+layers['thk'][0]) * 1000))
        fid.write('{}\n'.format(''))
        fid.write('{}\n'.format('create surface parallelogram vertex 5 6 7'))
        fid.write('{}\n'.format(''))
        fid.write('{} {} {} {}\n'.format('create vertex', coord1[0], coord1[1], 7000))
        fid.write('{}\n'.format('create curve vertex 5 9'))
        fid.write('{}\n'.format('sweep surface 2 along curve 9 # creation of voulme 2'))
        fid.write('{}\n'.format('volume all scale 0.001'))
        fid.write('{}\n'.format('curve 9 scale 0.001'))
        fid.write('{}\n'.format('webcut volume 2 with sheet body 1'))
        fid.write('{}\n'.format(''))
        fid.write('{}\n'.format('surface 13 copy move x 0 y 0 z 0       # Top Surface'))
        fid.write('{} {} {}\n'.format('surface 7 copy move x 0 y 0 z', 0, '# Bottom Surface'))
        fid.write('{}\n'.format(''))
        fid.write('{}\n'.format('delete curve 9'))
        fid.write('{}\n'.format('delete volume 2 3'))
        fid.write('{}\n'.format('delete volume 1'))
        fid.write('{}\n'.format(''))
        fid.write('{}\n'.format('compress ids'))
        fid.write('{}\n'.format(''))
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('## Making Velocity Layer Interfaces'))
    fid.write('{}\n'.format('## Shifting topography layer by depth of each 1D velocity layer'))
    fid.write('{}\n'.format(''))
    nlayers = len(layers['rho'])

    if topo == 'no':
        fid.write('{} {}\n'.format('surface 1 copy move x 0 y 0 z', -layers['thk'][0]))  # create first layer alone for consistency with topo procedure
        depth_layers0 = layers['thk'][0]
    if topo == 'yes':
        depth_layers0 = layers['depth_top_layer'] + layers['thk'][0]
    depth_layers = depth_layers0
    for i in range(1, nlayers):
        depth_layers += layers['thk'][i] 
        fid.write('{} {}\n'.format('surface 2 copy move x 0 y 0 z', -(depth_layers - depth_layers0)))
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('# create all vertical joining curves'))
    for j in range(4):
        fid.write('{}\n'.format(''))
        vertex_id_1 = j + 1
        for i in range(nlayers):
            vertex_id_2 = vertex_id_1 + 4
            if topo == 'no':
                num1 = vertex_id_1
                num2 = vertex_id_2
            else:
                num1 = rename_vertices(vertex_id_1)
                num2 = rename_vertices(vertex_id_2)
            fid.write('{} {} {}\n'.format('create curve vertex', num1, num2))
            vertex_id_1 = vertex_id_2
    num_max = vertex_id_2
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('# create all vertical joining surfaces'))
    curve_id_3 = num_max
    for j in range(4):
        curve_id_1 = j + 1
        fid.write('{}\n'.format(''))
        for i in range(nlayers):
            curve_id_2 = curve_id_1 + 4
            if j == 3:
                if i == 0:
                    curve_id_4 = num_max + 1
                else:
                    curve_id_4 = curve_id_4 + 1
            else:
                curve_id_4 = curve_id_3 + (nlayers + 1)
            curve_id_3 = curve_id_3 + 1
            if topo == 'no':
                num1 = curve_id_1
                num2 = curve_id_2
                num3 = curve_id_3
                num4 = curve_id_4
            else:
                num1 = rename_curves(curve_id_1)
                num2 = rename_curves(curve_id_2)
                num3 = rename_curves(curve_id_3)
                num4 = rename_curves(curve_id_4)
            fid.write('{} {} {} {} {}\n'.format('create surface curve', num1, num2, num3, num4))
            curve_id_1 = curve_id_2
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('# create all volumes from joining surfaces'))
    fid.write('{}\n'.format('#  heal option will attempt to close small gaps in the surface'))
    fid.write('{}\n'.format('# keep option preserves  original surfaces'))
    fid.write('{}\n'.format(''))
    for i in range(nlayers):
        num1 = i + 1
        num2 = num1 + 1
        num3 = num2 + nlayers
        num4 = num3 + nlayers
        num5 = num4 + nlayers
        num6 = num5 + nlayers
        fid.write('{} {} {} {} {} {} {} {}\n'.format('create volume surface', num1, num2, num3,
                                                     num4, num5, num6, 'heal keep'))
    fid.write('{}\n'.format(''))
    num_max2 = num6
    fid.write('{} {}\n'.format('del body 1 to', num_max2))
    fid.write('{}\n'.format('merge volume all'))
    fid.write('{}\n'.format('compress ids'))
    command_cub_volumes = 'save as "' + folder_mesh + '/Volumes.cub" overwrite'
    fid.write('{}\n'.format(command_cub_volumes))
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('#*************************************'))
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('# Initial Coarse Mesh'))
    mesh_size_from_vel = define_mesh_size(layers, computational_param)
    mesh_size_from_vel = mesh_size_from_vel/1000.
    if topo == 'yes':
        resolution_topo = 0.325  # m corrispondenti a 15 arcsec alle nostre latitudini
        if mesh_size_from_vel/3. > resolution_topo:
            mesh_size = resolution_topo * 3
            need_refinement = 1
        else:
            mesh_size = min(mesh_size_from_vel, resolution_topo)
            need_refinement = 0
    else:
        mesh_size = mesh_size_from_vel
    print('mesh_size= ', mesh_size)
    fid.write('{} {}\n'.format('volume all size', mesh_size))
    if topo == 'yes':
        # note: we will mesh first the topography surface, then sweep down the mesh
        # topography surface
        # cubit.cmd('control skew surface 2')
        fid.write('{}\n'.format('surface 2 submap smooth off'))
        fid.write('{}\n'.format('surface 2 scheme submap'))
        fid.write('{}\n'.format('mesh surface 2'))
        # propagates mesh down for whole volume
        fid.write('{}\n'.format('volume all redistribute nodes off'))
        fid.write('{}\n'.format('volume 1 scheme sweep source surface 2 target surface 3'))
        fid.write('{}\n'.format(''))
    num1 = 1
    for i in range(nlayers):
        num2 = math.ceil(layers['thk'][i] / mesh_size_from_vel)
        if topo == 'yes':
            if i > 0:
                fid.write('{} {} {} {}\n'.format('curve', num1, 'interval', num2))
        else:
            fid.write('{} {} {} {}\n'.format('curve', num1, 'interval', num2))
        if i == 0:
            num1 = num1 + 12
        else:
            num1 = num1 + 8
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('volume all scheme auto'))
    fid.write('{}\n'.format('mesh volume all'))
    fid.write('{}\n'.format(''))

    command_cub_mesh_before_refinement = 'save as "' + folder_mesh + '/Mesh_before_refinement.cub" overwrite'
    fid.write('{}\n'.format(command_cub_mesh_before_refinement))
    if topo == 'yes':
        if need_refinement == 1:
            refine_command = 'refine surface 2 numsplit 1 bias 1.0 depth 2' #suddivide il lato in 3 con numsplit 1
        fid.write('{}\n'.format(refine_command))
        fid.write('{}\n'.format('#*************************************'))
        fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('# Exporting'))
    for i in range(nlayers):
        num1 = i + 1
        fid.write('{} {} {} {}\n'.format('block', num1, 'volume', num1))
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('# Absorbing Boundaries at surfaces'))
    abso_block = nlayers + 1
    abso_list = []
    s1 = 1
    s2 = 4
    s3 = 5
    s4 = 6
    s5 = 3
    abso_list.append(s1)
    abso_list.append(s2)
    abso_list.append(s3)
    abso_list.append(s4)
    if nlayers == 1:
        abso_list.append(3)
    else:
        for i in range(1, nlayers):
            if i == 1:
                s1 = s1 + 6
            else:
                s1 = s1 + 5
            s2 = s2 + 5
            s3 = s3 + 5
            s4 = s4 + 5
            s5 = s5 + 5
            abso_list.append(s1)
            abso_list.append(s2)
            abso_list.append(s3)
            abso_list.append(s4)
            if i == nlayers - 1:
                abso_list.append(s5)

    abso_string = listToString(abso_list)
    fid.write('{} {} {} {}\n'.format('block', abso_block, 'surface', abso_string))
    fid.write('{} {}\n'.format('draw block', abso_block))
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('list model'))
    fid.write('{}\n'.format(''))
    command_cub_mesh = 'save as "' + folder_mesh + '/Mesh.cub" overwrite'
    fid.write('{}\n'.format(command_cub_mesh))
    fid.write('{}\n'.format('# Exporting Meshed Volume to Exodus Format'))
    fid.write('{}\n'.format('set exodus netcdf4 off'))
    fid.write('{}\n'.format('set large exodus file off'))
    file_exodus = folder_mesh + '/crustal_model.e'
    fid.write('{}\n'.format('transform mesh output scale 1000'))
    fid.write('{}{}{}\n'.format('export mesh "', file_exodus, '" dimension 3 overwrite'))

    if topo == 'yes':
        fid.write('{}\n'.format(''))
        fid.write('{}\n'.format('# Start topography'))
        fid.write('{}\n'.format('reset'))
        file_volumes = folder_mesh + '/Volumes.cub'
        command = 'open \"' + file_volumes + '\"'
        fid.write('{}\n'.format(command))
        fid.write('{}\n'.format('surface 2 copy move x 0 y 0 z 0 '))
        surf_topo = num_max2 + 1
        resolution_topo = math.ceil(mesh_size)
        fid.write('{} {} {} {}\n'.format('Surface', surf_topo, 'scheme trimesh size', resolution_topo))
        fid.write('{} {}\n'.format('mesh surf', surf_topo))
        block_topo = abso_block + 1
        fid.write('{} {} {} {}\n'.format('block', block_topo, 'surf', surf_topo))
        file_exodus_topo = folder_mesh + '/Topography.e'
        fid.write('{}\n'.format('transform mesh output scale 1000'))
        fid.write('{}{}{} {} {}\n'.format('export mesh "', file_exodus_topo, '" dimension 3 block', block_topo,
                                          'overwrite'))
    else:
        file_exodus_topo = ''

    fid.write('{}\n'.format('exit'))
    fid.close()
    return cubit_journal, file_exodus, file_exodus_topo


def create_input_mate(folder, computational_param, layers, fault):
    import numpy as np
    fid = open(folder + '/fault_layer.mate', 'w')
    for i in range(len(layers['thk'])):
        blockid = i + 1

        fid.write('{} {} {} {} {} {} {} {}\n'.format('MATE', blockid, computational_param['spectral_degree'],
                                                     layers['rho'][i] * 1000, layers['vs'][i] * 1000,
                                                     layers['vp'][i] * 1000,
                                                     layers['qs'][i], layers['qp'][i]))
    fid.write('\n')
    blockid_abso = blockid + 1
    fid.write('{} {}\n'.format('ABSO', blockid_abso))
    fid.write('\n')
    if computational_param['damping'] == 1:
        fid.write('{} {}\n'.format('FMAX', computational_param['fval_quality_factor']))
        fid.write('\n')
    IDfunc = blockid_abso + 1
    tau = fault['rise_time']
    if fault['IDx'] == 'exp':
        tau = fault['rise_time'] / 0.25 / (2 * np.pi)
        IDx = 14
        amplitude = 1.00
        fid.write('{} {} {} {} {}\n'.format('FUNC', IDfunc, IDx, tau, amplitude))
    fid.write('\n')
    if fault['IDx'] == 'Archuleta' or fault['IDx'] == 'Yoffe-DCF':
        IDx = 33
        num_timevalues, fileSTF = compute_STF(folder, fault)
        fid.write('{} {} {} {} {}\n'.format('FUNC', IDfunc, IDx, num_timevalues, fileSTF))
    fid.write('\n')
    if fault['fault_type'] == 'point':
        slip_vec = compute_slip_vec(fault['strike'], fault['rake'], fault['dip'])
        norm_vec = compute_norm_vec(fault['strike'], fault['rake'], fault['dip'])
        Trup = 0
        fid.write(f"SISM {IDfunc} 0 {fault['hypo_utm']['Y']} {fault['hypo_utm']['X']} {-fault['hypo_utm']['Z']} "
                  f"{fault['hypo_utm']['Y']} {fault['hypo_utm']['X']} {-fault['hypo_utm']['Z']} "
                  f"{fault['hypo_utm']['Y']} {fault['hypo_utm']['X']} {-fault['hypo_utm']['Z']} "
                  f"{fault['hypo_utm']['Y']} {fault['hypo_utm']['X']} {-fault['hypo_utm']['Z']} "
                  f"{slip_vec['X']} {slip_vec['Y']} {slip_vec['Z']} "
                  f"{norm_vec['X']} {norm_vec['Y']} {norm_vec['Z']} "
                  f"{Trup} {fault['Mo']} {tau}\n")

    if fault['fault_type'] == 'extended':
        fid.write('{}\n'.format('SLIP LOAD-SRCMOD2'))
        segdata, vertdata, nseg, segStk, segDip, HypGlo = \
            define_param_sism_commands(computational_param['ndiv'], folder, fault)

        for iseg in range(nseg):
            x = segdata[iseg]['x']
            y = segdata[iseg]['y']
            z = segdata[iseg]['z']
            slip = segdata[iseg]['slip']
            if fault['IDx'] == 'exp':
                tau = segdata[iseg]['tau'] / 0.25 / (2 * np.pi)
            else:
                tau = segdata[iseg]['tau']
            trup = segdata[iseg]['trup']
            rake = segdata[iseg]['rake']
            SeisMom = segdata[iseg]['Mo']

            xvert = vertdata[iseg]['xvert']
            yvert = vertdata[iseg]['yvert']
            zvert = vertdata[iseg]['zvert']

            npts = len(x)

            for ipt in range(npts):
                if slip[ipt] > 0:
                    SV = compute_slip_vec(segStk[iseg], rake[ipt], segDip[iseg])
                    NV = compute_norm_vec(segStk[iseg], rake[ipt], segDip[iseg])

                    # SISM IDfunc dummy [x y z]hypo [x y z]centroid [x y z]slip [x y z]norm Trup Mo Trise
                    fid.write(f"SISM  {IDfunc}  2  {HypGlo[0]:+13.7e}  {HypGlo[1]:+13.7e}  {HypGlo[2]:+13.7e}      "
                              f"{x[ipt]:+13.7e}  {y[ipt]:+13.7e}  {z[ipt]:+13.7e}      {SV['X']:+13.7e}  "
                              f"{SV['Y']:+13.7e}  {SV['Z']:+13.7e}      {NV['X']:+13.7e}  {NV['Y']:+13.7e}  "
                              f"{NV['Z']:+13.7e}     {trup[ipt]:+13.7e}    {SeisMom[ipt]:+13.7e}    "
                              f"{tau[ipt]:+13.7e} \n")

    fid.close()
    return


def create_input_LS(folder, sites, topo, path_data, fault, computational_param):
    from rapids.conversions import determine_utm_coord, utm_to_lon_lat
    import os
    import numpy as np
    nobs = len(sites['Z'])
    if topo == 'yes':
        path_topo = path_data+'/topo'
        coord1, coord2, coord3, coord4 = define_area(fault, sites, computational_param)
        temp1, temp2, zone, letter = determine_utm_coord(fault['hypo']['lon'], fault['hypo']['lat'])
        coord1_lon, coord1_lat = utm_to_lon_lat(coord1[0], coord1[1], zone)
        coord2_lon, coord2_lat = utm_to_lon_lat(coord2[0], coord2[1], zone)
        coord3_lon, coord3_lat = utm_to_lon_lat(coord3[0], coord3[1], zone)
        coord4_lon, coord4_lat = utm_to_lon_lat(coord4[0], coord4[1], zone)
        buffer = 0.2
        minlon = min(coord1_lon, coord2_lon, coord3_lon, coord4_lon) - buffer
        maxlon = max(coord1_lon, coord2_lon, coord3_lon, coord4_lon) + buffer
        minlat = min(coord1_lat, coord2_lat, coord3_lat, coord4_lat) - buffer
        maxlat = max(coord1_lat, coord2_lat, coord3_lat, coord4_lat) + buffer
        topo_xyz_file = folder + '/topo_reduced.xyz'
        command_extract = 'gmt blockmean ' + path_topo + '/bedrock.xyz -R' + str(minlon) + '/' + str(
                maxlon) + '/' + str(minlat) \
                + '/' + str(maxlat) + ' -I15s+e/15s+e > ' + topo_xyz_file
        os.system(command_extract)
        elevation = compute_elevation(sites, topo_xyz_file)
    else:
        elevation = np.zeros((nobs))
        for k in range(nobs):
            elevation[k] = sites['Z'][k]

    fid = open(folder + '/LS.input', 'w')
    fid.write('{}\n'.format(nobs))
    for k in range(nobs):
        fid.write('{} {} {} {}\n'.format(k + 1, sites['Y'][k], sites['X'][k], elevation[k]))
    fid.close()

    return


def compute_elevation(sites, file_topo):
    import pandas as pd
    import numpy as np
    from scipy.spatial import cKDTree

    topo_val = pd.read_csv(file_topo, sep='\t', names = ['lon','lat','elev'])

    sites_coords = np.vstack((sites['lat'], sites['lon'])).T
    topo_coords = np.vstack((topo_val['lat'], topo_val['lon'])).T
    
    # Build a KDTree for efficient nearest neighbor search
    tree = cKDTree(topo_coords)
    
    # Query the KDTree to find the index of the nearest topo point for each site
    _, nearest_idx = tree.query(sites_coords, k=1)
    
    # Get the elevation corresponding to the nearest topo point for each site
    elevation = topo_val['elev'].values[nearest_idx]

    return elevation


def haversine_np(lon1, lat1, lon2, lat2):
    import numpy as np
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    Reference:
        https://stackoverflow.com/a/29546836/7657658
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(
        dlat / 2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6371 * c
    return km

def create_input_speed_file(folder, computational_param, meshfile):
    fid = open(folder + '/SPEED.input', 'w')
    fid.write('{}  {}\n'.format('GRIDFILE', meshfile))
    fid.write('{}   {}\n'.format('MATFILE', 'fault_layer'))
    fid.write('{}   {}\n'.format('MPIFILE', 'FILES_MPI'))
    fid.write('{}   {}\n'.format('MONFILE', 'MONITOR'))
    fid.write('{}\n'.format(''))
    fid.write('{}  {}\n'.format('DAMPING', computational_param['damping']))
    fid.write('{}\n'.format(''))
    fid.write('{} {} {} {} {} {} {}\n'.format('OPTIOUT', int(computational_param['optiout'][0]),
                                              int(computational_param['optiout'][1]),
                                              int(computational_param['optiout'][2]),
                                              int(computational_param['optiout'][3]),
                                              int(computational_param['optiout'][4]),
                                              int(computational_param['optiout'][5])))
    fid.write('{}\n'.format(''))
    fid.write('{}  {}\n'.format('TIMESTEP', computational_param['dt_speed']))
    fid.write('{}  {}\n'.format('TMONITOR', computational_param['tmonitor']))
    fid.write('{}  {}\n'.format('STOPTIME', computational_param['stoptime']))
    fid.write('{}\n'.format(''))
    fid.write('{}   {} {}\n'.format('MLST', int(computational_param['mlst'][0]), int(computational_param['mlst'][1])))
    fid.write('{}  {}\n'.format('SETUPONL', computational_param['setuponl']))
    fid.close()
    return


def create_script_speed(folder, path_code_speed, sites):
    import numpy as np

    fid = open(folder + '/run_speed_local.sh', 'w')
    fid.write('{}\n'.format('#!/bin/bash'))
    fid.write('{}\n'.format(''))
    string = 'export SPEED_PROGRAM=' + path_code_speed + '"'
    fid.write('{}\n'.format(string))
    fid.write('{}\n'.format('export OUTPUT=speed.out'))
    fid.write('{}\n'.format('export TOTAL_MPI_PROCESSES=1'))
    fid.write('{}\n'.format('export THREADS=1'))
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('let NCORES=$TOTAL_MPI_PROCESSES'))
    fid.write('{}\n'.format('echo "Running on $NCORES cores"'))
    fid.write('{}\n'.format('echo "Job started at `date`"'))
    fid.write('{}\n'.format(''))
    fid.write('{}\n'.format('mkdir -p MONITOR'))
    fid.write('{}\n'.format('mkdir -p FILES_MPI'))
    fid.write('{}\n'.format('rm MONITOR/*'))
    fid.write('{}\n'.format('rm FILES_MPI/*'))
    fid.write('{}\n'.format(''))
    fid.write(
        '{}\n'.format('time mpirun -np $TOTAL_MPI_PROCESSES  -x OMP_NUM_THREADS=$THREADS  $SPEED_PROGRAM >& $OUTPUT'))
    fid.write('{}\n'.format('wait'))
    fid.close()
