def read_slip_xta(folder, fault, layers):
    from rapids.plot import plot_slip_from_file
    import numpy as np

    values = np.loadtxt(fault['file_xta'])
    nx = fault['nsubcells_length_xta']
    ny = fault['nsubcells_width_xta']
    fault['number_subfaults_strike'] = nx
    fault['number_subfaults_dip'] = ny

    moment_matrix = values[:, 3].reshape(ny, nx)
    moment_matrix = moment_matrix * 10 ** 20  # dyne*cm
    moment_matrix = moment_matrix * 10 ** (-7)  # Nm
    rupture_matrix = values[:, 2].reshape(ny, nx)

    fault['subfault_length'] = fault['length']/nx
    fault['subfault_width'] = fault['width']/ny

    moment_xta = []
    z_xta = []
    rake_xta = []
    rupture_time_xta = []
    rise_time_xta = []
    for j in range(ny):
        for k in range(nx):
            moment_xta.append(moment_matrix[j, k])
            depth_subsource = j * fault['subfault_width'] + fault['subfault_width'] / 2
            z_xta.append(depth_subsource * 1000)  # in m
            #dist_x_subsource = k * fault['subfault_length'] + fault['subfault_length'] / 2  # in km
            rake_xta.append(fault['rake'])
            rise_time_xta.append(fault['rise_time'])
            #dist_subsource = np.sqrt((fault['hypo_down_dip'] * fault['width'] - depth_subsource) ** 2 + (
            #        fault['hypo_along_strike'] * fault['length'] - dist_x_subsource) ** 2)
            #trup_subsource = dist_subsource / fault['rupture_velocity']
            rupture_time_xta.append(rupture_matrix[j, k])

    moment_xta = np.asarray(moment_xta)
    z_xta = np.asarray(z_xta)
    rake_xta = np.asarray(rake_xta)
    rise_time_xta = np.asarray(rise_time_xta)
    rupture_time_xta = np.asarray(rupture_time_xta)
    slip_xta = convert_ucsb_moment2slip(moment_xta, z_xta, layers, fault)
    plot_slip_from_file(folder, fault, slip_xta, rupture_time_xta)
    return slip_xta, rake_xta, rise_time_xta, rupture_time_xta


def convert_ucsb_moment2slip(moment, z_glo, layers, fault):
    import numpy as np
    b = layers['vs'] * 1000  # in m/s
    rho = layers['rho'] * 1000
    h = np.cumsum(layers['thk']) * 1000

    A = fault['subfault_length'] * fault['subfault_width'] * 1e6
    # Material Property at each point
    lay_num = np.zeros((len(moment)))
    for ipt in range(len(lay_num)):
        indx1 = np.where(h >= z_glo[ipt])
        if not indx1:
            lay_num[ipt] = len(h) - 1
        else:
            lay_num[ipt] = indx1[0][0]
    lay_num = np.asarray(lay_num, dtype=int)

    slip = np.zeros((len(lay_num)))
    mu = np.zeros((len(lay_num)))
    for ipt in range(len(lay_num)):
        mu[ipt] = rho[lay_num[ipt]] * b[lay_num[ipt]] * b[lay_num[ipt]]
        slip[ipt] = moment[ipt] / (mu[ipt] * A)
    return slip


def compute_derivative(u, dt):
    import numpy as np
    dudt = np.copy(u)
    npts = len(u)
    for j in range(1, npts - 1):
        dudt[j] = (u[j + 1] - u[j - 1]) / (2 * dt)
    dudt[0] = 0.0
    dudt[npts - 1] = dudt[npts - 2]
    return dudt


def find_rotated_coords(x_local, y_local, Theta):
    import numpy as np
    # Used to find the Global coordinates of (x_local,y_local).
    # where (x_local,y_local) are with respect to a local coordinate system with
    # Origin as (0,0) and local_x_axis makes an angle "Theta" with respect to
    # global_x_axis in anticlockwise direction

    if isinstance(x_local, float):
        x_local = np.array([x_local])
    if isinstance(y_local, float):
        y_local = np.array([y_local])

    x_local_orig = x_local
    y_local_orig = y_local

    x_local = x_local.flatten()
    y_local = y_local.flatten()

    xglobal = np.zeros((np.size(x_local)))
    yglobal = np.zeros((np.size(y_local)))

    for i in range(np.size(x_local)):
        Radi = np.sqrt((x_local[i] ** 2) + (y_local[i] ** 2))
        np.seterr(divide='ignore', invalid='ignore')  # suppress warning in divide by 0
        rto = y_local[i] / x_local[i]

        if np.sign(rto) == 1 and x_local[i] > 0:  # 1st Quadrant in Local coordinates
            phi = np.rad2deg(np.arctan(rto))
        elif np.sign(rto) == -1 and x_local[i] < 0:  # 2nd Quadrant in Local coordinates
            phi = 180 + np.rad2deg(np.arctan(rto))
        elif np.sign(rto) == 1 and x_local[i] < 0:  # 3rd Quadrant in Local coordinates
            phi = 180 + np.rad2deg(np.arctan(rto))
        elif np.sign(rto) == -1 and x_local[i] > 0:  # 4th Quadrant in Local coordinates
            phi = np.rad2deg(np.arctan(rto))
        elif np.sign(rto) == 0 and x_local[i] > 0:
            phi = 0
        elif np.sign(rto) == 0 and x_local[i] < 0:
            phi = 180
        elif np.sign(rto) == -1 and x_local[i] == 0:
            phi = 270
        elif np.sign(rto) == 1 and x_local[i] == 0:
            phi = 90
        elif Radi == 0:
            phi = 0

        xglobal[i] = Radi * np.cos(np.deg2rad(phi + Theta))
        yglobal[i] = Radi * np.sin(np.deg2rad(phi + Theta))

    x_global = np.reshape(xglobal, x_local_orig.shape)
    y_global = np.reshape(yglobal, y_local_orig.shape)

    return x_global, y_global


def create_vertex_utm(pbl, pbr, ptl, ptr):
    vertex_utm = {
        "pbl": create_point_utm(pbl['X'], pbl['Y'], pbl['Z']),
        "pbr": create_point_utm(pbr['X'], pbr['Y'], pbr['Z']),
        "ptl": create_point_utm(ptl['X'], ptl['Y'], ptl['Z']),
        "ptr": create_point_utm(ptr['X'], ptr['Y'], ptr['Z']),
    }
    return vertex_utm


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
        print(startend[iseg, 1])
        nlines = int(startend[iseg, 1] - startend[iseg, 0] - 1)
        print(nlines)
        srcmoddata.append(np.loadtxt(file_name, max_rows=nlines, skiprows=int(startend[iseg, 0] + 1)))
        print(srcmoddata)

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


def create_point(point_lon, point_lat, point_z):
    point = {
        "lon": point_lon,
        "lat": point_lat,
        "Z": point_z
    }
    return point


def create_vertex(pbl, pbr, ptl, ptr):
    vertex = {
        "pbl": create_point(pbl.x, pbl.y, pbl.z),
        "pbr": create_point(pbr.x, pbr.y, pbr.z),
        "ptl": create_point(ptl.x, ptl.y, ptl.z),
        "ptr": create_point(ptr.x, ptr.y, ptr.z),
    }
    return vertex


def determine_fault_coordinates_from_hypocentre(fault):
    import numpy as np
    from shapely.geometry import Point

    strike = fault['strike']
    dip = fault['dip']
    length = fault['length']
    width = fault['width']
    hypo_along_strike = fault['hypo_along_strike']
    hypo_down_dip = fault['hypo_down_dip']
    hypo_depth = fault['hypo']['Z']
    lon = fault['hypo']['lon']
    lat = fault['hypo']['lat']

    #
    #     .....     the dotted line is the hdist
    #     \      |
    #      \     |  this dashed vertical line is the height
    #       \    |
    #        \   |
    # rupture \  |
    #

    # https://docs.openquake.org/oq-engine/latest/reference/_modules/openquake/hazardlib/geo/surface/planar.html

    # if ztor is not None:
    #	depth = ztor + height / 2

    height = width * np.sin(np.radians(dip))
    hdist = width * np.cos(np.radians(dip))

    # Move hor. hypo_down_dip*hdist in direction -90
    hypo_top = point_at(lon, lat, strike - 90, hypo_down_dip * hdist)  # dist in km
    # Move hor. (1-hypo_down_dip)* hdist in direction +90
    hypo_bot = point_at(lon, lat, strike + 90, (1 - hypo_down_dip) * hdist)

    if strike <= 180:
        # compute corner points at the surface
        top_right = point_at(hypo_top[0], hypo_top[1], strike, (1 - hypo_along_strike) * length)
        top_left = point_at(hypo_top[0], hypo_top[1], strike + 180, hypo_along_strike * length)
        bot_right = point_at(hypo_bot[0], hypo_bot[1], strike, (1 - hypo_along_strike) * length)
        bot_left = point_at(hypo_bot[0], hypo_bot[1], strike + 180, hypo_along_strike * length)
    else:
        top_left = point_at(hypo_top[0], hypo_top[1], strike + 180, hypo_along_strike * length)
        top_right = point_at(hypo_top[0], hypo_top[1], strike, (1 - hypo_along_strike) * length)
        bot_left = point_at(hypo_bot[0], hypo_bot[1], strike + 180, hypo_along_strike * length)
        bot_right = point_at(hypo_bot[0], hypo_bot[1], strike, (1 - hypo_along_strike) * length)

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
    if lat > 0:
        hemisphere = 'north'
    else:
        hemisphere = 'south'
    utm_zone_num = math.ceil((lon + 180) / 6)
    utm_zone_letter = chr(ord('C') + int((lat + 80) / 8))
    if utm_zone_letter > 'H':
        zone_letter = chr(ord(utm_zone_letter) + 1)
    if zone_letter > 'N':
        utm_zone_letter = chr(ord(zone_letter) + 1)
    myProj = Proj("+proj=utm +zone=" + str(utm_zone_num) + " +" + hemisphere + " +datum=WGS84 +units=m +no_defs ")
    UTMx, UTMy = myProj(lon, lat, inverse=False)
    # x Easting
    # y Northing
    return UTMx, UTMy, utm_zone_num, utm_zone_letter


def utm_to_lon_lat(easting, northing, zone_number):
    import pyproj
    utm_coordinate_system = pyproj.CRS.from_string(f"+proj=utm +zone={zone_number} +ellps=WGS84")
    lon_lat_coordinate_system = pyproj.CRS.from_string("+proj=latlong +datum=WGS84")
    transformer = pyproj.Transformer.from_crs(utm_coordinate_system, lon_lat_coordinate_system, always_xy=True)
    lon, lat = transformer.transform(easting, northing)
    return lon, lat


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


def create_point_utm(point_x, point_y, point_z):
    point = {
        "X": point_x,
        "Y": point_y,
        "Z": point_z
    }
    return point


def create_receivers(receivers_X, receivers_Y, receivers_Z, receivers_ID):
    sites = {
        "X": receivers_X,
        "Y": receivers_Y,
        "Z": receivers_Z,
        "ID": receivers_ID
    }
    return sites


def convert_mesh(fileName, type_conv, folder):
    import numpy as np
    import os
    import re

    print('\n\t Starting make_MeshFile ...\n')

    fileNameRead = fileName
    (FilePath, Ext) = os.path.splitext(fileName)

    InfoMesh = {
        'FilePath': FilePath,
        'FileExt': Ext
    }

    if type_conv == 'mesh':
        head, tail = os.path.split(FilePath)
        fileNameWrite = folder + '/' + tail + '.mesh'

    if type_conv == 'xyz':
        fileNameWrite = folder+'/XYZ.out'

    chexa = 0
    ctetra = 0
    cpyram = 0
    cquad = 0
    ctria = 0

    num_chexa = 0
    num_cpyram = 0
    num_ctetra = 0

    num_cquad = 0
    num_ctria = 0

    with open(fileName, 'r') as fid:
        # read first 6 rows
        for _ in range(5):
            tline = fid.readline()
        # read nodes
        tline = fid.readline()
        num_nodes = int(tline.split()[2])
        # read elements
        tline = fid.readline()
        num_elem = int(tline.split()[2])
        # read block
        tline = fid.readline()
        num_el_blk = int(tline.split()[2])
        #skip one line
        tline = fid.readline()

        num_el_in_blk = []
        for i in range(num_el_blk):
            string = tline.split()[0][0:7]
            if string == 'num_el_':
                num_el_in_blk.append(int(tline.split()[2]))
                tline = fid.readline()
                string = tline.split()[0][0:7]
            if string == 'num_nod':
                tline = fid.readline()
                string = tline.split()[0][0:7]
            if string == 'num_att':
                tline = fid.readline()
                string = tline.split()[0][0:7]

        while (tline[2:9] != 'connect'):
            tline = fid.readline()
        # tipo_elementi
        # 1 HEXA
        # 2 TETRA
        # 3 PYRA
        # 4 QUAD
        # 5 TRI
        ##########################################
        # Number of elements in a block
        type_elem = []
        for i in range(num_el_blk):
            string = tline.split()[2][1:len(tline.split()[2])-1]
            if string == 'HEX8':
                type_elem.append(1)
                num_chexa += num_el_in_blk[i]
            elif string == 'TETRA':
                type_elem.append(2)
                num_ctetra += num_el_in_blk[i]
            elif string == 'PYRAMID5':
                type_elem.append(3)
                num_cpyram += num_el_in_blk[i]
            elif string == 'TRI3':
                type_elem.append(5)
                num_ctria += num_el_in_blk[i]
            else:
                type_elem.append(4)
                num_cquad += num_el_in_blk[i]
            tline = fid.readline()
            string = tline[8:14]
            if string == 'attrib':
                tline = fid.readline()
            tline = fid.readline()

        InfoMesh["NumHexa"] = num_chexa
        InfoMesh["NumQuad"] = num_cquad

        # Node Coordinates
        grid_id = [0] * num_nodes

        #print('\n\t BEGIN - Reading nodes coordinates ... ')
        grid = num_nodes
        grid_id = list(range(1,num_nodes+1))

        # Opening the whole file as a character array
        dum_ind = []
        with open(fileNameRead, 'r') as file:
            mesh_dat = file.read()
        for i in list(range(len(mesh_dat))):
            if mesh_dat[i] == ';':
                dum_ind.append(i)

        # X-coordinates
        ind_node_start = mesh_dat.find(' coordx =') + 9
        srch_ind = [i for i, val in enumerate(dum_ind) if val >= ind_node_start]
        ind_node_end = dum_ind[srch_ind[0]] - 1
        dum = re.sub(r'[\n\r]+','',mesh_dat[ind_node_start:ind_node_end].strip())
        grid_x = list(map(float, dum.split(',')))

        # Y-coordinates
        ind_node_start = mesh_dat.find(' coordy =') + 9
        srch_ind = [i for i, val in enumerate(dum_ind) if val >= ind_node_start]
        ind_node_end = dum_ind[srch_ind[0]] - 1
        dum = re.sub(r'[\n\r]+','',mesh_dat[ind_node_start:ind_node_end].strip())
        grid_y = list(map(float, dum.split(',')))

        # Z-coordinates
        ind_node_start = mesh_dat.find(' coordz =') + 9
        srch_ind = [i for i, val in enumerate(dum_ind) if val >= ind_node_start]
        ind_node_end = dum_ind[srch_ind[0]] - 1
        dum = re.sub(r'[\n\r]+','',mesh_dat[ind_node_start:ind_node_end].strip())
        grid_z = list(map(float, dum.split(',')))

        #print(f"\n\t END - Reading nodes coordinates in {time.process_time() - tini} sec.")
        l_tline = 8
        tline = fid.readline()

        # Connectivity volume matrices
        con_chexa = np.zeros((num_chexa, 8))
        con_ctetra = np.zeros((num_ctetra, 4))
        con_cpyram = np.zeros((num_cpyram, 5))

        # Connectivity surface matrices
        con_cquad = np.zeros((num_cquad, 4))
        con_ctria = np.zeros((num_ctria, 3))

        # Matrices with tag and id
        chexa_tag = np.zeros((num_chexa, 1))
        chexa_id = np.zeros((num_chexa, 1))

        ctetra_tag = np.zeros((num_ctetra, 1))
        ctetra_id = np.zeros((num_ctetra, 1))

        cpyram_tag = np.zeros((num_ctetra, 1))
        cpyram_id = np.zeros((num_ctetra, 1))

        cquad_tag = np.zeros((num_cquad, 1))
        cquad_id = np.zeros((num_cquad, 1))

        ctria_tag = np.zeros((num_ctria, 1))
        ctria_id = np.zeros((num_ctria, 1))

        chexa_1 = np.zeros((num_chexa, 1))
        chexa_2 = np.zeros((num_chexa, 1))
        chexa_3 = np.zeros((num_chexa, 1))
        chexa_4 = np.zeros((num_chexa, 1))
        chexa_5 = np.zeros((num_chexa, 1))
        chexa_6 = np.zeros((num_chexa, 1))
        chexa_7 = np.zeros((num_chexa, 1))
        chexa_8 = np.zeros((num_chexa, 1))

        ctetra_1 = np.zeros((num_ctetra, 1))
        ctetra_2 = np.zeros((num_ctetra, 1))
        ctetra_3 = np.zeros((num_ctetra, 1))
        ctetra_4 = np.zeros((num_ctetra, 1))

        cpyram_1 = np.zeros((num_cpyram, 1))
        cpyram_2 = np.zeros((num_cpyram, 1))
        cpyram_3 = np.zeros((num_cpyram, 1))
        cpyram_4 = np.zeros((num_cpyram, 1))
        cpyram_5 = np.zeros((num_cpyram, 1))

        cquad_1 = np.zeros((num_cquad, 1))
        cquad_2 = np.zeros((num_cquad, 1))
        cquad_3 = np.zeros((num_cquad, 1))
        cquad_4 = np.zeros((num_cquad, 1))

        ctria_1 = np.zeros((num_ctria, 1))
        ctria_2 = np.zeros((num_ctria, 1))
        ctria_3 = np.zeros((num_ctria, 1))

        for i in range(num_el_blk):
            if i == 0:
                tline = fid.readline()
            while tline[0:l_tline] != ' connect':
                tline = fid.readline()
                if len(tline) < 8:
                    l_tline = len(tline)
                elif not tline:
                    l_tline = 1
                else:
                    l_tline = 8

            if type_elem[i] == 1:
                for j in range(num_el_in_blk[i]):
                    tline = fid.readline()
                    chexa_id[chexa] = chexa + 1
                    chexa_tag[chexa] = i + 1
                    con_chexa[chexa, :] = list(map(int,tline[:-2].split(',')))
                    chexa += 1
            elif type_elem[i] == 2:
                for j in range(num_el_in_blk[i]):
                    tline = fid.readline()
                    ctetra_id[ctetra] = ctetra + 1
                    ctetra_tag[ctetra] = i + 1
                    con_ctetra[ctetra, :] = list(map(int,tline[:-2].split(',')))
                    ctetra += 1
            elif type_elem[i] == 3:
                for j in range(num_el_in_blk[i]):
                    tline = fid.readline()
                    cpyram_id[cpyram] = cpyram + 1
                    cpyram_tag[cpyram] = i + 1
                    con_cpyram[cpyram, :] = list(map(int,tline[:-2].split(',')))
                    cpyram += 1
            elif type_elem[i] == 4:
                for j in range(num_el_in_blk[i]):
                    tline = fid.readline()
                    cquad_id[cquad] = cquad + 1
                    cquad_tag[cquad] = i + 1
                    con_cquad[cquad, :] = list(map(int,tline[:-2].split(',')))
                    cquad += 1
            else:
                for j in range(num_el_in_blk[i]):
                    tline = fid.readline()
                    ctria_id[ctria] = ctria + 1
                    ctria_tag[ctria] = i + 1
                    con_ctria[ctria, :] = list(map(int,tline[:-2].split(',')))
                    ctria += 1

        #print("\n\t Storing data informations")

        if chexa > 0:
            chexa_1[0:chexa] = np.reshape(con_chexa[0:chexa, 0],(chexa,1))
            chexa_2[0:chexa] = np.reshape(con_chexa[0:chexa, 1],(chexa,1))
            chexa_3[0:chexa] = np.reshape(con_chexa[0:chexa, 2],(chexa,1))
            chexa_4[0:chexa] = np.reshape(con_chexa[0:chexa, 3],(chexa,1))
            chexa_5[0:chexa] = np.reshape(con_chexa[0:chexa, 4],(chexa,1))
            chexa_6[0:chexa] = np.reshape(con_chexa[0:chexa, 5],(chexa,1))
            chexa_7[0:chexa] = np.reshape(con_chexa[0:chexa, 6],(chexa,1))
            chexa_8[0:chexa] = np.reshape(con_chexa[0:chexa, 7],(chexa,1))

        if ctetra > 0:
            ctetra_1[0:ctetra] = np.reshape(con_ctetra[0:ctetra, 0],(ctetra,1))
            ctetra_2[0:ctetra] = np.reshape(con_ctetra[0:ctetra, 1],(ctetra,1))
            ctetra_3[0:ctetra] = np.reshape(con_ctetra[0:ctetra, 2],(ctetra,1))
            ctetra_4[0:ctetra] = np.reshape(con_ctetra[0:ctetra, 3],(ctetra,1))

        if cpyram > 0:
            cpyram_1[0:cpyram] = np.reshape(con_cpyram[0:cpyram, 0],(cpyram,1))
            cpyram_2[0:cpyram] = np.reshape(con_cpyram[0:cpyram, 1],(cpyram,1))
            cpyram_3[0:cpyram] = np.reshape(con_cpyram[0:cpyram, 2],(cpyram,1))
            cpyram_4[0:cpyram] = np.reshape(con_cpyram[0:cpyram, 3],(cpyram,1))
            cpyram_5[0:cpyram] = np.reshape(con_cpyram[0:cpyram, 4],(cpyram,1))

        if cquad > 0:
            cquad_1[0:cquad] = np.reshape(con_cquad[0:cquad, 0],(cquad,1))
            cquad_2[0:cquad] = np.reshape(con_cquad[0:cquad, 1],(cquad,1))
            cquad_3[0:cquad] = np.reshape(con_cquad[0:cquad, 2],(cquad,1))
            cquad_4[0:cquad] = np.reshape(con_cquad[0:cquad, 3],(cquad,1))

        if ctria > 0:
            ctria_1[0:ctria] = np.reshape(con_ctria[0:ctria, 0],(ctria,1))
            ctria_2[0:ctria] = np.reshape(con_ctria[0:ctria, 1],(ctria,1))
            ctria_3[0:ctria] = np.reshape(con_ctria[0:ctria, 2],(ctria,1))

    #print('\n\t BEGIN - Writing inp format')
    fid = open(fileNameWrite, 'w')

    if type_conv == 'mesh':
        fid.write('   %i   %i   %i   %i   %i\n' % (grid, chexa + ctetra + cpyram + cquad + ctria, 0, 0, 0))
    else:
        fid.write('    %i      %i\n' % (num_nodes, ctria))

    if grid > 0:
        for i in range(grid):
            fid.write('%i  %+13.7e  %+13.7e  %+13.7e\n' % (grid_id[i], grid_x[i], grid_y[i], grid_z[i]))

    if type_conv == 'mesh':

        if ctria > 0:
            for i in range(ctria):
                fid.write('%i %i tria %i %i %i\n' % (ctria_id[i], ctria_tag[i], ctria_1[i], ctria_2[i], ctria_3[i]))

        if cquad > 0:
            for i in range(cquad):
                fid.write('%i  %i  quad  %i  %i  %i  %i\n' % (
                cquad_id[i], cquad_tag[i], cquad_1[i], cquad_2[i], cquad_3[i], cquad_4[i]))

        if ctetra > 0:
            for i in range(ctetra):
                fid.write('%i  %i  tetra  %i  %i  %i  %i\n' % (
                ctetra_id[i], ctetra_tag[i], ctetra_1[i], ctetra_2[i], ctetra_3[i], ctetra_4[i]))

        if cpyram > 0:
            for i in range(cpyram):
                fid.write('%i  %i  pyram  %i  %i  %i  %i  %i\n' % (
                cpyram_id[i], cpyram_tag[i], cpyram_1[i], cpyram_2[i], cpyram_3[i], cpyram_4[i], cpyram_5[i]))

        if chexa > 0:
            for i in range(chexa):
                fid.write('%i  %i   hex  %i  %i  %i  %i  %i  %i  %i  %i\n' % (
                chexa_id[i], chexa_tag[i], chexa_1[i], chexa_2[i], chexa_3[i], chexa_4[i], chexa_5[i], chexa_6[i],
                chexa_7[i], chexa_8[i]))

    else:
        if ctria > 0:
            for i in range(ctria):
                fid.write('%i    %i  %i  %i    %i\n' % (ctria_id[i], ctria_1[i], ctria_2[i], ctria_3[i], ctria_tag[i]))

    fid.close()
    return


def read_file_monitor_info(folder_speed, folder_monitor):
    import numpy as np
    file_info = folder_speed + '/' + folder_monitor + '/' + "MONITOR.INFO"
    info = np.loadtxt(file_info)
    T = info[0]  # final time simulation
    dt_s = info[1]
    ndt_monit = info[2]
    dt = ndt_monit * dt_s
    mpi_num_proc = int(info[3])  # number of mpi proc
    mpi_mnt_id = info[4:4 + mpi_num_proc].astype(int)  # id mpi for monitors
    return T, dt, mpi_mnt_id


def define_values_opt(i):
    if i == 0:  # DISPLACEMENT
        ext_in = 'D'
        k_increment = 3
    if i == 1:  # VELOCITY
        ext_in = 'V'
        k_increment = 3
    if i == 2:  # ACCELERATION
        ext_in = 'A'
        k_increment = 3
    if i == 3:  # STRAIN
        ext_in = 'E'
        k_increment = 6
    if i == 4:  # STRESS
        ext_in = 'S'
        k_increment = 6
    if i == 5:  # OMEGA
        ext_in = 'O'
        k_increment = 3

    return ext_in, k_increment


def read_speed_input(folder_speed):
    import numpy as np
    # transform files MONITORXXXXX.D INTO siteXXXXX.d
    # 1 file for each site
    file_speed_input = folder_speed + "/SPEED.input"
    with open(file_speed_input) as fsi:
        line = fsi.readline()
        while line:
            if "MONFILE" in line:
                line_stripped = line.strip("MONFILE")
                folder_monitor_str = line_stripped.replace(' ', '').rstrip('\n')
            if "OPTIOUT" in line:
                line_stripped = line.strip("OPTIOUT").split()
                out_opt_array = np.array(line_stripped, dtype=int)
            line = fsi.readline()

    return folder_monitor_str, out_opt_array

def read_file_monitorX_info(filename2):
    import numpy as np
    info_monitor = np.loadtxt(filename2)
    num_of_mon = int(info_monitor[0])
    id_of_mon = info_monitor[1:1 + num_of_mon].astype(int)
    return num_of_mon, id_of_mon

# def write_monitor(time, values, datafilename):
#     file_id = open(datafilename, "w")
#     size = values.shape[1]
#     for h in range(len(time)):
#         if size == 3:
#             file_id.write('%10.8e   %10.8e   %10.8e  %10.8e \n' % (
#                 time[h], values[h, 0], values[h, 1], values[h, 2]))
#         if size == 6:
#             file_id.write('%10.8e   %10.8e   %10.8e  %10.8e %10.8e   %10.8e   %10.8e \n' % (
#                 time[h], values[h, 0], values[h, 1], values[h, 2], values[h, 3], values[h, 4], values[h, 5]))
#     file_id.close()

def chunks(lst, n):
    """ Yeld successive n-sized chunks from lst"""
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


def write_uscb_format(outfilename, npts, dt, values):
    with open(outfilename, 'w') as f:
        f.write('{:12d}{:17.8E}\n'.format(npts, dt))
        for chunk in chunks(values, 5):
            for value in chunk:
                f.write('{:14.5E}'.format(value))
            f.write('\n')
    return


def read_write_monitor(folder_speed, folder_monitor, mpi_mnt_id, ext_in, k_increment, sites, use_sac):
    import numpy as np
    import os
    from obspy.io.sac import SACTrace

    for i in range(len(mpi_mnt_id)):
        filename1 = folder_speed + '/' + folder_monitor + '/' + "MONITOR" + '{0:05d}'.format(i) + "." + ext_in
        filename2 = folder_speed + '/' + folder_monitor + '/' + "MONITOR" + '{0:05d}'.format(i) + ".INFO"
        if os.path.isfile(filename2):
            num_of_mon, id_of_mon = read_file_monitorX_info(filename2)
            val_monitor = np.loadtxt(filename1)
            k = 0
            for j in range(num_of_mon):
                time = val_monitor[:, 0]

                for ii in range(3):
                    if ii == 0:
                        comp_str = '090' #EW
                        disp = val_monitor[:, 1 + k]
                    elif ii == 1:
                        comp_str = '000' #NS
                        disp = val_monitor[:, 1 + k + 1]
                    else:
                        comp_str = 'ver' #Z
                        disp = val_monitor[:, 1 + k + 2]

                    dt = time[1] - time[0]
                    npts = len(time)
                    vel = compute_derivative(disp, dt)
                    values = vel

                    for iobs in range(len(sites['ID'])):
                        if iobs + 1 == id_of_mon[j]:
                            outfilename = folder_speed + "/" + sites['ID'][iobs] + "." + comp_str + ".gm3D." + \
                                          str(1).zfill(3)

                            if use_sac == 1:
                                sac = SACTrace(npts=npts, delta=dt, nzyear=1970, nzjday=1, nzhour=1, nzmin=0,
                                               nzsec=0, nzmsec=0, data=np.asarray(values))
                                sac.write(outfilename)

                            else:
                                write_uscb_format(outfilename, npts, dt, values)

                k = k + k_increment
    return


def speed2sac(folder, sites):
    folder = folder + '/SPEED'
    use_sac = 1 #sì sac
    folder_monitor, out_opt = read_speed_input(folder)
    T, dt, mpi_mnt_id = read_file_monitor_info(folder, folder_monitor)
    for iopt in range(len(out_opt)):
        ext_in, k_increment = define_values_opt(iopt)
        if out_opt[iopt] == 1:
            read_write_monitor(folder, folder_monitor, mpi_mnt_id, ext_in, k_increment, sites, use_sac)
    return


def speed2ascii(folder, sites):
    folder = folder + '/SPEED'
    use_sac = 0 #no sac
    folder_monitor, out_opt = read_speed_input(folder)
    T, dt, mpi_mnt_id = read_file_monitor_info(folder, folder_monitor)
    for iopt in range(len(out_opt)):
        ext_in, k_increment = define_values_opt(iopt)
        if out_opt[iopt] == 1:
            read_write_monitor(folder, folder_monitor, mpi_mnt_id, ext_in, k_increment, sites, use_sac)
    return









