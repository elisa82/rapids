#!/usr/bin/python3

def strtoint(sf):
    """

    """
    try:
        x = int(sf)
    except ValueError:
        return None
    return x


def strtofloat(sf):
    """
    """
    try:
        x = float(sf)
    except ValueError:
        return None
    return x

def to_utc_date_time(value):
    """
    """
    from obspy.core import UTCDateTime

    try:
        date, time = value.split('_')
    except ValueError:
        date = value

    year = int(date[0:4])
    month = int(date[4:6])
    day = int(date[6:8])

    hour = int(time[0:2])
    mins = int(time[2:4])
    secs = float(time[4:])

    return UTCDateTime(year, month, day, hour, mins) + secs

def read_esm(filename_in):

    # Import libraries
    from obspy.core import Stats
    import numpy as np

    headers = {}

    # read file
    fh = open(filename_in, 'rt')
    for j in range(64):
        key, value = fh.readline().strip().split(':', 1)
        headers[key.strip()] = value.strip()

    header = Stats()

    header['dyna'] = {}

    header['network'] = headers['NETWORK']
    header['station'] = headers['STATION_CODE']
    header['location'] = headers['LOCATION']
    header['channel'] = headers['STREAM']
    try:
        # use toUTCDateTime to convert from DYNA format
        header['starttime'] \
            = to_utc_date_time(headers
                                ['DATE_TIME_FIRST_SAMPLE_YYYYMMDD_HHMMSS'])
    except ValueError:
        header['starttime'] = to_utc_date_time('19700101_000000')
    header['sampling_rate'] = 1 / float(headers['SAMPLING_INTERVAL_S'])
    header['delta'] = float(headers['SAMPLING_INTERVAL_S'])
    header['npts'] = int(headers['NDATA'])
    header['calib'] = 1  # not in file header

    # DYNA dict float data
    header['dyna']['EVENT_LATITUDE_DEGREE'] = strtofloat(
        headers['EVENT_LATITUDE_DEGREE'])
    header['dyna']['EVENT_LONGITUDE_DEGREE'] = strtofloat(
        headers['EVENT_LONGITUDE_DEGREE'])
    header['dyna']['EVENT_DEPTH_KM'] = strtofloat(headers['EVENT_DEPTH_KM'])
    header['dyna']['HYPOCENTER_REFERENCE'] = headers['HYPOCENTER_REFERENCE']
    header['dyna']['MAGNITUDE_W'] = strtofloat(headers['MAGNITUDE_W'])
    header['dyna']['MAGNITUDE_L'] = strtofloat(headers['MAGNITUDE_L'])
    header['dyna']['STATION_LATITUDE_DEGREE'] = strtofloat(
        headers['STATION_LATITUDE_DEGREE'])
    header['dyna']['STATION_LONGITUDE_DEGREE'] = strtofloat(
        headers['STATION_LONGITUDE_DEGREE'])
    header['dyna']['VS30_M_S'] = strtofloat(headers['VS30_M/S'])
    header['dyna']['EPICENTRAL_DISTANCE_KM'] = strtofloat(
        headers['EPICENTRAL_DISTANCE_KM'])
    header['dyna']['EARTHQUAKE_BACKAZIMUTH_DEGREE'] = strtofloat(
        headers['EARTHQUAKE_BACKAZIMUTH_DEGREE'])
    header['dyna']['DURATION_S'] = strtofloat(headers['DURATION_S'])
    header['dyna']['INSTRUMENTAL_FREQUENCY_HZ'] = strtofloat(
        headers['INSTRUMENTAL_FREQUENCY_HZ'])
    header['dyna']['INSTRUMENTAL_DAMPING'] = strtofloat(
        headers['INSTRUMENTAL_DAMPING'])
    header['dyna']['FULL_SCALE_G'] = strtofloat(headers['FULL_SCALE_G'])

    # data type is acceleration
    if headers['DATA_TYPE'] == "ACCELERATION" \
            or headers['DATA_TYPE'] == "ACCELERATION RESPONSE SPECTRUM":
        header['dyna']['PGA_CM_S_2'] = strtofloat(headers['PGA_CM/S^2'])
        header['dyna']['TIME_PGA_S'] = strtofloat(headers['TIME_PGA_S'])
    # data type is velocity
    if headers['DATA_TYPE'] == "VELOCITY" \
            or headers['DATA_TYPE'] == "PSEUDO-VELOCITY RESPONSE SPECTRUM":
        header['dyna']['PGV_CM_S'] = strtofloat(headers['PGV_CM/S'])
        header['dyna']['TIME_PGV_S'] = strtofloat(headers['TIME_PGV_S'])
    # data type is displacement
    if headers['DATA_TYPE'] == "DISPLACEMENT" \
            or headers['DATA_TYPE'] == "DISPLACEMENT RESPONSE SPECTRUM":
        header['dyna']['PGD_CM'] = strtofloat(headers['PGD_CM'])
        header['dyna']['TIME_PGD_S'] = strtofloat(headers['TIME_PGD_S'])

    header['dyna']['LOW_CUT_FREQUENCY_HZ'] = strtofloat(
        headers['LOW_CUT_FREQUENCY_HZ'])
    header['dyna']['HIGH_CUT_FREQUENCY_HZ'] = strtofloat(
        headers['HIGH_CUT_FREQUENCY_HZ'])

    # DYNA dict int data
    header['dyna']['STATION_ELEVATION_M'] = strtoint(
        headers['STATION_ELEVATION_M'])
    header['dyna']['SENSOR_DEPTH_M'] = strtoint(headers['SENSOR_DEPTH_M'])
    header['dyna']['N_BIT_DIGITAL_CONVERTER'] = strtoint(
        headers['N_BIT_DIGITAL_CONVERTER'])
    header['dyna']['FILTER_ORDER'] = strtoint(headers['FILTER_ORDER'])

    # DYNA dict string data
    header['dyna']['EVENT_NAME'] = headers['EVENT_NAME']
    header['dyna']['EVENT_ID'] = headers['EVENT_ID']
    header['dyna']['EVENT_DATE_YYYYMMDD'] = headers['EVENT_DATE_YYYYMMDD']
    header['dyna']['EVENT_TIME_HHMMSS'] = headers['EVENT_TIME_HHMMSS']
    header['dyna']['MAGNITUDE_W_REFERENCE'] = headers[
        'MAGNITUDE_W_REFERENCE']
    header['dyna']['MAGNITUDE_L_REFERENCE'] = headers[
        'MAGNITUDE_L_REFERENCE']
    header['dyna']['FOCAL_MECHANISM'] = headers['FOCAL_MECHANISM']
    header['dyna']['STATION_NAME'] = headers['STATION_NAME']
    header['dyna']['SITE_CLASSIFICATION_EC8'] = headers[
        'SITE_CLASSIFICATION_EC8']
    header['dyna']['MORPHOLOGIC_CLASSIFICATION'] = headers[
        'MORPHOLOGIC_CLASSIFICATION']
    header['dyna']['DATE_TIME_FIRST_SAMPLE_PRECISION'] = headers[
        'DATE_TIME_FIRST_SAMPLE_PRECISION']
    header['dyna']['UNITS'] = headers['UNITS']
    header['dyna']['INSTRUMENT'] = headers['INSTRUMENT']
    header['dyna']['INSTRUMENT_ANALOG_DIGITAL'] = headers[
        'INSTRUMENT_ANALOG/DIGITAL']
    header['dyna']['BASELINE_CORRECTION'] = headers['BASELINE_CORRECTION']
    header['dyna']['FILTER_TYPE'] = headers['FILTER_TYPE']
    header['dyna']['LATE_NORMAL_TRIGGERED'] = headers[
        'LATE/NORMAL_TRIGGERED']
    header['dyna']['HEADER_FORMAT'] = headers['HEADER_FORMAT']
    header['dyna']['DATABASE_VERSION'] = headers['DATABASE_VERSION']
    header['dyna']['DATA_TYPE'] = headers['DATA_TYPE']
    header['dyna']['PROCESSING'] = headers['PROCESSING']
    header['dyna']['DATA_LICENSE'] = headers['DATA_LICENSE']
    header['dyna']['DATA_TIMESTAMP_YYYYMMDD_HHMMSS'] = headers[
        'DATA_TIMESTAMP_YYYYMMDD_HHMMSS']
    header['dyna']['DATA_CITATION'] = headers['DATA_CITATION']
    header['dyna']['DATA_CREATOR'] = headers['DATA_CREATOR']
    header['dyna']['ORIGINAL_DATA_MEDIATOR_CITATION'] = headers[
        'ORIGINAL_DATA_MEDIATOR_CITATION']
    header['dyna']['ORIGINAL_DATA_MEDIATOR'] = headers[
        'ORIGINAL_DATA_MEDIATOR']
    header['dyna']['ORIGINAL_DATA_CREATOR_CITATION'] = headers[
        'ORIGINAL_DATA_CREATOR_CITATION']
    header['dyna']['ORIGINAL_DATA_CREATOR'] = headers[
        'ORIGINAL_DATA_CREATOR']
    header['dyna']['USER1'] = headers['USER1']
    header['dyna']['USER2'] = headers['USER2']
    header['dyna']['USER3'] = headers['USER3']
    header['dyna']['USER4'] = headers['USER4']
    header['dyna']['USER5'] = headers['USER5']

    # read data
    acc_data = np.loadtxt(fh, dtype='float32')
    fh.close()

    time = []
    for j in range(0, header['npts']):
        t = j * header['delta']
        time.append(t)

    return time, np.asarray(acc_data)

def create_single_map(intensity_measure_all, minlon, maxlon, minlat, maxlat, comp, bounds, fault, label, desired_output, file_topo_plot):
    import matplotlib.pylab as plt
    import numpy as np
    from mpl_toolkits.basemap import Basemap
    from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
    from matplotlib.patches import Polygon
    from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
    from scipy.interpolate import griddata
    import rasterio

    with rasterio.open(file_topo_plot) as src:
        data = src.read(1)  # Legge il primo (e spesso unico) canale/raster
        # Estrai i limiti del raster (latitudine e longitudine)
        lon_min_t, lat_min_t, lon_max_t, lat_max_t = src.bounds
        crs = src.crs  # Sistema di riferimento spaziale
        data = data[::-1]  # Inverte i dati (se necessario)

    # Trova il range di valori nell'intero raster (minimo e massimo)
    min_value_t = -2000
    max_value_t = data.max()

    plt.figure()
    #plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = '10'
    plt.rcParams['mathtext.fontset'] = "stix"
    plt.rcParams['legend.title_fontsize'] = '16'

    aspect = 20
    pad_fraction = 0.5

    lon = intensity_measure_all[:, 1]
    lat = intensity_measure_all[:, 2]
    val = intensity_measure_all[:, 3 + comp]
    val = np.nan_to_num(val, nan=0)
    if desired_output == 2:
        val = val*100.

    # create map
    map = Basemap(projection='merc', resolution='i', llcrnrlon=minlon, llcrnrlat=minlat,
                  urcrnrlon=maxlon, urcrnrlat=maxlat)

    # Calcola i limiti dell'area visibile di Basemap
    map_llcrnrlat, map_llcrnrlon = map.llcrnrlat, map.llcrnrlon
    map_urcrnrlat, map_urcrnrlon = map.urcrnrlat, map.urcrnrlon

    # Estrai i lati in pixel delle coordinate del raster
    lon_t = np.linspace(lon_min_t, lon_max_t, data.shape[1])
    lat_t = np.linspace(lat_min_t, lat_max_t, data.shape[0])

    # Trova l'indice di inizio e fine per latitudine e longitudine
    lat_start = np.argmax(lat_t >= map_llcrnrlat)
    lat_end = np.argmax(lat_t >= map_urcrnrlat)
    lon_start = np.argmax(lon_t >= map_llcrnrlon)
    lon_end = np.argmax(lon_t >= map_urcrnrlon)

    # Ritaglia il raster in base agli indici
    data_cropped = data[lat_start:lat_end, lon_start:lon_end]
    norm_t = Normalize(vmin=min_value_t, vmax=max_value_t)
    map.imshow(data_cropped, origin='lower', cmap='terrain', extent=(lon_t[lon_start], lon_t[lon_end], lat_t[lat_start], lat_t[lat_end]), norm=norm_t, alpha=0.5)

    map.drawparallels(np.arange(-90, 91., 0.5), labels=[1, 0, 0, 1], dashes=[1, 1], linewidth=0.25, color='0.5')
    map.drawmeridians(np.arange(-180., 181., 0.5), labels=[0, 1, 0, 1], dashes=[1, 1], linewidth=0.25, color='0.5')
    map.drawcoastlines()
    map.drawcountries()

    colors = [
        "#FFFFFF",  # Bianco (valori bassissimi o vuoti)
        "#00BFFF",  # Blu chiaro (0.1)
        "#40E0D0",  # Turchese (0.2)
        "#7FFF00",  # Verde chiaro (0.5)
        "#ADFF2F",  # Verde giallastro (1)
        "#FFFF00",  # Giallo (2)
        "#FFD700",  # Oro (5)
        "#FFA500",  # Arancione (10)
        "#FF4500",  # Rosso aranciato (20)
        "#FF0000",  # Rosso (50)
        "#B22222",  # Rosso scuro (100)
        "#800000"   # Marrone rossastro (200)
        ]

    cmap = ListedColormap(colors)
    norm = BoundaryNorm(boundaries=bounds, ncolors=len(colors))
    x, y = map(np.asarray(lon), np.asarray(lat))
    # Bounds PGA e PGV based on Oliveti Faenza Michelini, per PGD ? PGV/2
    #questo fa sistemato. Ora ho usato Faenza e Michelini (2010, 2011) come per le shakemaps
    map.scatter(x, y, s=7, c=val, cmap=cmap, zorder=2, norm=norm, alpha=0.8)

    x1, y1 = map(fault['vertex']['pbl']['lon'], fault['vertex']['pbl']['lat'])
    x2, y2 = map(fault['vertex']['ptl']['lon'], fault['vertex']['ptl']['lat'])
    x3, y3 = map(fault['vertex']['ptr']['lon'], fault['vertex']['ptr']['lat'])
    x4, y4 = map(fault['vertex']['pbr']['lon'], fault['vertex']['pbr']['lat'])
    hypo_lon, hypo_lat = map(fault['hypo']['lon'], fault['hypo']['lat'])
    poly = Polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4)], facecolor='lightyellow', edgecolor='black', linewidth=2, zorder=5, alpha=0.3)
    plt.gca().add_patch(poly)
    map.plot(hypo_lon, hypo_lat, marker='*', markersize=10, c='black')

    ax = plt.gca()
    # Create a regular grid over the map extent
    lon_unique = np.linspace(lon.min(), lon.max(), num=500)
    lat_unique = np.linspace(lat.min(), lat.max(), num=500)
    grid_lon, grid_lat = np.meshgrid(lon_unique, lat_unique)
    grid_val = griddata(
        (lon, lat),   # Punti originali
        val,          # Valori da interpolare
        (grid_lon, grid_lat),  # Griglia di destinazione
        method='linear'  # Metodo di interpolazione: 'linear', 'nearest', 'cubic'
    )
    x_grid, y_grid = map(grid_lon, grid_lat)
    specific_value = 0.1*100 #voglio la linea di contorno a 0.1 g
    contour = ax.contour(x_grid, y_grid, grid_val, levels=[specific_value], colors='red', linewidths=1, zorder=4)
    ax.clabel(contour, inline=True, fontsize=8, fmt='%1.1f')
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1. / aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)
    ticks = bounds
    cb = plt.colorbar(cax=cax, cmap=cmap, boundaries=bounds, ticks=ticks, spacing='uniform', orientation='vertical')
    cb.set_label(label=label, size=16)

    # cb_ymin = min(val)  # minimum value to show on colobar
    # cb_ymax = max(val)  # maximum value to show on colobar
    # cb_xmin, cb_xmax = cb.ax.get_xlim()
    # cb.ax.set_ylim(cb_ymin, cb_ymax)
    # cb.outline.set_visible(False)  # hide the surrounding spines, which are too large after set_ylim
    # cb.ax.add_patch(plt.Rectangle((cb_xmin, cb_ymin), cb_xmax - cb_xmin, cb_ymax - cb_ymin,
    #             
    return


def create_figure_3plots(intensity_measure_all, fault, sites, label_map, desired_output, folder_plot, selected_code, file_topo_plot, ext_out):
    import matplotlib.pylab as plt
    import numpy as np


    #minlon, maxlon, minlat, maxlat = define_area_plot(fault, sites)
    minlon = 11.5
    maxlon = 14.25
    minlat = 45.5
    maxlat =  47.5
    bounds = np.array([0, 0.1, 0.2, 0.5, 1., 2., 5., 10., 20, 50, 100, 200])

    if desired_output == 1:
        #velocity
        #bounds = np.array([0.0178, 0.0939, 0.686, 2.08, 5.06, 10.9, 21.6, 40.3, 71.7, 123])
        bounds = bounds
    if desired_output == 2:
        #acceleration
        #bounds = np.array([0.000555, 0.00232, 0.0121, 0.0338, 0.0746, 0.145, 0.261, 0.444, 0.723, 1.138])
        #bounds = np.array([0, 0.06, 0.21, 0.81, 1.97, 4.82, 11.8, 28.7, 70.1, 171])
        bounds = bounds 
    if desired_output == 0:
        #displacement
        #bounds = [8.900e-03, 4.695e-02, 3.430e-01, 1.040e+00, 2.530e+00, 5.450e+00, 1.080e+01, 2.015e+01, 3.585e+01,
        #          61.5]
        bounds = bounds/2.
    if desired_output == 3:
        vmax = 1.
        bounds = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5])

    for k in range(3):
        if k == 0:
            label_comp = 'NS'
        if k == 1:
            label_comp = 'EW'
        if k == 2:
            label_comp = 'Z'

        plot_file_map = folder_plot + '/' + selected_code + "." + ext_out + '_' + label_comp + ".png"
        create_single_map(intensity_measure_all, minlon, maxlon, minlat, maxlat, k, bounds, fault, label_map, desired_output, file_topo_plot)
        plt.title(label_comp)
        plt.savefig(plot_file_map)
    plt.close()

    return


def define_area_plot(fault, sites):
    buffer = 0.2
    for i in range(2):
        if i == 0:
            str_coord = 'lon'
        else:
            str_coord = 'lat'
        min_val = min(sites[str_coord])
        max_val = max(sites[str_coord])
        min_val = min(min_val, fault['hypo'][str_coord])
        max_val = max(max_val, fault['hypo'][str_coord])

        if fault['vertex']['pbl'][str_coord] > max_val:
            max_val = fault['vertex']['pbl'][str_coord]
        if fault['vertex']['pbr'][str_coord] > max_val:
            max_val = fault['vertex']['pbr'][str_coord]
        if fault['vertex']['ptr'][str_coord] > max_val:
            max_val = fault['vertex']['ptr'][str_coord]
        if fault['vertex']['ptl'][str_coord] > max_val:
            max_val = fault['vertex']['ptl'][str_coord]

        if str_coord == 'lon':
            maxlon = max_val + buffer
        else:
            maxlat = max_val + buffer

        if fault['vertex']['pbl'][str_coord] < min_val:
            min_val = fault['vertex']['pbl'][str_coord]
        if fault['vertex']['pbr'][str_coord] < min_val:
            min_val = fault['vertex']['pbr'][str_coord]
        if fault['vertex']['ptr'][str_coord] < min_val:
            min_val = fault['vertex']['ptr'][str_coord]
        if fault['vertex']['ptl'][str_coord] < min_val:
            min_val = fault['vertex']['ptl'][str_coord]
        if str_coord == 'lon':
            minlon = min_val - buffer
        else:
            minlat = min_val - buffer

    return minlon, maxlon, minlat, maxlat


def plot_selected_code(selected_code, folder_simulation, iobs, desired_output, sites, label_code, col,
                       fmin, fmax, isource):
    import matplotlib.pylab as plt
    import numpy as np

    peaks = np.zeros(3)
    #arias_intensity = np.zeros(3)

    time, signal_filtered, fft, xf, peaks[0] = \
        define_selected_time_history(selected_code, folder_simulation, iobs, desired_output, 'NS', sites,
                                     fmin, fmax, isource)
    plt.subplot(321)
    line_name, = plt.plot(time, signal_filtered, col, label=label_code, linewidth=0.4)
    plt.subplot(322)

    plt.loglog(xf, fft, col, linewidth=0.4)

    time, signal_filtered, fft, xf, peaks[1] = \
        define_selected_time_history(selected_code, folder_simulation, iobs, desired_output, 'EW', sites,
                                     fmin, fmax, isource)
    plt.subplot(323)
    plt.plot(time, signal_filtered, col, linewidth=0.4)
    plt.subplot(324)
    plt.loglog(xf, fft, col, linewidth=0.4)

    time, signal_filtered, fft, xf, peaks[2] = \
        define_selected_time_history(selected_code, folder_simulation, iobs, desired_output, 'Z', sites,
                                     fmin, fmax, isource)
    plt.subplot(325)
    plt.plot(time, signal_filtered, col, linewidth=0.4)
    plt.subplot(326)
    plt.loglog(xf, fft, col, linewidth=0.4)

    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    return line_name, peaks


def compute_integral_2(A, dt):
    import numpy as np
    # Starting from a ground acceleration records calculate
    # ground velocities and displacement from the integrations
    # of the acceleration

    v = np.zeros(A.size)
    #u = np.zeros(A.size)

    # NB le formule implementate prima del 19/10/2005 utilizzavano beta =1/6
    # con le quali si ottengono risultati no stabili per dt maggiori di
    # dt/Tmin<di 0.5513 ==> questo influenza i risultati dell'integrazione
    # di equazioni del moto (linear acceleration method-Newamrk family)

    # %formule di integrazione di Newmark
    gamma = 1. / 2.
    beta = 1. / 4.
    for i in range(1, len(A), 1):
        v[i] = v[i - 1] + (1 - gamma) * dt * A[i - 1] + gamma * dt * A[i]
        #u[i] = u[i - 1] + dt * v[i - 1] + (1. / 2. - beta) * dt ** 2 * A[i - 1] + beta * dt ** 2 * A[i]

    vel = v
    #disp = u

    #return vel, disp
    return vel


def compute_Repi(sites, fault):
    import numpy as np
    Repi = []
    for i in range(len(sites['Z'])):
    #reduced for i in range(10):
        Repi.append(
            np.sqrt((sites['X'][i] - fault['hypo_utm']['X']) ** 2 + (sites['Y'][i] - fault['hypo_utm']['Y']) ** 2))
    return Repi


def compute_derivative(u, dt):
    import numpy as np
    dudt = np.copy(u)
    npts = len(u)
    for j in range(1, npts - 1):
        dudt[j] = (u[j + 1] - u[j - 1]) / (2 * dt)
    dudt[0] = 0.0
    dudt[npts - 1] = dudt[npts - 2]
    return dudt


def compute_integral(A, dt):
    import numpy as np
    gfun = np.copy(A)
    npts = len(A)
    dt05 = 0.5 * dt
    v1 = gfun[0]
    gfun[0] = v1 * dt
    for i in range(1, npts):
        v2 = gfun[i]
        gfun[i] = gfun[i - 1] + dt05 * (v1 + v2)
        v1 = v2
    return gfun


def define_selected_time_history(selected_code, folder_simulation, nobs, desired_output, comp, sites, fmin, fmax,
                                 isource):
    import numpy as np
    from obspy.core import read
    from rapids.conversions import write_uscb_format
    import glob
    from scipy.integrate import cumulative_trapezoid as cumtrapz
    import math

    output_type = 'vel'

    if selected_code == 'msdwn' or selected_code == 'ms' or selected_code == 'dwn':
        isource = 1
        if selected_code == 'msdwn':
            name_run = 'aaa'
        if selected_code == 'dwn':
            name_run = 'aaa'
        if selected_code == 'ms':
            name_run = 'aaa'
        if comp == 'NS':
            comp_str_out = '000'
            comp_str_in = 'sns'
        if comp == 'EW':
            comp_str_out = '090'
            comp_str_in = 'sew'
        if comp == 'Z':
            comp_str_out = 'ver'
            comp_str_in = 'rzz'
        type_msdwn = 'f1'
        filename = folder_simulation + "/R" + str(isource).zfill(4) + '.' + str(nobs+1).zfill(6) + "." + name_run + type_msdwn + '.' + comp_str_in + '.sac'

        tshift = 30

        #tshift = 30 #il 30 va sistemato facendolo leggere dal src del pulsyn
        st = read(filename)
        time_series = st[0]
        dt = time_series.stats.delta
        npts = time_series.stats.npts

        npts = npts - tshift
        time = np.arange(0, npts) * dt

        if comp == 'EW':
            sig_vel = -time_series.data[tshift : npts + tshift]  # cm/s2
        else:
            sig_vel = time_series.data[tshift : npts + tshift]  # cm/s2

        sig_disp = cumtrapz(sig_vel, time, initial=0)
        #sig_disp = compute_integral(sig_vel, dt)
        sig_acc = compute_derivative(sig_vel, dt)

        filename_out = folder_simulation + "/" + sites['ID'][nobs] + "." + comp_str_out + ".gm1D." + str(isource).zfill(3)
        write_uscb_format(filename_out, npts, dt, sig_vel)

    if selected_code == 'ucsb' or selected_code == 'speed' or selected_code == 'stitched-U' or selected_code == 'stitched-SU':
        if comp == 'NS':
            comp_str = '000'
        if comp == 'EW':
            comp_str = '090'
        if comp == 'Z':
            comp_str = 'ver'
        if selected_code == 'ucsb':
            filename = folder_simulation + "/" + sites['ID'][nobs] + "." + comp_str + ".gm1D." + str(isource).zfill(3)
        if selected_code == 'speed':
            filename = folder_simulation + "/" + sites['ID'][nobs] + "." + comp_str + ".gm3D." + str(isource).zfill(3)
        if selected_code == 'stitched-U' or selected_code == 'stitched-SU':
            filename = folder_simulation + "/" + sites['ID'][nobs] + "." + comp_str + ".gmBB." + str(isource).zfill(3)
        time_series = []
        with open(filename, 'r') as f:
            content = f.readlines()
            for x in range(len(content)):
                if x == 0:
                    val = content[x].split()
                    npts = int(val[0])
                    dt = float(val[1])
                else:
                    data = content[x].split()
                    for value in data:
                        a = float(value)
                        time_series.append(a)
        time_series = np.asarray(time_series)
        time = np.arange(0, npts) * dt
        #st = read(filename)
        #time_series = st[0]
        #dt = time_series.stats.delta
        #time = np.arange(0, time_series.stats.npts) * dt
        if output_type == 'acc':
            #sig_acc = time_series.data * 100
            sig_acc = time_series * 100
            sig_vel = cumtrapz(sig_acc, time, initial=0)
            #sig_vel = compute_integral(sig_acc, dt)
            sig_disp = cumtrapz(sig_vel, time, initial=0)
            #sig_disp = compute_integral(sig_vel, dt)
        if output_type == 'vel':
            #sig_vel = time_series.data * 100
            sig_vel = time_series * 100
            sig_disp = cumtrapz(sig_vel, time, initial=0)
            #sig_disp = compute_integral(sig_vel, dt)
            sig_acc = compute_derivative(sig_vel, dt)
        if output_type == 'dis':
            #sig_disp = time_series.data * 100
            sig_disp = time_series * 100
            sig_vel = compute_derivative(sig_disp, dt)
            sig_acc = compute_derivative(sig_vel, dt)

    if selected_code == 'hisada':
        if comp == 'NS':
            comp_str = 'x'
        if comp == 'EW':
            comp_str = 'y'
        if comp == 'Z':
            comp_str = 'z'
        filename = folder_simulation + "/" + comp_str + "_dat." + str(nobs + 1)
        time_series = np.loadtxt(filename)
        time = time_series[:, 0]
        dt = time[1] - time[0]
        sig_vel = time_series[:, 1] * 100
        if comp == 'Z':
            sig_vel = - sig_vel
        sig_acc = compute_derivative(sig_vel, dt)
        sig_disp = cumtrapz(sig_vel, time, initial=0)
        #sig_disp = compute_integral(sig_vel, dt)

    if selected_code == 'esm':
        event = '20240327_0000228'
        tshift = 20
        if comp == 'NS':
            filename_in = glob.glob(folder_simulation+'/*'+sites['ID'][nobs]+'..HHN.D.INT-*.ACC.MP.ASC')
        if comp == 'EW':
            filename_in = glob.glob(folder_simulation+'/*'+sites['ID'][nobs]+'..HHE.D.INT-*.ACC.MP.ASC')
        if comp == 'Z':
            filename_in = glob.glob(folder_simulation+'/*'+sites['ID'][nobs]+'..HHZ.D.INT-*.ACC.MP.ASC')
        time, sig_acc = read_esm(filename_in[0])
        dt = time[1]-time[0]
        for i in range(len(time)):
            time[i] = time[i] - 17
        sig_vel = cumtrapz(sig_acc, time, initial=0)
        #sig_vel = compute_integral(sig_acc, dt)
        sig_disp = cumtrapz(sig_vel, time, initial=0)
        #sig_disp = compute_integral(sig_vel, dt)

    if desired_output == 0:
        sig = sig_disp
        #sig = sig - np.mean(sig)
        #sig = signal.detrend(sig)
    if desired_output == 1:
        sig = sig_vel
    if desired_output == 2:
        sig = sig_acc / 981

    signal_filtered, fft, xf = prepare_signal(sig, dt, fmin, fmax)
    peak = np.max(np.abs(signal_filtered))

    #arias_intensity = 0
    #if desired_output == 2:
    #    signal_acc_cms2 = signal_filtered * 981
    #    integrand = 0
    #    for i in range(len(signal_acc_cms2)):
    #        integrand += signal_acc_cms2[i]**2*dt
    #    arias_intensity = math.pi /(2*981) * integrand

    return time, signal_filtered, fft, xf, peak


def do_fft(sig, delta):
    # https://sourcespec.readthedocs.io/en/latest/_modules/spectrum.html
    import numpy as np
    """Compute the complex Fourier transform of a signal."""
    npts = len(sig)
    # if npts is even, we make it odd
    # so that we do not have a negative frequency in the last point
    # (see numpy.fft.rfft doc)
    if not npts % 2:
        npts -= 1
    fft = np.abs(np.fft.rfft(sig, npts)) * delta
    fftfreq = np.fft.fftfreq(npts, delta)
    fftfreq = fftfreq[0:fft.size]
    return fft, fftfreq


def prepare_signal(sig, dt, fmin, fmax):
    from obspy.signal.filter import bandpass, lowpass, highpass
    order = 4
    signal_filtered = sig
    if fmin is None and fmax is not None:
        signal_filtered = lowpass(sig, fmax, df=1 / dt, corners=order, zerophase=True)
    if fmax is None and fmin is not None:
        signal_filtered = highpass(sig, fmin, df=1 / dt, corners=order, zerophase=True)
    if fmax is not None and fmin is not None:
        signal_filtered = bandpass(sig, freqmin=fmin, freqmax=fmax, df=1 / dt, corners=order, zerophase=True)
    signal_fft, xf = do_fft(signal_filtered, dt)
    return signal_filtered, signal_fft, xf


def create_file_intensity_measure_for_each_seismogram(iobs, sites, ext_out, selected_code, intensity_measure, output_folder):
    import numpy as np
    single_peaks = np.zeros((6))
    single_peaks[0] = iobs+1
    single_peaks[1] = sites['lon'][iobs]
    single_peaks[2] = sites['lat'][iobs]
    single_peaks[3: 6] = intensity_measure
    file_single_peaks = output_folder + '/' + str(iobs) + '_' + ext_out + '_' + selected_code + '.npy'
    np.save(file_single_peaks, single_peaks)
    return


def create_waveforms(folder_plot, plot_param, code, sites, fault, computational_param, path_data, output_folder, rank=0, nranks=1):
    import matplotlib.pylab as plt
    import numpy as np

    folder_hisada = output_folder + '/HISADA'
    folder_ucsb = output_folder + '/UCSB/HF'
    folder_speed = output_folder + '/SPEED'
    folder_stitched_ucsb = output_folder + '/UCSB/STITCHED'
    folder_stitched_speeducsb = output_folder + '/STITCHED'
    folder_msdwn = output_folder + '/MS-DWN'
    folder_ms = output_folder + '/MS'
    folder_dwn = output_folder + '/DWN'
    folder_esm = path_data + '/ESM'

    fmin = plot_param['fmin_filter']
    fmax = plot_param['fmax_filter']
    time_max_plot = plot_param['time_max_plot']
    #output_type = computational_param['output_type'] #ora output_type viene definito direttamente dentro il codice

    Repi = compute_Repi(sites, fault)

    Nsites = len(sites['Z'])
    #reduced Nsites = 10
    Site_Indexes=np.arange(Nsites)

    for iobs in Site_Indexes[rank::nranks]:
        print("rank %d works on site %d"  %(rank,iobs))
        for j in range(3):
            if j == 0:
                ext_out = 'pgd'
                unit_measure = 'cm'
            if j == 1:
                ext_out = 'pgv'
                unit_measure = 'cm/s'
            if j == 2:
                ext_out = 'pga'
                unit_measure = 'g'
                #ext_out2 = 'ai'
                #unit_measure = 'cm/2'
            if fault['IDx'] == 'Yoffe-DCF':
                num_realizations = 1
            else:
                num_realizations = computational_param['realizations']
            for isource in range(1, num_realizations+1):
                handles_tag = []
                plot_file = folder_plot + '/' + str(sites['ID'][iobs]) + "_" + code + "." + ext_out + '_' + \
                        str(isource).zfill(3) + ".png"
                plt.figure(100 * j + iobs)

                if 'hisada' in code:
                    selected_code = 'hisada'
                    label_code = 'HisadaBielak'
                    folder_simulation = folder_hisada
                    col = 'g'
                    if computational_param['fmax_hisada'] < fmax:
                        fmax_hisada = computational_param['fmax_hisada']
                    else:
                        fmax_hisada = fmax
                    line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                      sites, label_code, col, fmin, fmax_hisada, isource)
                    handles_tag.append(line_name)
                    create_file_intensity_measure_for_each_seismogram(iobs, sites, ext_out, selected_code, peaks, output_folder)


                if 'speed' in code:
                    selected_code = 'speed'
                    folder_simulation = folder_speed
                    label_code = 'SPEED'
                    col = 'b'
                    if fmax is not None and computational_param['fmax_speed'] < fmax:
                        fmax_speed = computational_param['fmax_speed']
                    else:
                        fmax_speed = fmax
                    line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                      sites, label_code, col, fmin, fmax_speed, isource)
                    handles_tag.append(line_name)
                    create_file_intensity_measure_for_each_seismogram(iobs, sites, ext_out, selected_code, peaks, output_folder)

                if 'ucsb' in code:
                    selected_code = 'ucsb'
                    folder_simulation = folder_ucsb
                    label_code = 'UCSB'
                    col = 'r'
                    if fmax is not None and computational_param['fmax_ucsb'] < fmax:
                        fmax_ucsb = computational_param['fmax_ucsb']
                    else:
                        fmax_ucsb = fmax
                    line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                      sites, label_code, col, fmin, fmax_ucsb, isource)
                    handles_tag.append(line_name)
                    create_file_intensity_measure_for_each_seismogram(iobs, sites, ext_out, selected_code, peaks, output_folder)

                if 'stitched-U' in code:
                    selected_code = 'stitched-U'
                    folder_simulation = folder_stitched_ucsb
                    label_code = 'STITCHED-UCSB'
                    col = 'k'
                    if fmax is not None and computational_param['fmax_ucsb'] < fmax:
                        fmax_stitched_ucsb = computational_param['fmax_ucsb'] #stitched è uguale ad ucsb
                    else:
                        fmax_stitched_ucsb = fmax
                    line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                      sites, label_code, col, fmin, fmax_stitched_ucsb, isource)
                    handles_tag.append(line_name)
                    create_file_intensity_measure_for_each_seismogram(iobs, sites, ext_out, selected_code, peaks, output_folder)

                if 'stitched-SU' in code:
                    selected_code = 'stitched-SU'
                    folder_simulation = folder_stitched_speeducsb
                    label_code = 'STITCHED-SPEEDUCSB'
                    col = 'k'
                    if fmax is not None and computational_param['fmax_ucsb'] < fmax:
                        fmax_stitched_speeducsb = computational_param['fmax_ucsb'] #stitched è uguale ad ucsb
                    else:
                        fmax_stitched_speeducsb = fmax
                    line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                      sites, label_code, col, fmin, fmax_stitched_speeducsb, isource)
                    handles_tag.append(line_name)
                    create_file_intensity_measure_for_each_seismogram(iobs, sites, ext_out, selected_code, peaks, output_folder)

                if 'hybridmd' in code:
                    selected_code = 'msdwn'
                    folder_simulation = folder_msdwn
                    label_code = 'MS-DWN'
                    col = 'k'
                    fmax_msdwn = fmax
                    line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                      sites, label_code, col, fmin, fmax_msdwn, isource)
                    handles_tag.append(line_name)
                    create_file_intensity_measure_for_each_seismogram(iobs, sites, ext_out, selected_code, peaks, output_folder)

                if 'ms' in code:
                    selected_code = 'ms'
                    folder_simulation = folder_ms
                    label_code = 'MS'
                    col = 'cyan'
                    fmax_msdwn = fmax
                    line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                      sites, label_code, col, fmin, fmax_msdwn, isource)
                    handles_tag.append(line_name)
                    create_file_intensity_measure_for_each_seismogram(iobs, sites, ext_out, selected_code, peaks, output_folder)

                if 'dwn' in code:
                    selected_code = 'dwn'
                    folder_simulation = folder_dwn
                    label_code = 'DWN'
                    col = 'm'
                    fmax_msdwn = fmax
                    line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                      sites, label_code, col, fmin, fmax_msdwn, isource)
                    handles_tag.append(line_name)
                    create_file_intensity_measure_for_each_seismogram(iobs, sites, ext_out, selected_code, peaks, output_folder)

                if 'esm' in code:
                    selected_code = 'esm'
                    label_code = 'recorded'
                    folder_simulation = folder_esm
                    col = 'r'
                    fmax_esm = fmax
                    line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                      sites, label_code, col, fmin, fmax_esm, isource)
                    handles_tag.append(line_name)
                    create_file_intensity_measure_for_each_seismogram(iobs, sites, ext_out, selected_code, peaks, output_folder)

                plt.subplot(321)
                plt.xlim(0, time_max_plot)
                plt.ylabel('NS (' + unit_measure + ')')
                title_lab = 'Repi = {:.2f} km'.format(Repi[iobs] / 1000)
                plt.title(title_lab)
                plt.grid(True)
                plt.subplot(323)
                plt.xlim(0, time_max_plot)
                plt.ylabel('EW (' + unit_measure + ')')
                plt.grid(True)
                plt.subplot(325)
                plt.xlim(0, time_max_plot)
                plt.xlabel("time (s)")
                plt.ylabel('Z (' + unit_measure + ')')
                plt.grid(True)
                plt.subplot(322)
                plt.xlim(fmin, fmax)
                plt.grid(True)
                plt.subplot(324)
                plt.xlim(fmin, fmax)
                plt.grid(True)
                plt.subplot(326)
                plt.xlim(fmin, fmax)
                plt.xlabel("frequency (Hz)")
                plt.legend(handles=handles_tag)
                plt.grid(True)
                plt.savefig(plot_file)
                plt.close('')
    return

def create_maps_for_each_code(sites, ext_out, selected_code, folder_plot, output_folder, fault, label_map, j, file_topo_plot):
    import os
    import numpy as np

    fmt = '%d', '%7.4f', '%7.4f', '%9.5f', '%9.5f', '%9.5f'

    intensity_measure_code = np.zeros((len(sites['ID']), 6))
    #reduced intensity_measure_code = np.zeros((10, 6))
    for iobs in range(len(sites['Z'])):
    #reduced for iobs in range(10):
        file_single_intensity_measure = output_folder + '/' + str(iobs) + '_' + ext_out + '_' + selected_code + '.npy'  
        single_site = np.load(file_single_intensity_measure)
        intensity_measure_code[iobs,:] = single_site 

    file_intensity_measure = folder_plot + '/' + ext_out + '_' +selected_code + '.txt'
    np.savetxt(file_intensity_measure, intensity_measure_code, fmt=fmt)
    create_figure_3plots(intensity_measure_code, fault, sites, label_map, j, folder_plot, selected_code, file_topo_plot, ext_out)
    command_rm_npy = 'rm -f *.npy'
    #os.system(command_rm_npy)
    
    return


def create_maps(folder_plot, sites, code, output_folder, fault, file_topo_plot):
    for j in range(3):
        if j == 0:
            ext_out = 'pgd'
            unit_measure = 'cm'
            label_map = 'PGD (' + unit_measure + ')'
        if j == 1:
            ext_out = 'pgv'
            unit_measure = 'cm/s'
            label_map = 'PGV (' + unit_measure + ')'
        if j == 2:
            ext_out = 'pga'
            unit_measure = '%g'
            label_map = 'PGA (' + unit_measure + ')'
        #if j == 4:
            #ext_out = 'ai'
            #unit_measure = 'cm/s'
            #label_map = 'AI (' + unit_measure + ')'

        #site_number lon lat peak_NS peak_EW peak_Z
        if 'hisada' in code:
            selected_code = 'hisada'
            create_maps_for_each_code(sites,ext_out,selected_code,folder_plot, output_folder, fault, label_map, j, file_topo_plot)

        if 'speed' in code:
            selected_code = 'speed'
            create_maps_for_each_code(sites,ext_out,selected_code,folder_plot, output_folder, fault, label_map, j, file_topo_plot)

        if 'ucsb' in code:
            selected_code = 'ucsb'
            create_maps_for_each_code(sites,ext_out,selected_code,folder_plot, output_folder, fault, label_map, j, file_topo_plot)

        if 'stitched-U' in code:
            selected_code = 'stitched-U'
            create_maps_for_each_code(sites,ext_out,selected_code,folder_plot, output_folder, fault, label_map, j, file_topo_plot)

        if 'stitched-SU' in code:
            selected_code = 'stitched-SU'
            create_maps_for_each_code(sites,ext_out,selected_code,folder_plot, output_folder, fault, label_map, j, file_topo_plot)

        if 'hybridmd' in code:
            selected_code = 'msdwn'
            create_maps_for_each_code(sites,ext_out,selected_code,folder_plot, output_folder, fault, label_map, j, file_topo_plot)

        if 'ms' in code:
            selected_code = 'ms'
            create_maps_for_each_code(sites,ext_out,selected_code,folder_plot, output_folder, fault, label_map, j, file_topo_plot)

        if 'dwn' in code:
            selected_code = 'dwn'
            create_maps_for_each_code(sites,ext_out,selected_code,folder_plot, output_folder, fault, label_map, j, file_topo_plot)

        if 'esm' in code:
            selected_code = 'esm'
            create_maps_for_each_code(sites,ext_out,selected_code,folder_plot, output_folder, fault, label_map, j, file_topo_plot)
    return

def post_processing(output_folder, plot_param, code, sites, fault, computational_param, path_data, comm=None):
    import os

    if comm is None:
        rank=0
        nranks=1
        isParallel = False
    else:
        rank= comm.Get_rank()
        nranks = comm.size
        isParallel = True

    folder_plot = output_folder + '/PLOT'
    folder_waveforms = folder_plot + '/WAVEFORMS'
    folder_maps = folder_plot + '/MAPS'

    if rank == 0:
        os.system("mkdir -p " + folder_plot)
        os.system("mkdir -p " + folder_waveforms)
        os.system("mkdir -p " + folder_maps)

    if isParallel : comm.Barrier()
    create_waveforms(folder_waveforms, plot_param, code, sites, fault, computational_param, path_data, output_folder, rank, nranks)
    if isParallel : comm.Barrier()
    if rank == 0 : 
        file_topo_plot = path_data + '/topo/exportImage.tiff'
        create_maps(folder_maps, sites, code, output_folder, fault, file_topo_plot)

    return
# https://sourcespec.readthedocs.io/en/latest/_modules/spectrum.html
