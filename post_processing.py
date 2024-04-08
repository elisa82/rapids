#!/usr/bin/python3
def create_single_map(peak_all, minlon, maxlon, minlat, maxlat, comp, bounds, fault, label):
    import matplotlib.pylab as plt
    import numpy as np
    from mpl_toolkits.basemap import Basemap
    from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
    from matplotlib.patches import Polygon
    import matplotlib.colors as colors

    aspect = 20
    pad_fraction = 0.5

    lon = peak_all[:, 1]
    lat = peak_all[:, 2]
    val = peak_all[:, 3 + comp]

    # create map
    map = Basemap(projection='merc', resolution='i', llcrnrlon=minlon, llcrnrlat=minlat,
                  urcrnrlon=maxlon, urcrnrlat=maxlat)
    map.drawparallels(np.arange(-90, 91., 0.1), labels=[1, 0, 0, 1], dashes=[1, 1], linewidth=0.25, color='0.5')
    map.drawmeridians(np.arange(-180., 181., 0.1), labels=[0, 1, 0, 1], dashes=[1, 1], linewidth=0.25, color='0.5')
    map.drawmapboundary(fill_color='lightblue')
    map.fillcontinents(color='white', lake_color='lightblue')
    map.drawcoastlines()
    cmap = plt.cm.get_cmap('gist_rainbow_r')

    x, y = map(np.asarray(lon), np.asarray(lat))
    # Bounds PGA e PGV da Oliveti Faenza Michelini, per PGD è PGV/2
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256, extend='both')
    map.scatter(x, y, s=50, c=val, cmap=cmap, zorder=2, norm=norm)
    x1, y1 = map(fault['vertex']['pbl']['lon'], fault['vertex']['pbl']['lat'])
    x2, y2 = map(fault['vertex']['ptl']['lon'], fault['vertex']['ptl']['lat'])
    x3, y3 = map(fault['vertex']['ptr']['lon'], fault['vertex']['ptr']['lat'])
    x4, y4 = map(fault['vertex']['pbr']['lon'], fault['vertex']['pbr']['lat'])
    hypo_lon, hypo_lat = map(fault['hypo']['lon'], fault['hypo']['lat'])
    map.plot(hypo_lon, hypo_lat, marker='*', markersize=15, c='black')
    poly = Polygon([(x1, y1), (x2, y2), (x3, y3), (x4, y4)], facecolor='lightyellow', edgecolor='black', linewidth=2)
    plt.gca().add_patch(poly)

    ax = plt.gca()
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1. / aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)
    cb = plt.colorbar(cax=cax, cmap=cmap, spacing='uniform', orientation='vertical', extend='both')
    cb.set_label(label=label, size=16)

    # cb_ymin = min(val)  # minimum value to show on colobar
    # cb_ymax = max(val)  # maximum value to show on colobar
    # cb_xmin, cb_xmax = cb.ax.get_xlim()
    # cb.ax.set_ylim(cb_ymin, cb_ymax)
    # cb.outline.set_visible(False)  # hide the surrounding spines, which are too large after set_ylim
    # cb.ax.add_patch(plt.Rectangle((cb_xmin, cb_ymin), cb_xmax - cb_xmin, cb_ymax - cb_ymin,
    #                                fc='none', ec='black', clip_on=False))
    return


def create_figure_3plots(peak_all, fault, sites, label_map, desired_output, plot_file_map):
    import matplotlib.pylab as plt
    import numpy as np

    plt.figure()
    plt.figure(figsize=(10, 10))
    #plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = '10'
    plt.rcParams['mathtext.fontset'] = "stix"
    plt.rcParams['legend.title_fontsize'] = '16'

    minlon, maxlon, minlat, maxlat = define_area_plot(fault, sites)

    if desired_output == 0:
        vmax = 50
        bounds = [8.900e-03, 4.695e-02, 3.430e-01, 1.040e+00, 2.530e+00, 5.450e+00, 1.080e+01, 2.015e+01, 3.585e+01,
                  61.5]
    if desired_output == 1:
        vmax = 100
        bounds = np.array([0.0178, 0.0939, 0.686, 2.08, 5.06, 10.9, 21.6, 40.3, 71.7, 123])
    if desired_output == 2:
        vmax = 1.
        bounds = np.array([0.000555, 0.00232, 0.0121, 0.0338, 0.0746, 0.145, 0.261, 0.444, 0.723, 1.138])

    for k in range(3):
        if k == 0:
            label_comp = 'NS'
        if k == 1:
            label_comp = 'EW'
        if k == 2:
            label_comp = 'Z'

        plt.subplot(2, 2, k + 1)
        create_single_map(peak_all, minlon, maxlon, minlat, maxlat, k, bounds, fault, label_map)
        plt.title(label_comp)
    plt.tight_layout()
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

    peaks = np.zeros(4)

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

    peaks[3] = np.sqrt(peaks[0]**2+peaks[1]**2)

    return line_name, peaks


def integr_accel(A, dt):
    import numpy as np
    # Starting from a ground acceleration records calculate
    # ground velocities and displacement from the integrations
    # of the acceleration

    v = np.zeros(A.size)
    u = np.zeros(A.size)

    # NB le formule implementate prima del 19/10/2005 utilizzavano beta =1/6
    # con le quali si ottengono risultati no stabili per dt maggiori di
    # dt/Tmin<di 0.5513 ==> questo influenza i risultati dell'integrazione
    # di equazioni del moto (linear acceleration method-Newamrk family)

    # %formule di integrazione di Newmark
    gamma = 1. / 2.
    beta = 1. / 4.
    for i in range(1, len(A), 1):
        v[i] = v[i - 1] + (1 - gamma) * dt * A[i - 1] + gamma * dt * A[i]
        u[i] = u[i - 1] + dt * v[i - 1] + (1. / 2. - beta) * dt ** 2 * A[i - 1] + beta * dt ** 2 * A[i]

    vel = v
    disp = u

    return vel, disp


def compute_Repi(sites, fault):
    import numpy as np
    Repi = []
    for i in range(len(sites['Z'])):
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

        sig_disp = compute_integral(sig_vel, dt)
        sig_acc = compute_derivative(sig_vel, dt)

        filename_out = folder_simulation + "/" + sites['ID'][nobs] + "." + comp_str_out + ".gm1D." + str(isource).zfill(3)
        write_uscb_format(filename_out, npts, dt, sig_vel)

    if selected_code == 'ucsb' or selected_code == 'speed' or selected_code == 'stitched-ucsb' or selected_code == 'stitched-speeducsb':
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
        if selected_code == 'stitched-ucsb' or selected_code == 'stitched-speeducsb':
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
            sig_vel = compute_integral(sig_acc, dt)
            sig_disp = compute_integral(sig_vel, dt)
        if output_type == 'vel':
            #sig_vel = time_series.data * 100
            sig_vel = time_series * 100
            sig_disp = compute_integral(sig_vel, dt)
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
        sig_disp = compute_integral(sig_vel, dt)

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

def post_processing(output_folder, plot_param, code, sites, fault, computational_param):
    import matplotlib.pylab as plt
    import numpy as np
    import os

    folder_plot = output_folder + '/plot'
    isExist = os.path.exists(folder_plot)
    if not isExist:
        os.makedirs(folder_plot)

    folder_hisada = output_folder + '/HISADA'
    folder_ucsb = output_folder + '/UCSB'
    folder_speed = output_folder + '/SPEED'
    folder_stitched_ucsb = output_folder + 'UCSB/STITCHED'
    folder_stitched_speeducsb = output_folder + '/STITCHED'
    folder_msdwn = output_folder + '/MS-DWN'
    folder_ms = output_folder + '/MS'
    folder_dwn = output_folder + '/DWN'

    fmin = plot_param['fmin_filter']
    fmax = plot_param['fmax_filter']
    time_max_plot = plot_param['time_max_plot']
    #output_type = computational_param['output_type'] #ora output_type viene definito direttamente dentro il codice

    Repi = compute_Repi(sites, fault)

    fmt = '%d', '%7.4f', '%7.4f', '%9.5f', '%9.5f', '%9.5f', '%9.5f'

    for j in range(3):
        peaks_hisada = np.zeros((len(sites['ID']), 7))
        peaks_speed = np.zeros((len(sites['ID']), 7))
        peaks_ucsb = np.zeros((len(sites['ID']), 7))
        peaks_stitched_ucsb = np.zeros((len(sites['ID']), 7))
        peaks_stitched_speeducsb = np.zeros((len(sites['ID']), 7))
        peaks_msdwn = np.zeros((len(sites['ID']), 7))
        peaks_ms = np.zeros((len(sites['ID']), 7))
        peaks_dwn = np.zeros((len(sites['ID']), 7))
        if j == 0:
            ext_out = 'd'
            unit_measure = 'cm'
            label_map = 'PGD (' + unit_measure + ')'
        if j == 1:
            ext_out = 'v'
            unit_measure = 'cm/s'
            label_map = 'PGV (' + unit_measure + ')'
        if j == 2:
            ext_out = 'a'
            unit_measure = 'g'
            label_map = 'PGA (' + unit_measure + ')'
        #for iobs in range(1):
        for iobs in range(len(sites['Z'])):
            for isource in range(1, computational_param['realizations']+1):
                handles_tag = []
                plot_file = folder_plot + '/' + str(sites['ID'][iobs]) + "_" + code + "." + ext_out + '_' + \
                        str(isource).zfill(3) + ".png"
                plt.figure(100 * j + iobs)
                for jj in range(8):
                    if jj == 0:
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
                            peaks_hisada[iobs, 0] = iobs+1
                            peaks_hisada[iobs, 1] = sites['lon'][iobs]
                            peaks_hisada[iobs, 2] = sites['lat'][iobs]
                            peaks_hisada[iobs, 3: 7] = peaks
                        else:
                            pass
                    if jj == 1:
                        if 'speed' in code and 'stitched' not in code:
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
                            peaks_speed[iobs, 3: 7] = peaks
                            peaks_speed[iobs, 0] = iobs+1
                            peaks_speed[iobs, 1] = sites['lon'][iobs]
                            peaks_speed[iobs, 2] = sites['lat'][iobs]
                        else:
                            pass
                    if jj == 2:
                        if 'ucsb' in code and 'stitched' not in code:
                            selected_code = 'ucsb'
                            folder_simulation = folder_ucsb
                            label_code = 'UCSB'
                            col = 'g'
                            if fmax is not None and computational_param['fmax_ucsb'] < fmax:
                                fmax_ucsb = computational_param['fmax_ucsb']
                            else:
                                fmax_ucsb = fmax
                            line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                              sites, label_code, col, fmin, fmax_ucsb, isource)
                            handles_tag.append(line_name)
                            peaks_ucsb[iobs, 3: 7] = peaks
                            peaks_ucsb[iobs, 0] = iobs+1
                            peaks_ucsb[iobs, 1] = sites['lon'][iobs]
                            peaks_ucsb[iobs, 2] = sites['lat'][iobs]
                        else:
                            pass
                    if jj == 3:
                        if 'stitched-ucsb' in code:
                            selected_code = 'stitched-ucsb'
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
                            peaks_stitched_ucsb[iobs, 3: 7] = peaks
                            peaks_stitched_ucsb[iobs, 0] = iobs+1
                            peaks_stitched_ucsb[iobs, 1] = sites['lon'][iobs]
                            peaks_stitched_ucsb[iobs, 2] = sites['lat'][iobs]
                        else:
                            pass
                    if jj == 4:
                        if 'stitched-speeducsb' in code:
                            selected_code = 'stitched-speeducsb'
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
                            peaks_stitched_speeducsb[iobs, 3: 7] = peaks
                            peaks_stitched_speeducsb[iobs, 0] = iobs+1
                            peaks_stitched_speeducsb[iobs, 1] = sites['lon'][iobs]
                            peaks_stitched_speeducsb[iobs, 2] = sites['lat'][iobs]
                        else:
                            pass
                    if jj == 5:
                        if 'hybridmd' in code:
                            selected_code = 'msdwn'
                            folder_simulation = folder_msdwn
                            label_code = 'MS-DWN'
                            col = 'k'
                            fmax_msdwn = fmax
                            line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                              sites, label_code, col, fmin, fmax_msdwn, isource)
                            handles_tag.append(line_name)
                            peaks_msdwn[iobs, 3: 7] = peaks
                            peaks_msdwn[iobs, 0] = iobs+1
                            peaks_msdwn[iobs, 1] = sites['lon'][iobs]
                            peaks_msdwn[iobs, 2] = sites['lat'][iobs]
                        else:
                            pass

                    if jj == 6:
                        if 'ms' in code:
                            selected_code = 'ms'
                            folder_simulation = folder_ms
                            label_code = 'MS'
                            col = 'cyan'
                            fmax_msdwn = fmax
                            line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                              sites, label_code, col, fmin, fmax_msdwn, isource)
                            handles_tag.append(line_name)
                            peaks_ms[iobs, 3: 7] = peaks
                            peaks_ms[iobs, 0] = iobs+1
                            peaks_ms[iobs, 1] = sites['lon'][iobs]
                            peaks_ms[iobs, 2] = sites['lat'][iobs]
                        else:
                            pass

                    if jj == 7:
                        if 'dwn' in code:
                            selected_code = 'dwn'
                            folder_simulation = folder_dwn
                            label_code = 'DWN'
                            col = 'm'
                            fmax_msdwn = fmax
                            line_name, peaks = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                              sites, label_code, col, fmin, fmax_msdwn, isource)
                            handles_tag.append(line_name)
                            peaks_dwn[iobs, 3: 7] = peaks
                            peaks_dwn[iobs, 0] = iobs+1
                            peaks_dwn[iobs, 1] = sites['lon'][iobs]
                            peaks_dwn[iobs, 2] = sites['lat'][iobs]
                        else:
                            pass


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

        for icode in range(8):
            #site_number lon lat peak_NS peak_EW peak_Z
            if icode == 0:
                if 'hisada' in code:
                    selected_code = 'hisada'
                    file_peaks = folder_plot + '/peaks_'+selected_code + "." + ext_out +'.txt'
                    np.savetxt(file_peaks, peaks_hisada, fmt=fmt)
                    plot_file_map = folder_plot + '/' + selected_code + "." + ext_out + ".png"
                    create_figure_3plots(peaks_hisada, fault, sites, label_map, j, plot_file_map)
                else:
                    pass

            if icode == 1:
                if 'speed' in code and 'stitched' not in code:
                    selected_code = 'speed'
                    file_peaks = folder_plot + '/peaks_'+selected_code + "." + ext_out +'.txt'
                    np.savetxt(file_peaks, peaks_speed, fmt=fmt)
                    plot_file_map = folder_plot + '/' + selected_code + "." + ext_out + ".png"
                    create_figure_3plots(peaks_speed, fault, sites, label_map, j, plot_file_map)
                else:
                    pass

            if icode == 2:
                if 'ucsb' in code and 'stitched' not in code:
                    selected_code = 'ucsb'
                    file_peaks = folder_plot + '/peaks_'+ selected_code + "." + ext_out + '.txt'
                    np.savetxt(file_peaks, peaks_ucsb, fmt=fmt)
                    plot_file_map = folder_plot + '/' + selected_code + "." + ext_out + ".png"
                    create_figure_3plots(peaks_ucsb, fault, sites, label_map, j, plot_file_map)
                else:
                    pass

            if icode == 3:
                if 'stitched-ucsb' in code:
                    selected_code = 'stitched-ucsb'
                    file_peaks = folder_plot + '/peaks_'+ selected_code + "." + ext_out + '.txt'
                    np.savetxt(file_peaks, peaks_stitched_ucsb, fmt=fmt)
                    plot_file_map = folder_plot + '/' + selected_code + "." + ext_out + ".png"
                    create_figure_3plots(peaks_stitched_ucsb, fault, sites, label_map, j, plot_file_map)
                else:
                    pass

            if icode == 4:
                if 'stitched-speeducsb' in code:
                    selected_code = 'stitched-speeducsb'
                    file_peaks = folder_plot + '/peaks_'+ selected_code + "." + ext_out + '.txt'
                    np.savetxt(file_peaks, peaks_stitched_speeducsb, fmt=fmt)
                    plot_file_map = folder_plot + '/' + selected_code + "." + ext_out + ".png"
                    create_figure_3plots(peaks_stitched_speeducsb, fault, sites, label_map, j, plot_file_map)
                else:
                    pass

            if icode == 5:
                if 'hybridmd' in code:
                    selected_code = 'msdwn'
                    file_peaks = folder_plot + '/peaks_'+ selected_code + "." + ext_out + '.txt'
                    np.savetxt(file_peaks, peaks_msdwn, fmt=fmt)
                    plot_file_map = folder_plot + '/' + selected_code + "." + ext_out + ".png"
                    create_figure_3plots(peaks_msdwn, fault, sites, label_map, j, plot_file_map)
                else:
                    pass

            if icode == 6:
                if 'ms' in code:
                    selected_code = 'ms'
                    file_peaks = folder_plot + '/peaks_'+ selected_code + "." + ext_out + '.txt'
                    np.savetxt(file_peaks, peaks_ms, fmt=fmt)
                    plot_file_map = folder_plot + '/' + selected_code + "." + ext_out + ".png"
                    create_figure_3plots(peaks_ms, fault, sites, label_map, j, plot_file_map)
                else:
                    pass

            if icode == 7:
                if 'dwn' in code:
                    selected_code = 'dwn'
                    file_peaks = folder_plot + '/peaks_'+ selected_code + "." + ext_out + '.txt'
                    np.savetxt(file_peaks, peaks_dwn, fmt=fmt)
                    plot_file_map = folder_plot + '/' + selected_code + "." + ext_out + ".png"
                    create_figure_3plots(peaks_dwn, fault, sites, label_map, j, plot_file_map)
                else:
                    pass


# https://sourcespec.readthedocs.io/en/latest/_modules/spectrum.html
