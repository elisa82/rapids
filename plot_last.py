def create_figure_3plots(selected_code, folder_simulation, j, sites, folder_plot,
                         fmin, fmax, label_code, ext_out, label_map, fault):
    import matplotlib.pylab as plt
    plot_file_map = folder_plot + '/' + label_code + "." + ext_out + ".png"
    plt.figure()
    plt.figure(figsize=(10, 10))
    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = '13'
    plt.rcParams['mathtext.fontset'] = "stix"
    plt.rcParams['legend.title_fontsize'] = '16'
    for k in range(3):
        peak_all = []
        if k == 0:
            label_comp = 'NS'
        if k == 1:
            label_comp = 'EW'
        if k == 2:
            label_comp = 'Z'
        for iobs in range(len(sites['Z'])):
            peaks = retrieve_peaks(selected_code, folder_simulation, iobs, j, sites, folder_plot, fmin, fmax)
            peak_all.append(peaks[k])

        plt.subplot(2, 2, k + 1)
        create_single_map(peak_all, sites, fault, label_map, ext_out)
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


def create_single_map(peak_all, sites, fault, label, ext_out):
    import matplotlib.pylab as plt
    import numpy as np
    from mpl_toolkits.basemap import Basemap
    from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
    from matplotlib.patches import Polygon
    import matplotlib.colors as colors

    aspect = 20
    pad_fraction = 0.5

    lon = sites['lon']
    lat = sites['lat']
    val = peak_all

    minlon, maxlon, minlat, maxlat = define_area_plot(fault, sites)
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
    if ext_out == 'd':
        vmax = 50
        bounds = [8.900e-03, 4.695e-02, 3.430e-01, 1.040e+00, 2.530e+00, 5.450e+00, 1.080e+01, 2.015e+01, 3.585e+01, 61.5]
    if ext_out == 'v':
        vmax = 100
        bounds = np.array([0.0178, 0.0939, 0.686, 2.08, 5.06, 10.9, 21.6, 40.3, 71.7, 123])
    if ext_out == 'a':
        vmax = 1.
        bounds = np.array([0.000555, 0.00232, 0.0121, 0.0338, 0.0746, 0.145, 0.261, 0.444, 0.723, 1.138])
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


def retrieve_peaks(selected_code, folder_simulation, iobs, j, sites, folder_plot, fmin, fmax):
    import numpy as np

    peaks = np.zeros(3)
    time, signal_filtered, fft, xf, peaks[0] = \
        define_selected_time_history(selected_code, folder_simulation, iobs, j, 'NS', sites, folder_plot, fmin, fmax)
    time, signal_filtered, fft, xf, peaks[1] = \
        define_selected_time_history(selected_code, folder_simulation, iobs, j, 'EW', sites, folder_plot, fmin, fmax)
    time, signal_filtered, fft, xf, peaks[2] = \
        define_selected_time_history(selected_code, folder_simulation, iobs, j, 'Z', sites, folder_plot, fmin, fmax)
    return peaks


def plot_selected_code(selected_code, folder_simulation, iobs, desidered_output, sites, folder_plot, label_code, col, fmin, fmax):
    import matplotlib.pylab as plt

    time, signal_filtered, fft, xf, temp = \
        define_selected_time_history(selected_code, folder_simulation, iobs, desidered_output, 'NS', sites, folder_plot, fmin, fmax)
    plt.subplot(321)
    line_name, = plt.plot(time, signal_filtered, col, label=label_code)
    plt.subplot(322)

    plt.loglog(xf, fft, col)

    time, signal_filtered, fft, xf, temp = \
        define_selected_time_history(selected_code, folder_simulation, iobs, desidered_output, 'EW', sites, folder_plot, fmin, fmax)
    plt.subplot(323)
    plt.plot(time, signal_filtered, col)
    plt.subplot(324)
    plt.loglog(xf, fft, col)

    time, signal_filtered, fft, xf, temp = \
        define_selected_time_history(selected_code, folder_simulation, iobs, desidered_output, 'Z', sites, folder_plot, fmin, fmax)
    plt.subplot(325)
    plt.plot(time, signal_filtered, col)
    plt.subplot(326)
    plt.loglog(xf, fft, col)
    return line_name


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


def define_selected_time_history(selected_code, folder_simulation, nobs, desired_output, comp, sites, fmin, fmax):
    import numpy as np
    from obspy.core import read
    from scipy import signal

    if selected_code == 'ucsb' or 'speed':
        if comp == 'NS':
            comp_str = '000'
        if comp == 'EW':
            comp_str = '090'
        if comp == 'Z':
            comp_str = 'ver'
            if selected_code == 'ucsb':
                filename = folder_simulation + "/" + sites['ID'][nobs] + "." + comp_str + ".gm1D." + str(1).zfill(3)
            elif selected_code == 'speed':
                filename = folder_simulation + "/" + sites['ID'][nobs] + "." + comp_str + ".gm3D." + str(1).zfill(3)
            else:
                'Bla'
        st = read(filename)
        time_series = st[0]
        dt = time_series.stats.delta
        time = np.arange(0, time_series.stats.npts) * dt
        sig_vel = time_series.data * 100
        sig_disp = compute_integral(sig_vel, dt)
        sig_acc = compute_derivative(sig_vel, dt)
        # if selected_code == 'speed': questo era x speed in spostamento
        #     sig_disp = time_series.data * 100
        #     sig_vel = compute_derivative(sig_disp, dt)
        #     sig_acc = compute_derivative(sig_vel, dt)

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

    if desired_output == 'dis':
        sig = sig_disp
        sig = sig - np.mean(sig)
        sig = signal.detrend(sig)
    if desired_output == 'vel':
        sig = sig_vel
    if desired_output == 'acc':
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
    order = 2
    signal_filtered = sig
    if fmin is None and fmax is not None:
        signal_filtered = lowpass(sig, fmax, df=1 / dt, corners=order, zerophase=True)
    if fmax is None and fmin is not None:
        signal_filtered = highpass(sig, fmin, df=1 / dt, corners=order, zerophase=True)
    if fmax is not None and fmin is not None:
        signal_filtered = bandpass(sig, freqmin=fmin, freqmax=fmax, df=1 / dt, corners=order, zerophase=True)
    signal_fft, xf = do_fft(signal_filtered, dt)
    return signal_filtered, signal_fft, xf


def plot_seis(folder_speed, folder_hisada, folder_ucsb, folder_plot, plot_param, code, sites, fault, desidered_output):
    import matplotlib.pylab as plt

    fmin = plot_param['fmin_filter']
    fmax = plot_param['fmax_filter']
    time_max_plot = plot_param['time_max_plot']

    Repi = compute_Repi(sites, fault)

    if desidered_output == 'dis':
        ext_out = 'd'
        unit_measure = 'cm'
    if desidered_output == 'vel':
        ext_out = 'v'
        unit_measure = 'cm/s'
    if desidered_output == 'acc':
        ext_out = 'a'
        unit_measure = 'g'
    for iobs in range(len(sites['Z'])):
        handles_tag = []
        plot_file = folder_plot + '/' + str(sites['ID'][iobs]) + "." + ext_out + ".png"
        plt.figure(100 * j + iobs)
        for jj in range(3):
            if jj == 0:
                if 'hisada' in code:
                    selected_code = 'hisada'
                    label_code = 'HisadaBielak'
                    folder_simulation = folder_hisada
                    col = 'r'
                    line_name = plot_selected_code(selected_code, folder_simulation, iobs, j,
                                                       sites, folder_plot, label_code, col, fmin, fmax)
                    handles_tag.append(line_name)
                else:
                    pass
            if jj == 1:
                if 'speed' in code:
                    selected_code = 'speed'
                    folder_simulation = folder_speed
                    label_code = 'SPEED'
                    col = 'b'
                    line_name = plot_selected_code(selected_code, folder_simulation, iobs, desired_output,
                                                       sites, folder_plot, label_code, col, fmin, fmax)
                        handles_tag.append(line_name)
                else:
                    pass
            if jj == 2:
                if 'ucsb' in code:
                    selected_code = 'ucsb'
                    folder_simulation = folder_ucsb
                    label_code = 'UCSB'
                    col = 'g'
                    line_name = plot_selected_code(selected_code, folder_simulation, iobs, desired_output,
                                                       sites, folder_plot, label_code, col, fmin, fmax)
                    handles_tag.append(line_name)
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
        plt.tight_layout()
        plt.savefig(plot_file)
        plt.close('')


def plot_map(folder_speed, folder_hisada, folder_ucsb, folder_plot, plot_param, code, sites, fault, desired_output):
    fmin = plot_param['fmin_filter']
    fmax = plot_param['fmax_filter']

    if desired_output == 'dis':
        ext_out = 'd'
        unit_measure = 'cm'
        label_map = 'PGD (' + unit_measure + ')'
    if desired_output == 'vel':
        ext_out = 'v'
        unit_measure = 'cm/s'
        label_map = 'PGV (' + unit_measure + ')'
    if desired_output == 'acc':
        ext_out = 'a'
        unit_measure = 'g'
        label_map = 'PGA (' + unit_measure + ')'
    for jj in range(3):
        if jj == 0:
            if 'hisada' in code:
                selected_code = 'hisada'
                label_code = 'HisadaBielak'
                folder_simulation = folder_hisada
                create_figure_3plots(selected_code, folder_simulation, desired_output, sites, folder_plot,
                                         fmin, fmax, label_code, ext_out, label_map, fault)
            else:
                pass
        if jj == 1:
            if 'speed' in code:
                selected_code = 'speed'
                folder_simulation = folder_speed
                label_code = 'SPEED'
                create_figure_3plots(selected_code, folder_simulation, desired_output, sites, folder_plot,
                                         fmin, fmax, label_code, ext_out, label_map, fault)
            else:
                pass
        if jj == 2:
            if 'ucsb' in code:
                    selected_code = 'ucsb'
                    folder_simulation = folder_ucsb
                    label_code = 'UCSB'
                    create_figure_3plots(selected_code, folder_simulation, desired_output, sites, folder_plot,
                                         fmin, fmax, label_code, ext_out, label_map, fault)
                else:
                    pass

# https://sourcespec.readthedocs.io/en/latest/_modules/spectrum.html
