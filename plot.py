def plot_slip_from_file(folder, fault, slip_vector, rptm):
    import matplotlib.pyplot as plt
    import numpy as np

    x = (np.arange(0, fault['number_subfaults_strike']) * fault['subfault_length'])
    y = (np.arange(0, fault['number_subfaults_dip']) * fault['subfault_width'])
    X, Y = np.meshgrid(x, y)

    f, ax = plt.subplots(1, 1)

    # plotting
    f, ax = plt.subplots(1, 1)
    slip_matrix = slip_vector.reshape(fault['number_subfaults_dip'], fault['number_subfaults_strike'])
    rptm = rptm.reshape(fault['number_subfaults_dip'], fault['number_subfaults_strike'])  # rupture time in s

    A = ax.imshow(slip_matrix[::-1], cmap='gist_heat_r', extent=(0, fault['length'], 0, fault['width']), interpolation = 'nearest')
    B = plt.colorbar(A, orientation='vertical', shrink=fault['width']/fault['length'] , pad=0.0)
    B.set_label('Slip [m]')
    C = ax.contour(X, Y, rptm, 5, colors='blue')
    plt.clabel(C, fontsize=12, fmt='%2.1f', inline=1)
    ax.plot(fault['hypo_along_strike']*fault['length'], fault['hypo_down_dip']*fault['width'], '*g', markersize=20)
    plt.xlabel('Distance Along Strike [km]')
    plt.ylabel('Distance Down Dip [km]')
    plt.gca().invert_yaxis()
    plt.gca().set_aspect('equal')
    file_slip = folder + "/Slip.pdf"
    plt.savefig(file_slip)
    return
