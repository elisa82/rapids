def plot_slip(folder, fault, is_moment):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib

    #plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 16
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'

    fid = open('ffsp.inp', 'r+')
    fid.readline()
    line = fid.readline()
    c = line.split()
    lx = float(c[0])
    ly = float(c[1])
    line = fid.readline()
    c = line.split()
    cxp = float(c[0])
    cyp = float(c[1])
    hypo_depth = float(c[2])
    line = fid.readline()
    c = line.split()
    xref = float(c[0])
    yref = float(c[1])
    line = fid.readline()
    c = line.split()
    Mo = float(c[0])
    if fault['IDx'] == 'Yoffe-DCF':
        fc1 = float(c[1])
        fc2 = float(c[2])
        rv = float(c[3])
    else:
        fc = float(c[1])
        rv = float(c[2])
    line = fid.readline()
    c = line.split()
    pt_rt = float(c[0])
    line = fid.readline()
    c = line.split()
    strike = float(c[0])
    dip = float(c[1])
    rake = float(c[2])
    fid.readline()
    line = fid.readline()
    c = line.split()
    nx = int(c[0])
    ny = int(c[1])
    line = fid.readline()
    c = line.split()
    top = int(c[0])
    right = int(c[1])
    bottom = int(c[2])
    left = int(c[3])
    fid.close()
    dx = lx / (nx - 1)
    dy = ly / (ny - 1)
    A = dx * dy * 1e6
    x = (np.arange(0, nx) * dx) - lx / 2
    y = (np.arange(0, ny) * dy)

    up = np.array([left * dx - lx / 2,
                   top * dy])  # ,[lx/2-right*dx,top*dy],[lx/2-right*dx,ly-bottom*dy],[left*dx-lx/2,ly-bottom*dy]])
    down = np.array([[left * dx - lx / 2, ly - bottom * dy], [lx / 2 - right * dx, ly - bottom * dy]])

    for num in range(1, 2):
        index = '{0:03}'.format(num)
        sequence = str(index)
        if fault['IDx'] == 'Yoffe-DCF':
            file = "Source.bst"
            prefix_file_slip = "Slip.bst"
        else:
            file = "Source." + sequence
            prefix_file_slip = "Slip." + sequence
        data = np.loadtxt(file, skiprows=1)
        if is_moment == 1:
            moment = np.transpose(data[:, 3].reshape(nx, ny))  # moment in Nm
        if is_moment == 3:
            slip = np.transpose(data[:, 3].reshape(nx, ny))  # moment in Nm
        rptm = np.transpose(data[:, 4].reshape(nx, ny))  # rupture time in s
        rstm = np.transpose(data[:, 5].reshape(nx, ny))  # rise time in s

        # print np.mean(slip),np.mean(rstm)
        X, Y = np.meshgrid(x, y)

        if is_moment == 1:
            data2 = np.loadtxt('model.vel', skiprows=1)
            b = data2[:, 1] * 1000  # in m/s
            rho = data2[:, 2] * 1000
            h = np.cumsum(data2[:, 3]) * 1000
            mu = rho * b * b

            # Material Property at each point
            lay_num = np.zeros((len(data[:, 3])))
            z_glo = data[:, 2]
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
                slip[ipt] = data[:, 3][ipt] / (mu[ipt] * A)

            slip = np.transpose(slip.reshape(nx, ny))

        # plotting
        f, ax = plt.subplots(1, 1)
        A = ax.imshow(slip[::-1], cmap='gist_heat_r',
                      extent=(-lx / 2 - dx / 2, lx / 2 + dx / 2, -dy / 2, ly + dy / 2),
                      interpolation='nearest')
        B = plt.colorbar(A, orientation='vertical', shrink=ly / lx, pad=0.0)
        B.set_label('Slip [m]')
        C = ax.contour(X, Y, rptm, 5, colors='blue')
        plt.clabel(C, fontsize=12, fmt='%2.1f', inline=1)
        ax.plot(cxp - lx / 2, cyp, '*g', markersize=20)
        ax.plot([left * dx - lx / 2 - dx / 2, lx / 2 - right * dx + dx / 2], [top * dy - dy / 2, top * dy - dy / 2],
                'g')
        ax.plot([left * dx - lx / 2 - dx / 2, lx / 2 - right * dx + dx / 2],
                [ly - bottom * dy + dy / 2, ly - bottom * dy + dy / 2], 'g')
        ax.plot([left * dx - lx / 2 - dx / 2, lx / 2 - right * dx + dx / 2], [top * dy - dy / 2, top * dy - dy / 2],
                'g')
        ax.plot([left * dx - lx / 2 - dx / 2, left * dx - lx / 2 - dx / 2],
                [ly - bottom * dy + dy / 2, top * dy - dy / 2], 'g')
        ax.plot([lx / 2 - right * dx + dx / 2, lx / 2 - right * dx + dx / 2],
                [ly - bottom * dy + dy / 2, top * dy - dy / 2], 'g')

        plt.xlabel('Distance Along Strike [km]')
        plt.ylabel('Distance Down Dip [km]')
        plt.gca().invert_yaxis()
        plt.gca().set_aspect('equal')
        file_slip = folder + "/" + prefix_file_slip + ".pdf"
        plt.savefig(file_slip)
        # plt.show()
        return
