def number_grid_sites(receiver_grid_step_km, receiver_maximum_dist):
    import math

    nstep_grid = math.floor(receiver_maximum_dist / receiver_grid_step_km) + 1
    count = 0
    for i in range(2 * nstep_grid + 1):
        for j in range(2 * nstep_grid + 1):
            count += 1
    return count

def compute_magnitudes_nsites_matrix(start,stop,step):
    import numpy as np
    from rapids.read_input_data import compute_maximum_distance
    receiver_grid_step_km = 2.5

    magnitudes = np.arange(start, stop, step)
    ngrid = np.zeros((1,len(magnitudes)))
    for i in range(len(magnitudes)):
        receiver_maximum_dist = compute_maximum_distance(magnitudes[i])
        ngrid[0,i] = number_grid_sites(receiver_grid_step_km, receiver_maximum_dist)

    data = np.column_stack((magnitudes.T, ngrid.T))
    np.savetxt('magnitudes_nsites.txt', data, fmt='%0.1f', header='Magnitude Nsites')




