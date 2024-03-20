def create_input_ucsb_run(folder, layers, fault, computational_param, sites, path_code, mode_ucsb, green):
    import os
    import subprocess
    from pbs.create_input_ucsb_files import create_input_ffsp
    from pbs.create_input_ucsb_files import create_syn1D
    from pbs.create_input_ucsb_files import create_model
    from pbs.create_input_ucsb_files import create_Green
    from pbs.create_input_ucsb_files import create_stations
    from pbs.fault_slip import plot_slip

    folder = folder + '/UCSB'

    if not os.path.exists(folder):
        os.makedirs(folder)
    os.chdir(folder)
    create_model(folder, layers)
    is_moment = 1
    create_input_ffsp(folder, computational_param, fault, is_moment)
    if fault['IDx'] == 'Yoffe':
        subprocess.call([path_code + '/ffsp_dcf_v2'])
    else:
        subprocess.call([path_code + '/ffsp_v2'])
    plot_slip(folder, fault, is_moment)
    if mode_ucsb == 'full':
        if green == 1:
            create_Green(folder, computational_param, fault, sites)
            os.system('mpirun /Users/elisa/Documents/Programmi/UCSB_broadband/FK/FK_MPI_LF/gfbank_mpi')
            os.system('mv model.green.inf Green_Bank.inf')
        create_stations(folder, sites, fault)
        create_syn1D(folder, computational_param)
        subprocess.call([path_code + '/syn_1d'])
    return