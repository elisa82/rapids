def create_input_ucsb_run(folder, layers, fault, computational_param, sites, path_code, path_code_ucsb_green_HF, 
        path_code_ucsb_green_LF,calculation_mode, green, band_freq, path_data):
    import os
    import subprocess
    from rapids.create_input_ucsb_files import create_input_ffsp
    from rapids.create_input_ucsb_files import create_syn1D
    from rapids.create_input_ucsb_files import create_model
    from rapids.create_input_ucsb_files import create_model_HF
    from rapids.create_input_ucsb_files import create_Green
    from rapids.create_input_ucsb_files import create_stations
    from rapids.fault_slip import plot_slip

    folder = folder + '/UCSB'

    if not os.path.exists(folder):
        os.makedirs(folder)
    os.chdir(folder)
    if calculation_mode == '--source' or calculation_mode == '--run':
        create_model(folder, layers)
        is_moment = 1
        create_input_ffsp(folder, computational_param, fault, is_moment)
        if fault['IDx'] == 'Yoffe-DCF':
            subprocess.call([path_code + '/ffsp_dcf_v2'])
        else:
            subprocess.call([path_code + '/ffsp_v2'])
        plot_slip(folder, fault, is_moment)
    if calculation_mode == '--seis' or calculation_mode == '--run':
        create_stations(folder, sites, fault)
        create_syn1D(folder, computational_param)
        if 'LF' in band_freq:
            folder_LF = folder + '/HF' #L'ho chiamato HF ma sarebbe il run con la LF
            if not os.path.exists(folder_LF):
                os.makedirs(folder_LF)
            os.chdir(folder_LF)
            os.system('cp ../stations.xy .')
            os.system('cp ../model.vel model_lf.vel')
            if green == 'green':
                create_Green(folder_LF, computational_param, fault, sites, 'LF')
                command = 'mpirun ' + path_code_ucsb_green_LF + '/gfbank_mpi'
                os.system(command)
                os.system('mv model.green_LF.inf Green_Bank.inf')
            os.system('cp ../syn_1d.inp .')
            os.system('cp ../source_model.list .')
            if fault['IDx'] == 'Yoffe-DCF':
                os.system('cp ../Source.bst .')
            else:
                os.system('cp ../Source.001 .')
            subprocess.call([path_code + '/syn_1d'])
        if 'HF' in band_freq:
            folder_HF = folder + '/HF'
            if not os.path.exists(folder_HF):
                os.makedirs(folder_HF)
            os.chdir(folder_HF)
            os.system('cp ../stations.xy .')
            os.system('cp ../model.vel target.vel')
            if green == 'green':
                create_Green(folder_HF, computational_param, fault, sites, 'HF')
                create_model_HF(folder_HF, layers)
                command = 'mpirun -np 1 ' + path_code_ucsb_green_HF + '/gfbank_mpi'
                os.system(command)
                os.system('mv model.green_HF.inf Green_Bank.inf')
            os.system('cp ../syn_1d.inp .')
            os.system('cp ../source_model.list .')
            if fault['IDx'] == 'Yoffe-DCF':
                os.system('cp ../Source.bst .')
            else:
                os.system('cp ../Source.001 .')
            subprocess.call([path_code + '/syn_1d'])
    return
