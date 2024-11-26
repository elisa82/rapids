def create_input_ucsb_run(folder, layers, fault, computational_param, sites, path_code, path_code_ucsb_green_HF, 
        path_code_ucsb_green_LF, calculation_mode, green, band_freq, path_data):
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


    if calculation_mode == '--gf' or calculation_mode == '--run':
        if green == 'nogreen':
            if fault['Mw'] == 4.2:
                GF_precomputed_label = 'M4.2_2024_03_27'
            if fault['Mw'] == 6.4:
                GF_precomputed_label = 'M6.4_1976_05_06'
        if 'LF' in band_freq:
            folder_LF = folder + '/HF' #L'ho chiamato HF ma sarebbe il run con la LF
            if not os.path.exists(folder_LF):
                os.makedirs(folder_LF)
            os.chdir(folder_LF)
            os.system('cp ../model.vel model_lf.vel')
            if green == 'green':
                create_Green(folder_LF, computational_param, fault, sites, 'LF')
                command = 'mpirun -np ' + str(computational_param['nproc_gf']) + ' '  + path_code_ucsb_green_LF + '/gfbank_mpi'
                os.system(command)
                os.system('mv model.green_LF.inf Green_Bank.inf')
            if green == 'nogreen':
                command_cp_GF = 'cp '+path_data+'/GF/model.green_LF_'+GF_precomputed_label+' model.green_LF'
                command_cp_GF_info = 'cp '+path_data+'/GF/Green_Bank_LF.inf_'+GF_precomputed_label+' Green_Bank.inf'
                os.system(command_cp_GF)
                os.system(command_cp_GF_info)
        if 'HF' in band_freq:
            folder_HF = folder + '/HF'
            if not os.path.exists(folder_HF):
                os.makedirs(folder_HF)
            os.chdir(folder_HF)
            os.system('cp ../model.vel target.vel')
            if green == 'green':
                create_Green(folder_HF, computational_param, fault, sites, 'HF')
                create_model_HF(folder_HF, layers)
                command = 'mpirun -np ' + str(computational_param['nproc_gf']) + ' '  + path_code_ucsb_green_HF + '/gfbank_mpi'
                os.system(command)
                os.system('mv model.green_HF.inf Green_Bank.inf')
            if green == 'nogreen':
                command_cp_GF = 'cp DATA/GF/model.green_HF_',GF_precomputed_label,' model.green_HF'
                command_cp_GF_info = 'cp DATA/GF/Green_Bank_HF.inf_',GF_precomputed_label,' Green_Bank.inf'
                os.system(command_cp_GF)
                os.system(command_cp_GF_info)


    if calculation_mode == '--seis' or calculation_mode == '--run':
        os.chdir(folder)
        create_stations(folder, sites, fault)
        create_syn1D(folder, computational_param)
        if 'LF' in band_freq:
            folder_LF = folder + '/HF' #L'ho chiamato HF ma sarebbe il run con la LF
            os.chdir(folder_LF)
            os.system('cp ../stations.xy .')
            os.system('cp ../syn_1d.inp .')
            os.system('cp ../source_model.list .')
            if fault['IDx'] == 'Yoffe-DCF':
                os.system('cp ../Source.bst .')
            else:
                os.system('cp ../Source.001 .')
            command = 'mpirun -np ' + str(computational_param['nproc_seis']) + ' '  + path_code + '/syn_1d'
            os.system(command)

        if 'HF' in band_freq:
            folder_HF = folder + '/HF'
            os.chdir(folder_HF)
            os.system('cp ../stations.xy .')
            os.system('cp ../syn_1d.inp .')
            os.system('cp ../source_model.list .')
            if fault['IDx'] == 'Yoffe-DCF':
                os.system('cp ../Source.bst .')
            else:
                os.system('cp ../Source.001 .')
            command = 'mpirun -np ' + str(computational_param['nproc_seis']) + ' '  + path_code + '/syn_1d'
            os.system(command)
    return
