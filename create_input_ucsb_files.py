import numpy as np

def create_input_ffsp(folder, computational_param, fault, is_moment):
    import math
    import sys

    file_ffsp = folder + '/ffsp.inp'
    fid = open(file_ffsp, 'w')

    if fault['IDx'] == 'exp':
        id_sf_type = 1
    if fault['IDx'] == 'Archuleta':
        id_sf_type = 8
    if fault['IDx'] == 'Yoffe-DCF':
        id_sf_type = 8

    if fault['IDx'] == 'Yoffe-DCF':
        fid.write('{} {} {} {} \n'.format(id_sf_type, computational_param['fmin_ucsb'], computational_param['fmax_ucsb'],
                                      ' ! index of slip rate function'))
    else:
        fid.write('{} {} {} \n'.format(id_sf_type, computational_param['fmax_ucsb'], ' ! index of slip rate function'))
    fid.write('{} {} {}\n'.format(fault['length'], fault['width'], 'length and width of main fault'))

    HYPO_DOWN_DIP = fault['hypo_down_dip'] * fault['width']  # The reference is top center of fault
    cyp = HYPO_DOWN_DIP
    HYPO_ALONG_STK = fault['hypo_along_strike'] * fault['length'] - fault[
        'length'] / 2.  # The reference is top center of fault
    cxp = fault[
              'length'] / 2. + HYPO_ALONG_STK  # è come se cxp fosse semplicemente fault['hypo_along_strike']*fault['length']
    fid.write('{} {} {} {}\n'.format(cxp, cyp, fault['hypo']['Z'],
                                     ' ! distance in the strike and down dip direction, and the depth of the hypocenter'))
    # top center of fault is in the origin
    # 	deg2rad=np.pi/180.
    # controllare queste relazioni. non so se c'è un errore perchè in faultGlobal sembrerebbe che al posto di HYPO_ALONG_STK e HYPO_DOWN_DIP c'è cxp e xyp.
    # Da verificare. Per ora uso xref_hypc = 0 e yref_hypc = 0 come in 3ptool
    # xref_hypc = np.cos(strike*deg2rad)*HYPO_ALONG_STK - np.cos(dip*deg2rad)*np.sin(strike*deg2rad)*HYPO_DOWN_DIP
    # yref_hypc = np.sin(strike*deg2rad)*HYPO_ALONG_STK + np.cos(dip*deg2rad)*np.cos(strike*deg2rad)*HYPO_DOWN_DIP
    # change center of fault in the hypocentre
    xref_hypc = 0.
    yref_hypc = 0.
    fid.write('{} {} {}\n'.format(xref_hypc, yref_hypc, ' ! xref_hypc,yref_hypc'))

    # if M <= 6.5
    #     Corner_Freq = 0.5;
    # elseif
    # M <= 7
    # Corner_Freq = 0.4;
    # elseif
    # M <= 7.5
    # Corner_Freq = 0.08;
    # end

    if fault['IDx'] == 'Yoffe-DCF':
        fc_main_1 = 1 / (np.pi * fault['rupture_duration'])
        fc_main_2 = 0.8 / fault['rise_time']
        # JA19_2S (Non-Self-Similar)
        #if fault['Mw'] <= 3.3:
        #    sys.exit('JA19_2S relationship not defined for M<3.3')
        #elif 3.3 < fault['Mw'] <= 5.3:
        #    fc_main_1 = 10 ** (1.474 - 0.415 * fault['Mw'])
        #elif 5.3 < fault['Mw'] < 7.3:
        #    fc_main_1 = 10 ** (2.375 - 0.585 * fault['Mw'])
        #else:
        #    sys.exit('JA19_2S relationship not defined for M>=7.3')
        # JA19_2S (Self-Similar)
        # fc_main_1 = 10 ** (1.754 - 0.5 * fault['Mw'])
        # fc_main_2 is the same for both models
        #fc_main_2 = 10 ** (3.250 - 0.5 * fault['Mw'])

        fid.write('{} {} {} {} {}\n'.format(fault['Mo'], fc_main_1, fc_main_2,
                                            fault['rupture_velocity'], '! Moment_o,fc_main_1,fc_main_2,rv_avg'))
    else:
        fc = 10 ** (2.502 - 0.5 * fault['Mw'])
        fid.write('{} {} {} {}\n'.format(fault['Mo'], fc, fault['rupture_velocity'],
                                      '! Moment_o,fc,rv_avg'))
    fid.write('{} {}\n'.format(fault['Tp_Tr'], ' ! ratio of time between slip-rate increase and rise time'))
    fid.write('{} {} {}\n'.format(fault['strike'], fault['dip'], fault['rake']))
    Dip_rand = 10
    Rake_rand = 10
    fid.write('{} {} {}\n'.format(Dip_rand, Rake_rand, '! max perturbations of the dip and rake'))

    fid.write('{} {} {}\n'.format(fault['number_subfaults_strike'], fault['number_subfaults_dip'],
                                  '!subfaults number along strike and dip'))
    if fault['number_subfaults_strike'] > 50:
        SubFaults_Right_taper = 5
        SubFaults_Left_taper = 5
    else:
        SubFaults_Right_taper = 1
        SubFaults_Left_taper = 1
    if fault['number_subfaults_dip'] > 50:
        SubFaults_Top_taper = 5
        SubFaults_Bottom_taper = 5
    else:
        SubFaults_Top_taper = 1
        SubFaults_Bottom_taper = 1

    fid.write('{} {} {} {} {}\n'.format(SubFaults_Top_taper, SubFaults_Right_taper, SubFaults_Bottom_taper,
                                        SubFaults_Left_taper, '! numbers of subfaults to be tapered at each side'))
    fid.write('{} {} {} {}\n'.format(computational_param['seed'][0], computational_param['seed'][1],
                                     computational_param['seed'][2], '! seed for generating random numbers'))
    max_num_realizations = computational_param['realizations']

    fid.write('{} {} {}\n'.format(1, max_num_realizations, '! first and last index of random generations'))
    fid.write('{}\n'.format('model.vel'))
    North_Angle = 0.0
    fid.write('{} {}\n'.format(North_Angle, '! angle north to x'))
    fid.write('{} {}\n'.format(is_moment, ' ! moment (1), slip-area (2), slip(3)'))
    fid.write('{}\n'.format('Source'))
    fid.close()


def create_syn1D(folder, computational_param):
    file_syn = folder + '/syn_1d.inp'
    fid = open(file_syn, 'w')
    fid.write('1 1                    ! # of point source for each subfault\n')
    fid.write('60.0,    40.0,    20.0   ! #Perturbation on strike, rake, and dip\n')
    fid.write('0.1,  3.0                  ! #two frequency for perturbation\n')
    duration = computational_param['dt_ucsb'] * computational_param['npts_ucsb']
    fid.write('{} {} {}\n'.format(computational_param['kappa'], duration, ' ! # duration of outputing GM'))
    fid.write('source_model.list\n')
    station_name = 'stations.xy'
    fid.write('{}\n'.format(station_name))
#    if computational_param['output_type'] == 'dis':
#        type_output_int = 1
#    if computational_param['output_type'] == 'vel':
    type_output_int = 2 # Ho messo di default ucsb in velocità
#    if computational_param['output_type'] == 'acc':
#        type_output_int = 3
    fid.write('{} {}\n'.format(type_output_int,'       ! #1 for Displacement, 2 for Vel., 3 for Acc'))
    fid.write('2       ! #1 for SAC;  2 for TXT;  3 for Binary\n')
    fid.close()


def create_model(folder, layers):
    fid = open(folder + '/model.vel', 'w')
    nlayers = len(layers['thk'])
    fid.write('{}    {}\n'.format(nlayers, 1.0))
    thk_ucsb = list(layers['thk'])
    thk_ucsb[nlayers - 1] = 0
    for i in range(len(layers['thk'])):
        fid.write('{:7.2f}{:7.2f}{:7.2f}{:7.1f}{:8.1f}{:8.1f}\n'.format(layers['vp'][i], layers['vs'][i],
                                                                        layers['rho'][i], thk_ucsb[i], layers['qp'][i],
                                                                        layers['qs'][i]))
    fid.close()


def create_model_HF(folder, layers):
    fid = open(folder + '/model_hf.vel', 'w')
    nlayers = len(layers['thk'])
    thk_ucsb = list(layers['thk'])
    thk_ucsb[nlayers - 1] = 0

    vp_moho = 7.8

    nlayers = 0
    thick_crust = 0
    index_moho = -999
    #vp_crust = 0
    #vs_crust = 0
    #rho_crust = 0
    #qp_crust = 0
    #qs_crust = 0
    for i in range(len(layers['thk'])):
        if layers['vp'][i] < 7.8:
            thick_crust += layers['thk'][i] 
            #vp_crust += layers['vp'][i]*layers['thk'][i]
            #vs_crust += layers['vs'][i]*layers['thk'][i] 
            #rho_crust += layers['rho'][i]*layers['thk'][i] 
            #qp_crust += layers['qp'][i]*layers['thk'][i] 
            #qs_crust += layers['qs'][i]*layers['thk'][i] 
        else:
            nlayers += 1
            if index_moho < 0:
                index_moho = i
    #vp_crust = vp_crust/thick_crust 
    #vs_crust = vs_crust/thick_crust
    #rho_crust = rho_crust/thick_crust
    #qp_crust = qp_crust/thick_crust
    #qs_crust = qs_crust/thick_crust
    nlayers = nlayers + 1
    fid.write('{}    {}\n'.format(nlayers, 1.0))

    #fid.write('{:7.2f}{:7.2f}{:7.2f}{:7.1f}{:8.1f}{:8.1f}\n'.format(vp_crust, vs_crust,
    #                                                            rho_crust, thick_crust, qp_crust,qs_crust))
    fid.write('{:7.2f}{:7.2f}{:7.2f}{:7.1f}{:8.1f}{:8.1f}\n'.format(layers['vp'][index_moho - 1], layers['vs'][index_moho - 1],
                                                            layers['rho'][index_moho - 1], thick_crust, layers['qp'][index_moho-1],
                                                             layers['qs'][index_moho - 1]))
    for i in range(nlayers-1):
            fid.write('{:7.2f}{:7.2f}{:7.2f}{:7.1f}{:8.1f}{:8.1f}\n'.format(layers['vp'][index_moho + i], layers['vs'][index_moho + i],
                                                                    layers['rho'][index_moho + i], thk_ucsb[index_moho + i], layers['qp'][index_moho + i],
                                                                    layers['qs'][index_moho + i]))
    fid.close()


def create_Green(folder, computational_param, fault, sites, type_green):
    import numpy as np
    import math

    fid = open(folder + '/Green.in', 'w')
    if type_green == 'LF':
        file_model = 'model_lf.vel'
    if type_green == 'HF':
        file_model = 'model_hf.vel'
    fid.write('{}\n'.format(file_model))

    dep_step = 0.5

    dep_max = fault['max_fault_depth'] + dep_step * 2
    dep_min = fault['Ztor'] - 1
    if dep_min < 0:
        dep_min = 0.0005

    fid.write('{} {} {}\n'.format(dep_max, dep_min, dep_step))
    fault['number_subfaults_strike'] = int(2 ** math.ceil(math.log(fault['length'] / fault['subfault_length'] + 1, 2)))
    fault['number_subfaults_dip'] = int(2 ** math.ceil(math.log(fault['width'] / fault['subfault_width'] + 1, 2)))
    dist = []
    for k in range(len(sites['ID'])):
        statX = sites['X'][k] - fault['hypo_utm']['X']
        statY = sites['Y'][k] - fault['hypo_utm']['Y']
        dx = fault['length'] / fault['number_subfaults_strike']
        dy = fault['width'] / fault['number_subfaults_dip']
        for i in range(fault['number_subfaults_strike']):
            posx = (dx * i + dx / 2)
            #origin in the hypocentre
            xcoord = (posx - fault['hypo_along_strike'] * fault['length'])*1000
            for j in range(fault['number_subfaults_dip']):
                posy = (dy * j + dy / 2)
                ycoord = (posy - fault['hypo_down_dip'] * fault['width'])*1000
                deg2rad = np.pi / 180.
                xcoord_rot = np.cos(fault['strike'] * deg2rad) * xcoord - np.cos(fault['dip'] * deg2rad) * \
                             np.sin(fault['strike'] * deg2rad) * ycoord
                ycoord_rot = np.sin(fault['strike'] * deg2rad) * xcoord + np.cos(fault['dip'] * deg2rad) * \
                             np.cos(fault['strike'] * deg2rad) * ycoord
                d = np.sqrt((statX - xcoord_rot)**2 + (statY - ycoord_rot)**2)
                dist.append(d)
    dist_max = (np.max(dist)+1000) / 1000
    dist_min = (np.min(dist)-1000) / 1000
    if dist_min < 0:
        dist_min = 0.0005
    d_step = 2
    # N.B. dist max and dist min are only approssimative, nn tengono conto dell'effettiva orientazione della faglia e posizione ipocentro sul piano di faglia
    fid.write('{} {} {}\n'.format(dist_max, dist_min, d_step))
    t_cor = 3.0
    # number of time steps, time increment, seconds to be saved before first arrival. This should never be set to 0 (because of wrap‐ around artifacts!!!
    fid.write('{} {} {}\n'.format(computational_param['npts_ucsb'], computational_param['dt_ucsb'], t_cor))
    if type_green == 'HF':
        fid.write('{}\n'.format('model.green_HF'))  # "The name of file to store Green Bank"
        fid.write('{} {}\n'.format(computational_param['qzero'], computational_param['alpha']))
    if type_green == 'LF':
        fid.write('{}\n'.format('model.green_LF'))  # "The name of file to store Green Bank"
    fid.write('{}\n'.format('0'))
    fid.write('{}\n'.format(''))
    fid.close()


def create_stations(folder, sites, fault):
    import sys
    # x = north
    # y = east
    # z = down
    # check if coordinates system is ok!!!!
    nobs = len(sites['Z'])
    station_name = folder + '/stations.xy'
    fid = open(station_name, 'w')
    epiX = 0.0
    epiY = 0.0
    placeholder = 0.0
    fid.write('{} {} {} {}\n'.format(nobs, epiX, epiY, placeholder))
    for k in range(nobs):
        fid.write('{}\n'.format(sites['ID'][k]))
        statX = sites['X'][k] - fault['hypo_utm']['X']
        statY = sites['Y'][k] - fault['hypo_utm']['Y']
        statZ = sites['Z'][k]
        if statZ != 0.0:
            sys.exit('Error: depth of station should be 0')
        comp1 = 0.000000
        comp2 = 90.000000
        fid.write('{} {} {} {} {}\n'.format(statX, statY, statZ, comp1, comp2))
    fid.close()

    fid = open(folder + '/stations.ll', 'w')
    nobs = len(sites['Z'])
    fid.write('{}\n'.format(nobs))
    for k in range(nobs):
        fid.write('{} {} {}\n'.format(sites['lon'][k], sites['lat'][k], sites['ID'][k]))
    fid.close()
    return

def create_script_ucsb_source(folder, path_code_hisada):
    fid = open(folder +'/script_ucsb_source.sh', 'w')
    fid.write('#!/bin/bash\n')
    #fid.write('"'+path_code_ucsb+'/phs3sQ_slip_vr.out" < ph.dat\n')
    #fid.write('"'+path_code_ucsb+'/grflt12s_slip_vr_subf.out" < gr.dat\n')
    #fid.write('"'+path_code_ucsb+'/grfftspmm.out" < grspm.dat\n')
    #fid.close()
    return
