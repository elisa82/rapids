#!/usr/bin/env python3

import sys
import os
import subprocess
from rapids.create_input_hisada import create_input_hisada_run
from rapids.create_input_speed import create_input_speed_run
from rapids.create_input_ucsb import create_input_ucsb_run
from rapids.read_input_data import read_input_data
from rapids.read_settings import read_settings
from rapids.post_processing import post_processing
from rapids.define_missing_parameters import define_missing_parameters
from rapids.stitch_seismograms import stitch
from rapids.conversions import speed2sac, speed2ascii
import json
from rapids.read_input_data import read_folder
from rapids.compute_magnitudes_nsites_matrix import compute_magnitudes_nsites_matrix

#Di default SPEED è in spostamento e MS-DWS in accelerazione!!!
#Tutto il resto in velocità!!!!

#freq_band_gf = [HF/LF/LFHF/HFLF]
#mesh = [yes/no]
#gf = [yes/no]
#np number of processors for ucsb computations

def load_inputs(fileini):
    import numpy as np
    import json

    folder = read_folder(fileini)
    file_fault_parameters = folder + '/fault_parameters.json'
    file_layers = folder + '/layers.json'
    file_sites = folder + '/sites.json'
    file_plot_param = folder + '/plot_param.json'
    file_computational_param = folder + '/computational_param.json'

    f = open(file_layers)
    layers = json.load(f)
    layers['vp'] = np.asarray(layers['vp'])
    layers['vs'] = np.asarray(layers['vs'])
    layers['qp'] = np.asarray(layers['qp'])
    layers['qs'] = np.asarray(layers['qs'])
    layers['rho'] = np.asarray(layers['rho'])
    layers['thk'] = np.asarray(layers['thk'])
    if layers['fqp'] is not None:
        layers['fqp'] = np.asarray(layers['fqp'])
    if layers['fqs'] is not None:
        layers['fqs'] = np.asarray(layers['fqs'])
    f.close()


    f = open(file_fault_parameters)
    fault = json.load(f)
    f.close()

    f = open(file_computational_param)
    computational_param = json.load(f)
    computational_param['optiout'] = np.asarray(computational_param['optiout'])
    computational_param['seed'] = np.asarray(computational_param['seed'])
    computational_param['mlst'] = np.asarray(computational_param['mlst'])
    f.close()

    f = open(file_sites)
    sites = json.load(f)
    sites['X'] = np.asarray(sites['X'])
    sites['Y'] = np.asarray(sites['Y'])
    sites['Z'] = np.asarray(sites['Z'])
    sites['ID'] = np.asarray(sites['ID'])
    sites['lon'] = np.asarray(sites['lon'])
    sites['lat'] = np.asarray(sites['lat'])
    f.close()

    f = open(file_plot_param)
    plot_param = json.load(f)
    f.close()

    return layers, fault, computational_param, sites, plot_param, folder

if __name__ == '__main__':

    # %% Initial setup

    try:
        fileini = sys.argv[1]
        code = sys.argv[2]
        calculation_mode = sys.argv[3]
    except IndexError:
        sys.exit('usage:\n'
            'python3 -m rapids #input_file [code] [mode]')

    if fileini == 'none' and code == 'none' and calculation_mode == 'table':
        start = 4.0
        stop = 10.0
        step = 0.1
        compute_magnitudes_nsites_matrix(start,stop,step)
        sys.exit("Created the table magnitudes-nsites for UrgentShake")

        
    #read settings
    settings = read_settings('settings.ini', code)
    path_data = settings['path_data']
    #sistemare la lettura dei settings

    if calculation_mode == '--input' or calculation_mode == '--run' or calculation_mode == '--source': 
        # read input file
        [layers, fault, computational_param, sites, plot_param, topo, folder] = \
            read_input_data(fileini, code)

        if not os.path.exists(folder):
            os.makedirs(folder)
        fault, layers = define_missing_parameters(code, layers, fault, computational_param, path_data, topo, folder, sites)

        file_fault_parameters = folder + '/fault_parameters.json'
        file_layers = folder + '/layers.json'
        file_sites = folder + '/sites.json'
        file_plot_param = folder + '/plot_param.json'
        file_computational_param = folder + '/computational_param.json'

        if not os.path.exists(file_fault_parameters):
            with open(file_fault_parameters, 'w') as f:
                json.dump(fault, f, ensure_ascii=False)

        if not os.path.exists(file_layers):
            layers['vp'] = list(layers['vp'])
            layers['vs'] = list(layers['vs'])
            layers['qp'] = list(layers['qp'])
            layers['qs'] = list(layers['qs'])
            layers['rho'] = list(layers['rho'])
            layers['thk'] = list(layers['thk'])
            if layers['fqp'] is not None:
                layers['fqp'] = list(layers['fqp'])
            if layers['fqs'] is not None:
                layers['fqs'] = list(layers['fqs'])
            with open(file_layers, 'w') as f:
                json.dump(layers, f, ensure_ascii=False)

        if not os.path.exists(file_sites):
            sites['X'] = list(sites['X'])
            sites['Y'] = list(sites['Y'])
            sites['Z'] = list(sites['Z'])
            sites['ID'] = list(sites['ID'])
            sites['lon'] = list(sites['lon'])
            sites['lat'] = list(sites['lat'])
            with open(file_sites, 'w') as f:
                json.dump(sites, f, ensure_ascii=False)

                
        if not os.path.exists(file_plot_param):
            with open(file_plot_param, 'w') as f:
                json.dump(plot_param, f, ensure_ascii=False)

        if not os.path.exists(file_computational_param):
            computational_param['optiout'] = list(computational_param['optiout'])
            computational_param['seed'] = list(computational_param['seed'])
            if computational_param['mlst'] is not None:
                computational_param['mlst'] = list(computational_param['mlst'])
            else:
                computational_param['mlst'] = []
            with open(file_computational_param, 'w') as f:
                json.dump(computational_param, f, ensure_ascii=False)
            
    layers, fault, computational_param, sites, plot_param, folder = load_inputs(fileini)

    if 'ucsb' in code:
        path_code_ucsb_green_HF = settings['path_code_ucsb_green_HF']
        path_code_ucsb_green_LF = settings['path_code_ucsb_green_LF']
        if fault['IDx'] == 'Yoffe-DCF':
            path_code_ucsb = settings['path_code_ucsb_Yoffe']
        else:
            path_code_ucsb = settings['path_code_ucsb']
        if computational_param['gf'] == 'yes':
            green = 'green'
        else:
            green = computational_param['gf']
        freq_band = computational_param['freq_band_gf']

    if code == 'hisada':
        if calculation_mode == '--input':
            create_input_hisada_run(folder, layers, fault, computational_param, sites, settings['path_code_hisada'], path_data)

    if code == 'ucsb':
        if calculation_mode != '--post':
            create_input_ucsb_run(folder, layers, fault, computational_param, sites, path_code_ucsb,
                                path_code_ucsb_green_HF, path_code_ucsb_green_LF, calculation_mode, green, freq_band, 
                                path_data)

        if calculation_mode == '--run':
            if 'LF' in computational_param['freq_band_gf'] and 'HF' in computational_param['freq_band_gf']:
                num_sm = 1
                stitch(folder, path_code_ucsb, computational_param, fault, sites, num_sm, code)
                
                try:
                    from mpi4py import MPI
                    comm  = MPI.COMM_WORLD
                    rank  = comm.Get_rank()
                    nranks =comm.size

                except:
                    rank   = 0
                    nranks = 1
                    comm = None

                post_processing(folder, plot_param, 'stitched-U', sites, fault, computational_param, path_data, comm)

        if calculation_mode == '--post':
            #questo va rivisto x le due condizioni HF e LF che ora non sono contemplate
            try:
                from mpi4py import MPI
                comm  = MPI.COMM_WORLD
                rank  = comm.Get_rank()
                nranks =comm.size

            except:
                rank   = 0
                nranks = 1
                comm = None

            post_processing(folder, plot_param, 'ucsb', sites, fault, computational_param, path_data, comm)

    if 'ucsb' in code and calculation_mode == '--source':
            create_input_ucsb_run(folder, layers, fault, computational_param, sites, path_code_ucsb,
                                path_code_ucsb_green_HF, path_code_ucsb_green_LF, calculation_mode, green, freq_band, 
                                path_data)


    if code == 'speed':
        if calculation_mode == '--input':
            if fault['slip_mode'] == 'Archuleta':
                settings = read_settings('settings.ini', 'ucsb')
                if fault['IDx'] == 'Yoffe-DCF':
                    path_code_ucsb = settings['path_code_ucsb_Yoffe']
                else:
                    path_code_ucsb = settings['path_code_ucsb']
                if not os.path.exists(folder + '/UCSB'):
                    path_code_ucsb_green_HF = settings['path_code_ucsb_green_HF']
                    path_code_ucsb_green_LF = settings['path_code_ucsb_green_LF']
                    create_input_ucsb_run(folder, layers, fault, computational_param, sites,
                                      path_code_ucsb, path_code_ucsb_green_HF, path_code_ucsb_green_LF, '--source', '', 
                                      [], path_data)
            settings = read_settings('settings.ini', 'speed')
            path_cubit = settings['path_cubit']
            create_input_speed_run(folder, layers, fault, computational_param, sites, settings['path_code_speed'],
                                   topo, path_cubit, path_data)

        if calculation_mode == '--post':
            speed2ascii(folder, sites)

            try:
                from mpi4py import MPI
                comm  = MPI.COMM_WORLD
                rank  = comm.Get_rank()
                nranks =comm.size

            except:
                rank   = 0
                nranks = 1
                comm = None


            post_processing(folder, plot_param, code, sites, fault, computational_param, path_data, comm)


    #with open(folder + '/sites.obs', 'w') as f:
    #    f.write("lon lat\n")
    #    for line in range(len(sites['lon'])):
    #            f.write("%10.4f %10.4f\n" % (sites['lon'][line], sites['lat'][line]))

    if calculation_mode == '--stitch':

        if 'speed' in code:
            speed2ascii(folder, sites)

        if 'ucsb' in code:

            if fault['IDx'] == 'Yoffe-DCF':
                path_code_ucsb = settings['path_code_ucsb_Yoffe']
            else:
                path_code_ucsb = settings['path_code_ucsb']
        
        num_sm = 1
        stitch(folder, path_code_ucsb, computational_param, fault, sites, num_sm, code)

        #code = 'stitched-'+code
        #post_processing(folder, plot_param, code, sites, fault, computational_param, path_data)

    if calculation_mode == '--post' and code != 'ucsb' and code != 'speed':
        try:
            from mpi4py import MPI
            comm  = MPI.COMM_WORLD
            rank  = comm.Get_rank()
            nranks =comm.size

        except:
            rank   = 0
            nranks = 1
            comm = None

        settings = read_settings('settings.ini', code)
        post_processing(folder, plot_param, code, sites, fault, computational_param, path_data, comm)


# 	HYPO_DOWN_DIP=position_along_width*width
# 	ZTOR=depth_hypc-HYPO_DOWN_DIP*np.sin(np.radians(dip))
# #	if(ZTOR<min_depth):
# #		ZTOR=min_depth
# #		depth_hypc=ZTOR+HYPO_DOWN_DIP*np.sin(np.radians(dip))
# #		depth_max_fault=ZTOR+width*np.sin(np.radians(dip))
# #		if(depth_max_fault>max_depth):
# #			print('Error',depth_max_fault,max_depth,ZTOR,width,mag,depth_hypc)
# #

#     slip_srcmod = np.zeros((NumSubFault_RW, NumSubFault_RLD))
#     rupt_srcmod = np.zeros((NumSubFault_RW, NumSubFault_RLD))
#     rise_srcmod = np.zeros((NumSubFault_RW, NumSubFault_RLD))
#     rake_srcmod = np.zeros(NumSubFault_RW, NumSubFault_RLD))
#
#     inc = 0
#     for i in range(NumSubFault_RW):
#         m = inc + 1
#         for k in range(NumSubFault_RLD):
#             slip_srcmod[j, k] = Slip[m]
#             rupt_srcmod[j, k] = RupTime[m]
#             rise_srcmod[j, k] = RiseTime[m]
#             rake_srcmod[j, k] = RakeAng[m]
#             m = m + NumSubFault_RW
#         inc = inc + 1
#
# for id_slip = counter: counter + counter_end - 1
#
# if id_slip < 10
#     file_name = [Output_FileName, '.00', num2str(id_slip)];
# elseif
# id_slip < 100
# file_name = [Output_FileName, '.0', num2str(id_slip)];
# else
# file_name = [Output_FileName, '.', num2str(id_slip)];
#
# end
# % file_name = 'tmp.txt.001';
# fid = fopen(file_name, 'r');
# Attr = str2num(fgetl(fid));
# Num_lines = Attr(2);
# Slip_fun = Attr(3);
# TpTr = Attr(4);
#
# for j = 1: Num_lines
# Attr = str2num(fgetl(fid));
# North(j) = Attr(1);
# East(j) = Attr(2);
# Down(j) = Attr(3);
# Slip(j) = Attr(4);
# RupTime(j) = Attr(5);
# RiseTime(j) = Attr(6);
# dummy = Attr(7);
# StrikeAng(j) = Attr(8);
# DipAng(j) = Attr(9);
# RakeAng(j) = Attr(10);
# end
#
# fclose(fid);
