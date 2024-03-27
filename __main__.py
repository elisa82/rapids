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

#Di default SPEED è in spostamento e MS-DWS in accelerazione!!!
#Tutto il resto in velocità!!!!

if __name__ == '__main__':

    # %% Initial setup
    try:
        fileini = sys.argv[1]
        code = sys.argv[2]
        calculation_mode = sys.argv[3]
    except IndexError:
        sys.exit('usage:\n'
                 'python3 -m rapids #input_file [code] [mode]')

    #read settings
    settings = read_settings('settings.ini')

    # read input file
    [layers, fault, computational_param, sites, plot_param, topo, folder, cineca] = \
        read_input_data(fileini, code, calculation_mode)
    fault, layers = define_missing_parameters(code, layers, fault, computational_param)
    print(fault)

    if code == 'hisada':
        if calculation_mode == '--input':
            create_input_hisada_run(folder, layers, fault, computational_param, sites, settings['path_code_hisada'])

    if code == 'ucsb':
        if calculation_mode == '--run' or calculation_mode == '--run-nogreen':
            mode_ucsb = 'full'
            if fault['IDx'] == 'Yoffe':
                path_code_ucsb = settings['path_code_ucsb_Yoffe']
            else:
                path_code_ucsb = settings['path_code_ucsb']
            if calculation_mode == '--run':
                green = 1
            else:
                green = 0
            create_input_ucsb_run(folder, layers, fault, computational_param, sites, path_code_ucsb,
                                  mode_ucsb, green)
            post_processing(folder, plot_param, code, sites, fault, computational_param)

    if code == 'speed':
            if calculation_mode == '--input':
                if fault['slip_mode'] == 'Archuleta':
                    mode_ucsb = 'slip'
                    if fault['IDx'] == 'Yoffe':
                        path_code_ucsb = settings['path_code_ucsb_Yoffe']
                    else:
                        path_code_ucsb = settings['path_code_ucsb']
                    if not os.path.exists(folder + '/UCSB'):
                        print(fault['width'], fault['length'])
                        create_input_ucsb_run(folder, layers, fault, computational_param, sites,
                                      path_code_ucsb, mode_ucsb, 0) #0 significa non calcolare le funzioni di Green
                create_input_speed_run(folder, layers, fault, computational_param, sites, settings['path_code_speed'],
                                   topo, settings['path_cubit'], cineca)
            if calculation_mode == '--post':
                speed2ascii(folder, sites)
                post_processing(folder, plot_param, code, sites, fault, computational_param)


    with open(folder + '/sites.obs', 'w') as f:
        f.write("lon lat\n")
        for line in range(len(sites['lon'])):
                f.write("%10.4f %10.4f\n" % (sites['lon'][line], sites['lat'][line]))



    if calculation_mode == '--stitch':

        if 'speed' in code and 'ucsb' in code:
            speed2ascii(folder, sites)

            if fault['IDx'] == 'Yoffe':
                path_code_ucsb = settings['path_code_ucsb_Yoffe']
            else:
                path_code_ucsb = settings['path_code_ucsb']
            num_sm = 1
            stitch(folder, path_code_ucsb, computational_param, fault, sites, num_sm)
        else:
            sys.exit('Both ucsb and speed codes must be specified')

        code = code+'stitched'

        post_processing(folder, plot_param, code, sites, fault, computational_param)


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
