def write_stitch_VMname(folder, sites):
    fid = open(folder + '/VMname.list', 'w')
    for i in range(len(sites['ID'])):
        fid.write('{}\n'.format('model.vel'))
    fid.close()
    return


def write_stitch_input(folderLF, folderHF, folder, computational_param, fault, num_sm):
    fid = open(folder + '/stitch.inp', 'w')
    fid.write('{}\n'.format('stations.xy')) #Enter the name of file listing the stations
    fid.write('{}\n'.format('VMname.list')) #Enter the name of file listng 1D VM for each stations
    fid.write('{}\n'.format(folderHF+'/')) #Enter the name of folder with SPEED synthetics
    fid.write('{}\n'.format(folderLF+'/'))  # Enter the name of folder with UCSB synthetics
    fid.write('{} {}\n'.format(computational_param['freq_join'], computational_param['fmax_ucsb']))  # Enter the joint frequency and fmax
    fid.write('{}\n'.format(fault['hypo']['Z']))  # Enter the depth of hypocenter (km)
    fid.write('{}\n'.format(num_sm))  # Enter the number of source model
    #if computational_param['output_type'] == 'dis':
    #    id_motion = 1
    #if computational_param['output_type'] == 'vel':
    #    id_motion = 2
    #if computational_param['output_type'] == 'acc':
    #    id_motion = 3
    id_motion = 2 #deciso che è in vel
    fid.write('{}\n'.format(id_motion))  # Output Displacement (=1), Vel. (2), or Acc (3)
    fid.close()
    return


def stitch(folder, path_code, computational_param, fault, sites, num_sm, code):
    import subprocess
    import os

    if 'speed' in code:
        folderLF = folder + '/SPEED'
    if 'ucsb' in code:
        folderHF = folder + '/UCSB'
    elif 'msdwn' in code:
        folderHF = folder + '/MS-DWN'

    if not os.path.exists(folder):
        os.makedirs(folder)

    os.chdir(folderHF)
    subprocess.call([path_code + '/conv3Comp'])

    os.chdir(folder)
    copy_station_list = 'cp ' + folderHF + '/stations.xy .'
    os.system(copy_station_list)
    copy_station_list = 'cp ' + folderHF + '/stations.ll .'
    os.system(copy_station_list)
    copy_vel_model = 'cp ' + folderHF + '/model.vel .'
    os.system(copy_vel_model)
    write_stitch_VMname(folder, sites)
    write_stitch_input(folderLF, folderHF, folder, computational_param, fault, num_sm)
    subprocess.call([path_code + '/stitch'])
    return


