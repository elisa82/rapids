def create_input_speed_run(folder, layers, fault, computational_param, sites, path_code, topo, path_cubit, cineca, path_data):
    import os
    from rapids.create_input_speed_files import create_input_mate
    from rapids.create_input_speed_files import create_input_LS
    from rapids.create_input_speed_files import create_input_speed_file
    from rapids.create_input_speed_files import create_script_speed
    from rapids.create_input_speed_files import create_input_extended
    from rapids.create_input_speed_files import redefine_layers
    from rapids.create_input_speed_files import create_mesh
    from rapids.conversions import convert_mesh

    folder_ucsb = folder + '/UCSB'
    folder = folder + '/SPEED'

    if not os.path.exists(folder):
        os.makedirs(folder)
    os.chdir(folder)
    create_input_LS(folder, sites)

    layers = redefine_layers(layers, fault, computational_param)

    if fault['fault_type'] == 'extended':
        create_input_extended(folder, fault, layers, folder_ucsb)

    create_input_mate(folder, computational_param, layers, fault)

    cubit_journal, file_exodus_mesh, file_exodus_topo = create_mesh(folder, computational_param, layers, fault, sites,
                                                                    topo, path_cubit, path_data)
    command_cubit_jou = path_cubit + ' -nographics python3 ' + cubit_journal
    os.system(command_cubit_jou)

    fileNameExodus = file_exodus_mesh
    (FilePathMesh, ExtMesh) = os.path.splitext(fileNameExodus)
    fileNameTxtMesh = FilePathMesh + '.txt'
    command_ncdump_mesh = 'ncdump ' + file_exodus_mesh + ' > ' + fileNameTxtMesh
    os.system(command_ncdump_mesh)
    convert_mesh(fileNameTxtMesh, 'mesh', folder)

    if topo == 1:
        fileNameExodus = file_exodus_topo
        (FilePathTopo, ExtTopo) = os.path.splitext(fileNameExodus)
        fileNameTxtTopo = FilePathTopo + '.txt'
        command_ncdump_topo = 'ncdump ' + file_exodus_topo + ' > ' + fileNameTxtTopo
        os.system(command_ncdump_topo)
        convert_mesh(fileNameTxtTopo, 'xyz', folder)

    file_mesh_final = FilePathMesh.split('/')[len(FilePathMesh.split('/'))-1]
    create_input_speed_file(folder, computational_param, file_mesh_final)

    create_script_speed(folder, path_code, sites, cineca)
    return
