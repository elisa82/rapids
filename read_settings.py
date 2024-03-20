def create_settings():
    path_code_hisada = None
    path_code_speed = None
    path_code_ucsb = None
    path_cubit = None
    path_topography = None
    path_code_ucsb_Yoffe = None
    settings = {
        "path_code_hisada": path_code_hisada,
        "path_code_speed": path_code_speed,
        "path_code_ucsb": path_code_ucsb,
        "path_cubit": path_cubit,
        "path_topography": path_topography,
        "path_code_ucsb_Yoffe": path_code_ucsb_Yoffe
    }
    return settings

def read_settings(file_settings):

    settings = create_settings()

    input = {}
    with open(file_settings) as fp:
        line = fp.readline()
        while line:
            if line.strip().find('=') >= 0:
                key, value = line.strip().split('=', 1)
                input[key.strip()] = value.strip()
            line = fp.readline()

    settings['path_code_hisada'] = input['path_code_hisada']
    settings['path_code_speed'] = input['path_code_speed']
    settings['path_code_ucsb'] = input['path_code_ucsb']
    settings['path_cubit'] = input['path_cubit']
    settings['path_topography'] = input['path_topography']
    settings['path_code_ucsb_Yoffe'] = input['path_code_ucsb_Yoffe']


    return settings
