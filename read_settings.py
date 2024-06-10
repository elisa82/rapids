def create_settings():
    path_data = None
    path_code_hisada = None
    path_code_speed = None
    path_code_ucsb = None
    path_cubit = None
    path_code_ucsb_Yoffe = None
    path_code_ucsb_green_HF = None
    path_code_ucsb_green_LF = None
    settings = {
        "path_data": path_data,
        "path_code_hisada": path_code_hisada,
        "path_code_speed": path_code_speed,
        "path_code_ucsb": path_code_ucsb,
        "path_cubit": path_cubit,
        "path_code_ucsb_Yoffe": path_code_ucsb_Yoffe,
        "path_code_ucsb_green_HF": path_code_ucsb_green_HF,
        "path_code_ucsb_green_LF": path_code_ucsb_green_LF
    }
    return settings

def read_settings(file_settings, code):

    settings = create_settings()

    input = {}
    with open(file_settings) as fp:
        line = fp.readline()
        while line:
            if line.strip().find('=') >= 0:
                key, value = line.strip().split('=', 1)
                input[key.strip()] = value.strip()
            line = fp.readline()
    
    if 'speed' in code:
        settings['path_data'] = input['path_data']
        settings['path_cubit'] = input['path_cubit']
        settings['path_code_speed'] = input['path_code_speed']

    if 'hisada' in code:
        settings['path_data'] = input['path_data']
        settings['path_code_hisada'] = input['path_code_hisada']

    if 'ucsb' in code:
        settings['path_data'] = input['path_data']
        settings['path_code_ucsb'] = input['path_code_ucsb']
        settings['path_code_ucsb_Yoffe'] = input['path_code_ucsb_Yoffe']
        settings['path_code_ucsb_green_HF'] = input['path_code_ucsb_green_HF']
        settings['path_code_ucsb_green_LF'] = input['path_code_ucsb_green_LF']

    return settings
