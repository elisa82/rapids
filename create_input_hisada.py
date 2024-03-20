def create_input_hisada_run(folder, layers, fault, computational_param, sites, path_code):
    from rapids.create_input_hisada_files import create_input_grflt12s
    from rapids.create_input_hisada_files import create_script_hisada
    from rapids.create_input_hisada_files import create_answers_grflt12s
    from rapids.create_input_hisada_files import create_answers_grfftspmm
    from rapids.create_input_hisada_files import create_answers_phs3sQ
    import os

    folder = folder + '/HISADA'

    if not os.path.exists(folder):
        os.makedirs(folder)

    create_input_grflt12s(folder, computational_param, layers, fault, sites)
    create_script_hisada(folder, path_code)
    create_answers_grflt12s(folder)
    create_answers_grfftspmm(folder)
    create_answers_phs3sQ(folder)
    return