class parameters(object):

    quenchNL_file           = "./data/Quench5.root"
    Ncer_file               = "./data/Ncer_local.root"

    #electron parameters :
    Ysct = 1400.
    kB = 6.5e-3
    kC = 1.0

    a = 0.98
    b = 0.044
    n = 1.4

    # gamma samples :
    nGammaSamples = 5000

    ### detector config :
    detConfig = "Det1"

    Y_det1 = 3134.078 / 2.223
    Ge68mean_det1 = 1.30880e3
    Ge68sigma_det1 = 38.0672
    


    ###### Getter Property #####
    @staticmethod
    def get_quenchNL_file():
        return parameters.quenchNL_file


    @staticmethod
    def get_Ncer_file():
        return parameters.Ncer_file

    @staticmethod
    def get_Ngamma_samples():
        return parameters.nGammaSamples

    @staticmethod
    def get_Y():
        if parameters.detConfig == "Det1":
            return parameters.Y_det1

    @staticmethod
    def get_Ge68mean():
        if parameters.detConfig == "Det1":
            return parameters.Ge68mean_det1


    @staticmethod
    def get_Ge68sigma():
        if parameters.detConfig == "Det1":
            return parameters.Ge68sigma_det1


    @staticmethod
    def get_Ysct():
        return parameters.Ysct

    @staticmethod
    def get_kB():
        return parameters.kB

    @staticmethod
    def get_kC():
        return parameters.kC

    @staticmethod
    def get_res_a():
        return parameters.a


    @staticmethod
    def get_res_b():
        return parameters.b

    @staticmethod
    def get_res_n():
        return parameters.n

    @staticmethod
    def get_gamma_spectrum(name):
        filename = "/junofs/users/miaoyu/simulation/LS_Sim/jobs/" + parameters.detConfig + "/gamma/" + name + "_new.root"
        return filename

    @staticmethod
    def get_gamma_samples(name):
        filename = "/junofs/users/miaoyu/energy_model/fitter/energyModel_Fit/new_fitter/data/gamma/" + name + "_J19.root"
        return filename

    ###### Setter Property #####
    @staticmethod
    def set_Y(val):
        parameters.Y = val

    @staticmethod
    def set_Ysct(val):
        parameters.Ysct = val

    @staticmethod
    def set_kB(val):
        parameters.kB = val

    @staticmethod
    def set_kC(val):
        parameters.kC = val

    @staticmethod
    def set_res_a(val):
        parameters.a = val

    @staticmethod
    def set_res_b(val):
        parameters.b = val

    @staticmethod
    def set_res_n(val):
        parameters.n = val
