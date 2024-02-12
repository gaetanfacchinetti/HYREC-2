import copy

from .wrapperhyrec import init_INPUT_COSMOPARAMS


class HyRecParamsDefault:

    def __init__(self,  new_params):

        self.create_params(new_params)
        
        for key, value in self._params.items():
            self.__dict__[key] = value


    # define the params dictionary
    def create_params(self, new_params):
        self._params = self._defaults.copy()

        if new_params is None:
            return

        for key, value in new_params.items():
            # can only set parameters already present in the default dictionnary
            if key in self._defaults:
                self._params[key] = value
            else:
                raise ValueError("Trying to initialise parameter " + key + " that is not in the default list")

    # modify the __setattr__ to prevent any unwanted modifications
    def __setattr__(self, name, value):

        self.__dict__[name] = value

        if name in self._defaults:
            raise Exception("Attributes are read only! Use update() to modify them.")

    # reinitialise to the default values
    def set_defaults(self):

         for key, value in self._defaults.items():
            if key in self._defaults:
                self._params[key] = value
                self.__dict__[key] = value


    def update(self, **new_params):
        """ update parameters """

        for key, value in new_params.items():
            if key in self._defaults:
                self._params[key] = value
                self.__dict__[key] = value
            else:
                raise ValueError("Trying to modify a parameter not in default")
    
    def __str__(self):
        return "HyRec object with parameters: " + str(self._params)
    
    def __call__(self):
        return self._params
    
    


class HyRecCosmoParams(HyRecParamsDefault):

    def __init__(self, new_params = None) -> None:
        
        self._defaults = {
        "h" : 6.735837e-01, 
        "T0" : 2.7255,
        "Omega_b" : 0.0494142797907188,
        "Omega_cb" : 0.31242079216478097,
        "Omega_k" : 0.0,
        "w0" : -1.0,
        "wa" : 0.0,
        "Neff" : 3.046,
        "Nmnu" : 1,
        "mnu1" : 0.06,
        "mnu2" : 0.0 ,
        "mnu3" : 0.0, 
        "YHe"  : 0.245,
        "fsR"  : 1.0,
        "meR"  : 1.0,}

        super().__init__(new_params)

        self.check()

    def check(self) -> None: 
        
        if (self.mnu3 > self.mnu2) or (self.mnu3 > self.mnu1) or (self.mnu2 > self.mnu1):
            raise ValueError("The neutrinos masses should be ordered from highest to lowest")

        return None 


class HyRecInjectionParams(HyRecParamsDefault):

    def __init__(self, new_params = None) -> None:
        
        self._defaults = {
        "pann" : 0.0, 
        "pann_halo" : 0.0,
        "ann_z" : 0.0,
        "ann_zmax" : 0.0,
        "ann_zmin" : 0.0,
        "ann_var" : 0.0,
        "ann_z_halo" : 0.0,
        "decay" : 0.0,
        "on_the_spot" : 0.0,
        "Mpbh" : 1.0,
        "fpbh" : 0.0,}
        
        super().__init__(new_params)