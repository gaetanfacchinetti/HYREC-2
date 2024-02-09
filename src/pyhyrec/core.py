# UNTESTED CODE
class HyRecParamsDefault:

    def __init__(self,  new_params):
        self._default = None

        for key, value in self._default.iteritems():
            setattr(self, key, new_params.get(key, value))

    

class HyRecCosmoParams(HyRecParamsDefault):

    def __init__(self, new_params):
        self._defaults = {"h" : 0.6776}

        super().__init__(new_params)

