from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy as np
import spotpy
import MONICA_adapter
import re

class spot_setup(object):
    def __init__(self, user_params, exp_maps, host):
        self.user_params = user_params
        self.params = []
        for par in user_params:
            parname = par["name"]
            if re.search(r'\d', par["array"]): #check if par["array"] contains numbers
                parname += "_" + par["array"] #spotpy does not allow two parameters to have the same name
            if "derive_function" not in par.keys(): #spotpy does not care about derived params
                self.params.append(spotpy.parameter.Uniform(parname, par["low"], par["high"], par["stepsize"], par["optguess"], par["minbound"], par["maxbound"]))
        self.monica_model = MONICA_adapter.monica_adapter(exp_maps, host)

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        #the vector comes from spotpy, self.user_params holds the information coming from csv file
        simulation = self.monica_model._run(vector, self.user_params)        
        return simulation

    def evaluation(self):
        #trigger a dummy simulation to populate the expected outcome list
        vector = []
        for par in self.user_params:
            if "optguess" in par.keys():
                vector.append(par["optguess"])
        simulation = self.monica_model._run(vector, self.user_params) 
        return self.monica_model.expected_outcome

    def objectivefunction(self,simulation,evaluation):
        objectivefunction= -spotpy.objectivefunctions.rmse(evaluation,simulation)
        return objectivefunction

    def tests_metadata(self):
        return self.monica_model.exp_ids, self.monica_model.checked_rules
