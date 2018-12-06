import json
import sys
import monica_io
import zmq
import csv
import os
from datetime import date
import collections
import threading
from threading import Thread
from collections import defaultdict
from copy import deepcopy


class monica_adapter(object):
    def __init__(self, exp_maps, host):

        self.host = host

        #for multi-experiment: create a M-2 relationship between exp_IDs and param files
        self.IDs_paramspaths = {}
        for exp_map in exp_maps:
            self.IDs_paramspaths[exp_map["exp_ID"]] = {}
            self.IDs_paramspaths[exp_map["exp_ID"]]["species"] = exp_map["species_file"]
            self.IDs_paramspaths[exp_map["exp_ID"]]["cultivar"] = exp_map["cultivar_file"]

        #load rules (for MONICA output)
        with open("calibration_rules.json") as _:
            self.calibration_rules = json.load(_)

        #data structures for spotpy
        self.expected_outcome = [] #former "observations"; stores success scores
        self.simulation_outcome = [] #former "simulations"; stores outcome after checking rules
        
        self.exp_ids = [] #mirrors spotpy data structure to provide the user with info on the exp id (success/failure)
        self.checked_rules = [] #mirrors spotpy data structure to provide the user with info on the checked rule (success/failure)

        self.species_params={} #map to store different species params sets avoiding repetition
        self.cultivar_params={} #map to store different cultivar params sets avoiding repetition

        #create envs
        self.envs = []
        for exp_map in exp_maps:
            with open(exp_map["sim_file"]) as simfile:
                sim = json.load(simfile)
                sim["crop.json"] = exp_map["crop_file"]
                sim["site.json"] = exp_map["site_file"]
                #sim["climate.csv"] = exp_map["climate_file"]

            with open(exp_map["site_file"]) as sitefile:
                site = json.load(sitefile)

            with open(exp_map["crop_file"]) as cropfile:
                crop = json.load(cropfile)
                mycrop = exp_map["crop_ID"]
                crop["crops"][mycrop]["cropParams"]["species"][1] = exp_map["species_file"]
                crop["crops"][mycrop]["cropParams"]["cultivar"][1] = exp_map["cultivar_file"]

            env = monica_io.create_env_json_from_json_config({
                "crop": crop,
                "site": site,
                "sim": sim,
                "climate": ""
            })
            
            #customize monica output
            env["events"] = []
            for rule_id, rule_specs in self.calibration_rules.iteritems():
                env["events"].append(rule_specs["custom event"][0]) #customize output
                env["events"].append(rule_specs["custom event"][1])
            
            #climate is read by the server
            env["csvViaHeaderOptions"] = sim["climate.csv-options"]
            env["csvViaHeaderOptions"]["start-date"] = sim["climate.csv-options"]["start-date"]
            env["csvViaHeaderOptions"]["end-date"] = sim["climate.csv-options"]["end-date"]
            env["pathToClimateCSV"] = []
            env["pathToClimateCSV"].append(exp_map["climate_file"])

            #monica_io.add_climate_data_to_env(env, sim) this might not be supported anymore

            for ws in env["cropRotation"][0]["worksteps"]:
                if ws["type"] == "Seed" or ws["type"] == "Sowing":
                    self.species_params[exp_map["species_file"]] = ws["crop"]["cropParams"]["species"]
                    self.cultivar_params[exp_map["cultivar_file"]] = ws["crop"]["cropParams"]["cultivar"]
                    break
            
            env["customId"] = exp_map["exp_ID"]
            self.envs.append(env)

        self.context = zmq.Context()
        self.socket_producer = self.context.socket(zmq.PUSH)
        self.socket_producer.connect("tcp://" + self.host["server"] + ":" + self.host["push-port"])

    def run(self,args):
        return self._run(*args)

    def _run(self,vector, user_params):

        #clean data structures for spotpy and info to the user
        self.expected_outcome = [] #former "observations"; stores success scores
        self.simulation_outcome = [] #former "simulations"; stores outcome after checking rules
        self.exp_ids = []
        self.checked_rules = []

        def seek_set_param(par, p_value, model_params):
            p_name = par["name"]
            array = par["array"]
            add_index = False
            if isinstance(model_params[p_name], int) or isinstance(model_params[p_name], float):
                add_index = False
            elif len(model_params[p_name]) > 1 and isinstance(model_params[p_name][1], basestring):
                add_index = True #the param contains text (e.g., units)
            if array.upper() == "FALSE":
                if add_index:
                    model_params[p_name][0] = p_value
                else:
                    model_params[p_name] = p_value
            else: #param is in an array (possibly nested)
                array = array.split("_") #nested array
                if add_index:
                    array = [0] + array
                if len(array) == 1:
                    model_params[p_name][int(array[0])] = p_value
                elif len(array) == 2:
                    model_params[p_name][int(array[0])][int(array[1])] = p_value
                elif len(array) == 3:
                    model_params[p_name][int(array[0])][int(array[1])][int(array[2])] = p_value
                else:
                    print "param array too nested, contact developers"
            

        #set params according to spotpy sampling. Update all the species/cultivar available
        for i in range(len(user_params)):                        #loop on the user params
            for s in self.species_params:               #loop on the species
                if user_params[i]["name"] in self.species_params[s]:
                    seek_set_param(user_params[i],
                    user_params[i]["derive_function"](vector, self.species_params[s]) if "derive_function" in user_params[i] else vector[i],
                    self.species_params[s])
                else:
                    break                                   #break loop on species if the param is not there
            for cv in self.cultivar_params:                 #loop on the cultivars
                if user_params[i]["name"] in self.cultivar_params[cv]:
                    seek_set_param(user_params[i],
                    user_params[i]["derive_function"](vector, self.cultivar_params[cv]) if "derive_function" in user_params[i] else vector[i],
                    self.cultivar_params[cv])
                else:
                    break
        
        #launch parallel thread for the collector
        collector = Thread(target=self.collect_results)
        collector.daemon = True
        collector.start()

        #send jobs to the MONICA server
        for env in self.envs:
            species = self.species_params[self.IDs_paramspaths[env["customId"]]["species"]]
            cultivar = self.cultivar_params[self.IDs_paramspaths[env["customId"]]["cultivar"]]
            for ws in env["cropRotation"][0]["worksteps"]:
                if ws["type"] == "Seed" or ws["type"] == "Sowing":
                    ws["crop"]["cropParams"]["species"] = species
                    ws["crop"]["cropParams"]["cultivar"] = cultivar
                    break

            self.socket_producer.send_json(env)

        #wait until the collector finishes
        collector.join()
                 
        return self.simulation_outcome

        
    def collect_results(self):

        socket_collector = self.context.socket(zmq.PULL)
        socket_collector.connect("tcp://" + self.host["server"] + ":" + self.host["pull-port"])
        
        received_results = 0
        leave = False
        while not leave:
            try:
                rec_msg = socket_collector.recv_json()
            except:
                continue            
            
            results_rec = []
            for res in rec_msg["data"]:
                try:
                    rule_name = res["outputIds"][0]["displayName"]
                    event_out = res["results"][0]
                    rule_math = self.calibration_rules[rule_name]["rule math"]
                    func_obj = lambda x: eval(rule_math)

                    success_score = self.calibration_rules[rule_name]["success score"]
                    failure_score = self.calibration_rules[rule_name]["failure score"]
                    
                    for event_val in event_out:
                        #info for the user
                        self.exp_ids.append(rec_msg["customId"])
                        self.checked_rules.append(rule_name)
                        #append success score to expected outcome
                        self.expected_outcome.append(success_score)

                        #append result from rule check (either success or failure score)                        
                        if func_obj(event_val):
                            self.simulation_outcome.append(success_score)
                        else:
                            self.simulation_outcome.append(failure_score) 

                    if len(self.simulation_outcome) != len(self.expected_outcome):
                        print("len simulation_outcome != len expected_outcome while adding res from exp id " + rec_msg["customId"])                      
                    
                except:
                    print("missing data in exp id " + rec_msg["customId"])
                    exit()

            received_results += 1
            #print("total received: " + str(received_results))

            if received_results == len(self.envs):                
                leave = True
