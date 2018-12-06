from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import os
import spotpy
import spotpy_setup_MONICA
import csv
from datetime import date
import numpy as np
import json
from collections import defaultdict

config = {
    "host": {
        "server": "localhost",
        "push-port": "6666",
        "pull-port": "7777",
    },    
    "rep": 20,
    "cal-method": "SCE-UA"
}

def make_lambda(excel):
    return lambda v, p: eval(excel)

crop_sim_site_MAP = "crop_sim_site_MAP.csv"

#read general settings
exp_maps = []
all_exps = []
basepath = os.path.dirname(os.path.abspath(__file__))
with open(crop_sim_site_MAP) as exp_mapfile:
    dialect = csv.Sniffer().sniff(exp_mapfile.read(), delimiters=';,\t')
    exp_mapfile.seek(0)
    reader = csv.reader(exp_mapfile, dialect)
    next(reader, None)  # skip the header
    for row in reader:
        all_exps.append(row[0])
        exp_map = {}
        exp_map["exp_ID"] = row[0]
        exp_map["sim_file"] = basepath+"\\sim_files\\"+row[1]
        exp_map["crop_file"] = basepath+"\\crop_files\\"+row[2]
        exp_map["site_file"] = basepath+"\\site_files\\"+row[3]
        exp_map["climate_file"] = basepath+"\\climate_files\\"+row[4]
        exp_map["species_file"] = basepath+"\\param_files\\"+row[5]
        exp_map["cultivar_file"] = basepath+"\\param_files\\"+row[6]
        exp_map["crop_ID"] = row[7]
        exp_maps.append(exp_map)

#read params to be calibrated
params = []
with open('calibratethese.csv') as paramscsv:
    dialect = csv.Sniffer().sniff(paramscsv.read(), delimiters=';,\t')
    paramscsv.seek(0)
    reader = csv.reader(paramscsv, dialect)
    next(reader, None)  # skip the header
    for row in reader:
        p={}
        if len(row) == 9 and row[8] != "":
            p["name"] = row[0]
            p["derive_function"] = make_lambda(row[8])
            p["array"] = row[1]
        else:
            p["name"] = row[0]
            p["array"] = row[1]
            p["low"] = float(row[2])
            p["high"] = float(row[3])
            p["stepsize"] = float(row[4])
            p["optguess"] = float(row[5])
            p["minbound"] = float(row[6])
            p["maxbound"] = float(row[7])
        
        params.append(p)

spot_setup = spotpy_setup_MONICA.spot_setup(params, exp_maps, config["host"])
results = []


if config["cal-method"] == "SCE-UA":
    sampler = spotpy.algorithms.sceua(spot_setup, dbname='SCEUA', dbformat='ram')
    sampler.sample(repetitions=config["rep"], ngs=len(params)+1, kstop=10, pcento=10, peps=10)
elif config["cal-method"] == "MLE":
    sampler = spotpy.algorithms.mle(spot_setup, dbname='MLE_CMF', dbformat='ram')
    sampler.sample(repetitions=config["rep"])

results.append(sampler.getdata())

best = sampler.status.params

with open('optimizedparams.csv', 'wb') as _:
    writer = csv.writer(_)        
    for i in range(len(best)):
        outrow=[]
        arr_pos = ""
        if params[i]["array"].upper() != "FALSE":
            arr_pos = params[i]["array"]        
        outrow.append(params[i]["name"]+arr_pos)
        outrow.append(best[i])
        writer.writerow(outrow)
    if len(params) > len(best):
        reminder = []
        reminder.append("Don't forget to calculate and set derived params!")
        writer.writerow(reminder)
    print('optimized parameters saved!')

#run sim with optimized params
best_out = spot_setup.simulation(best)
expected_out = spot_setup.evaluation()
exp_ids, checked_rules = spot_setup.tests_metadata()

with open('summary_rules_bestrun.csv', 'wb') as _:
    writer = csv.writer(_)
    header = ["exp_id", "rule", "expected", "simulated"]
    writer.writerow(header)
    for i in range(len(best_out)):
        outrow =[]
        outrow.append(exp_ids[i])
        outrow.append(checked_rules[i])
        outrow.append(expected_out[i])
        outrow.append(best_out[i])
        writer.writerow(outrow)
    print('summary best run saved!')

print("finished!")



