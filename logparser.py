import pandas as pd
from enum import Enum
import os
import csv
import shutil
import numpy as np
import math 
import matplotlib.pyplot as plt





def extract_scip(entry_, file_path):
    ls = open(file_path).readlines()
    #print(ls)
    #print(file)
    stat_keys = ["Total Time",  "Dual Bound", "Primal Bound",  "Gap", "nodes"]
    stat_dict = {}
    entry = entry_
    for l in ls:
        for stat_key in stat_keys:
            if  l.split(":")[0].strip() == stat_key:
                stat_dict[stat_key] = l
    #print(stat_dict["Dual Bound"].split()[2])
    entry["total_time"] = float(stat_dict["Total Time"].split()[3])
    #print(stat_dict["Total Time"].split()[3])
    #print(stat_dict["intermis"].split())
    entry["primal"] = -float(stat_dict["Primal Bound"].split()[3])
    entry["dual"] = -float("NAN") if stat_dict["Dual Bound"].split()[3] == "-" else -float(stat_dict["Dual Bound"].split()[3])
    entry["gap"] = float(stat_dict["Gap"].split()[2])
    #print(stat_dict["nodes"].split())
    entry["nodes"] = int(stat_dict["nodes"].split()[2])
    return entry





def Stat(name):
    return {"setting": name, "solved": 0, "gap": 0, "dual": 0.0, "primal": 0, "nodes": 0} 


#print(opt_dict)


log_path = os.getcwd() + "/logs"
logs = os.listdir(log_path)

entries = []


for log in logs:
    log_ = log
    log = log[0:-4]
    #print(log.split("_"))
    setting = log.split("_")[-1]
    if setting == "normal":
        instance = log[0:-6]
    else:
        instance = log[0:-6]
    insclass = instance.split("_")[0]
    entry={}
    entry["instance"] = instance
    entry["setting"] = setting
    entry["class"] =  insclass
    entry = extract_scip(entry, log_path + "/"+log_)
    entries.append(entry)


    

def Stat(sname, pname):
    return {"setting": sname, "pclass": pname, "total": 0, "nodes": 0, "gap": 0,  "dual": 0, "primal": 0, "nodes_lst": [], "gap_lst": [], "dual_lst": [], "primal_lst": []} 

display_keys = ["setting", "pclass", "dual", "primal", "gap", "nodes"]

numerical_keys = ["dual", "primal", "gap", "nodes"]

details = ""

settings = ["scip1", "scip2", "scip3" , "scip4" , "scip5" , "scip6"]
pclasses = ['block2', 'normal']

classstats = {}




for setting in settings:
    for pclass in pclasses:    
        classstats[(setting, pclass)] = []
        #print(pclass,  subpclass)
        for entry in entries:
            if entry["class"] == pclass and entry["setting"] == setting:
                #print("1")
                instance = entry["instance"]
                classstats[(setting, pclass)].append(entry)
        #print(classstats[(setting, pclass)])

def add(stat, entry):
    #print(stat, entry, entry["closed"] is float("nan"))
    #print(entry)
    #stat["solved"] += entry["issolved"] 
    stat["total"] += 1
    stat["nodes_lst"].append(entry["nodes"])
    stat["gap_lst"].append(entry["gap"]) 
    stat["dual_lst"].append(entry["dual"])
    stat["primal_lst"].append(entry["primal"])


def SGM(lst, total, bias):
    return np.exp(np.sum([np.log(ele + bias) for ele in lst ]) / total) - bias


def avgStat(stat):
    stat["nodes"] =  SGM(stat["nodes_lst"], stat["total"], 1)
    stat["primal"] = SGM(stat["primal_lst"], stat["total"], 1)
    stat["dual"] = SGM(stat["dual_lst"] , stat["total"], 1)
    stat["gap"] = SGM(stat["gap_lst"] , stat["total"], 1)


def printStat(stat):
    #print(setting, pclass, stat)
    print([(display_key, round(stat[display_key], 2) if display_key in numerical_keys else stat[display_key] ) for display_key in display_keys])


stats = {}
allstat = {}


for setting in settings:
    allstat[setting] = Stat(setting, "all")
    for pclass in pclasses:
            stats[(pclass,setting)] = Stat(setting, pclass)
            for entry in classstats[(setting, pclass)]:
                add(stats[(pclass,setting)], entry)   
                add(allstat[setting], entry)  
            #print("\n", pclass, stats[(pclass,setting)])           
            avgStat(stats[(pclass,setting)])
            #stats[(pclass, setting)]["relative"] = (stats[(pclass, setting)]["closed"] + 1e-9) /  (stats[(pclass, "default")]["closed"] + 1e-9)
            printStat(stats[(pclass, setting)])  
    avgStat(allstat[setting])
    #allstat[setting]["relative"] = (allstat[setting]["closed"] + 1e-9) /  (allstat["default"]["closed"] + 1e-9)       
    printStat(allstat[setting])

