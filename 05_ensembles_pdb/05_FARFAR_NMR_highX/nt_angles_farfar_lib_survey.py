import os
import sys
import time
import json
import pandas as pd
import numpy as np
import learnna_json as lna_json
from pdblib.base import *
from commontool import read, readchar
import linecache
folder_name = 'json_out' 
outname = '05_FARFAR_NMR_highX'
#outname = folder_name[12:]
json_out = './json_out'
info = './sugar_pucker_info_%s.csv'%outname

pairlist = [['pdbid', 'nt_code', 'nt_resnum',  'alpha', 'beta',  'gamma',
             'delta', 'epsilon', 'zeta', 'chi', 'sugar_class','puckering']]

json_list = []
for file in os.listdir(json_out):
    if file.endswith(".json"):
        json_list.append(file)

json_list.sort(key=lambda x: int(x.split('_')[-1][:-5]))
print(json_list)
correct_uucg_list = []
tic = time.time()  # timer: start
for file in json_list:
    json_f1 = os.path.join(json_out, file)
    najson1 = lna_json.NA_JSON()  #initialize class objects
    with open(json_f1) as json_data:  #read each json file
        data1 = json.load(json_data)

    najson1.set_json(data1)  #pass json file to class pbject
    najson1.read_idx()  #set index from own json file

    nts = najson1.json_file['nts']
    for nt in nts[1:-1]:
        nt_code = nt['nt_code']
        nt_resnum =  nt['nt_resnum']
        alpha = nt['alpha']
        beta = nt['beta']
        gamma = nt['gamma']
        delta = nt['delta']
        epsilon = nt['epsilon']
        zeta = nt['zeta']
        chi = nt['chi']
        sugar_class = nt['sugar_class']
        puckering = nt['puckering']
        pairlist.append([file, nt_code, nt_resnum,  alpha, beta,  gamma,
             delta, epsilon, zeta, chi, sugar_class, puckering])

toc = time.time()  # timer: end
print(">>>>>> Total time used %f >>>>>>" % (toc - tic))

pair_df = pd.DataFrame(pairlist[1:], columns=pairlist[0])
pair_df.to_csv(info, index=False)

