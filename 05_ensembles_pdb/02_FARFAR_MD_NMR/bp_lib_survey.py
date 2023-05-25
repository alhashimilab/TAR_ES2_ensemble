'''Ainan's script for extract all RNA pdbID --> read json file
Need to specify which bp is the bp of interest
'''
# !/usr/bin/python

import os
import sys
import time
import json
import pandas as pd
import numpy as np
import learnna_json as lna_json
from pdblib.base import *
from commontool import read, readchar

'''specify pool'''
tier = 'FARFAR-MD'
#jsondir = '/home/ainan/tmp/AGs-pair_07152021/three_letter_Rename_new_AGs-pair/json_out'  # 3DNA calculated results
#jsondir = '/home/ainan/tmp/new_G27-A36imino/three_letter_Rename_G27-A36imino/json_out'
jsondir = './json_out'

'''specify base pair'''
#bp_interest = ['X.ADE27','X.GUA36', 'X.GUA28', 'X.ADE35']
oupname = 'bp_surveyLib_%s_bps.csv'%tier

outdir = 'bp_survey_lib'
if not os.path.exists(outdir):
    os.mkdir(outdir)
jsonlist = []  # list to store existing json file names

pdblist = [jsonf.split('.')[0] for jsonf in os.listdir(jsondir) if jsonf.startswith('far')]
pdblist.sort(key=lambda x: (x.split('_')[-1]))
blacklist = []

pairlist = [['pdbid', 'nt1', 'nt1_num',  'nt2', 'nt2_num',  'LW','C1C1_dist',
             'hbonds_num', 'hbonds_desc', 'in_helix', 'in_stem','chi1','chi2','sugar1','sugar2']]

tic = time.time()  # timer: start

for idx, pdbid in enumerate(pdblist):
    print("--- Working on [%s] (%d of %d) ---" % (pdbid, idx + 1, len(pdblist)))

    jsonf = pdbid + '.json'
    ijsonf = os.path.join(jsondir, jsonf)

    if os.path.isfile(ijsonf) is False:
        print("WARNING(): Cannot find Json file: %s" % jsonf)
        sys.exit()
        # continue

    # Step1: Read Json file and check if there is hairpins in the json file
    najson = lna_json.NA_JSON()  # initialize class objects
    with open(ijsonf) as json_data:
        data = json.load(json_data)  # read each json file
    najson.set_json(data)  # pass json file to class object
    najson.read_idx()  # set index from own json file

    if najson.json_file.has_key('pairs') == False:  # there is no bps in the json file
        continue

    # Step2: Loop the pairs in the json file:
    bps = najson.json_file['pairs']
    for i, bp in enumerate(bps):
        if bp['nt1'] not in (['G17',"G18"]):
        #if bp['bp'] == 'G-A':
            basepair = bp['bp']
            nt1 = bp['nt1']
            nt1_num = nt1[-2:]
            nt2 = bp['nt2']
            nt2_num = nt2[-2:]
            LW = bp['LW']
            chi1 = bp['chi1']
            chi2 = bp['chi2']
            sugar1 =bp['pucker1']
            sugar2 = bp['pucker2']
            hbonds_num = bp['hbonds_num']
            hbonds_desc = bp['hbonds_desc']
            C1C1_dist = bp['C1C1_dist']
            in_helix = False
            in_stem = False
            # in helices
            if najson.json_file.has_key('helices') == False:
                in_helix = False
            else:
                for helix in najson.json_file['helices']:
                    for hbp in helix['pairs']:
                        if hbp['nt1'] == nt1 and hbp['nt2'] == nt2:
                            in_helix = True

            # in stems
            if najson.json_file.has_key('stems') == False:
                in_stem = False
            else:
                for stem in najson.json_file['stems']:
                    for sbp in stem['pairs']:
                        if sbp['nt1'] == nt1 and sbp['nt2'] == nt2:
                            in_stem = True


            pairlist.append([pdbid, nt1, nt1_num, nt2, nt2_num, LW,C1C1_dist,
                             hbonds_num, hbonds_desc, in_helix, in_stem,chi1,chi2,sugar1,sugar2])
        else:
            continue

toc = time.time()  # timer: end
print(">>>>>> Total time used %f >>>>>>" % (toc - tic))

pair_df = pd.DataFrame(pairlist[1:], columns=pairlist[0])
pair_df.to_csv(os.path.join(outdir,oupname), index=False)



