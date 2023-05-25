#!/usr/bin/env python2.7

import pandas as pd
import os
import sys
reload(sys)
sys.setdefaultencoding('utf8')


expt = '47_wtTAR_30C_sw_MM_1.4mM'
title = "wtTAR 30$^{\circ}$C (sw)"
peaks = ['U38', 'G28', 'U42','G18','G36', 'G21', 'G43', 'G26']
outdier = './'
path = './'
name = expt
        
outname = expt+'_' + 'kex.csv'
inpath =[x for x in os.listdir(path) if x.endswith('_fitparm.csv')][0]
df = pd.read_csv(inpath)
df.iloc[[0],:]= df.iloc[[0],:].round(4)
df.iloc[1:,:]= df.iloc[1:,:].round(2)

df.astype(str)
dfo = pd.DataFrame()
dfo['para'] = df.iloc[:,0]

for i in range(len(peaks)):
    posi = i+1
    df.iloc[:, posi * 2 - 1] = df.iloc[:, posi * 2 - 1].astype(str)
    df.iloc[:, posi * 2 ] = df.iloc[:, posi * 2 ].astype(str)

    dfo[peaks[i]] = df.iloc[:,posi*2-1].str.cat(df.iloc[:,posi*2], sep=u"\u00B1")
    dfo.iloc[[0,1,3,5,7,9,10],:].to_csv(os.path.join(outdier, outname))
