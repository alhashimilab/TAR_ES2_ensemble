import sys
import os
import pandas as pd
peaks = ['U38', 'G28', 'U42','G18', 'G36', 'G21', 'G43','G26']
inpath =  sys.argv[1]
outdir = 'formated' + inpath
df = pd.read_csv(inpath)
df.iloc[[1],:]= df.iloc[[1],:].round(4)
df.iloc[1:,:]= df.iloc[1:,:].round(2)
df.astype(str)
dfo = pd.DataFrame()
dfo['para'] = df.iloc[:,0]
for i in range(len(peaks)):
    posi = i+1
    df.iloc[[1], posi * 2 - 1] = df.iloc[[1], posi * 2 - 1].astype(float).round(4)*100
    df.iloc[[1], posi * 2 - 1] = df.iloc[[1], posi * 2 - 1].astype(str)
    df.iloc[[1], posi * 2 ] = df.iloc[[1], posi * 2 ].astype(float).round(4)*100
    df.iloc[[1], posi * 2 ] = df.iloc[[1], posi * 2 ].astype(str)
    
    df.iloc[[3,5,8,9,10], posi * 2 - 1] = df.iloc[[3,5,8,9,10], posi * 2 - 1].astype(float).round(2).astype(str)
    df.iloc[[3,5,8,9,10], posi * 2 ] = df.iloc[[3,5,8,9,10], posi * 2 ].astype(float).round(2).astype(str)
    dfo[peaks[i]] = df.iloc[:,posi*2-1].str.cat(df.iloc[:,posi*2], sep=u"\u00B1")
dfo.iloc[[1],:] = dfo.iloc[[1],:] + str(" %")
dfo.iloc[[0,1,3,5,8,9,10],:].to_csv(outdir)
print(dfo)
