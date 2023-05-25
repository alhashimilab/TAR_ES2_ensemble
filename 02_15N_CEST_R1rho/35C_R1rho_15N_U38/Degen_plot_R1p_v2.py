#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import sys
import csv
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'
sns.set_context('poster', font_scale=1.0)
# matplotlib.rcParams['legend.title_fontsize'] = 'xx-small'
sns.set_style('ticks')

# In[6]:


def read_param(path_home, pathtest, num):
    f = open(path_home + '/' + pathtest + '/2-state_Formatted_LocalFits_{}.csv'.format(num))
    csv_reader = csv.reader(f)
    next(csv_reader)
    col = next(csv_reader)
    df = pd.DataFrame(columns=col)
    i=0
    for row in csv_reader:
        if row[0] == 'Red. Chi-sq':
            row.append(0)
        if len(row) == 3:
            df.loc[i] = row
        i += 1
    f.close()
    pB = float(df[df['%s Pars' % num]=='pB (%)']['Fit Value'].values[0])
    dwB = float(df[df['%s Pars' % num]=='dwB (ppm)']['Fit Value'].values[0])
    kexAB = float(df[df['%s Pars' % num]=='kexAB (s^-1)']['Fit Value'].values[0])
    rchi2 = float(df[df['%s Pars' % num]=='Red. Chi-sq']['Fit Value'].values[0])
    return [pB, dwB, kexAB, rchi2]

path = os.getcwd()
#for file in os.listdir(path):
#    if file.endswith('-mc.csv'):
num = "U38N3_wtTAR_35C"
pb_value = []
dw_value = []
kex_value = []

foldername = [name for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]
kex_foldername = [name for name in foldername if len(name.split('_')) == 3 and name.split('_')[1] == 'kex']
pb_foldername = [name for name in foldername if len(name.split('_')) == 3 and name.split('_')[1] == 'pb']
dw_foldername = [name for name in foldername if len(name.split('_')) == 3 and name.split('_')[1] == 'dw']

print (pb_foldername)
print (num)

test_opt_value = read_param(path, 'test', num)
for i in pb_foldername:
    pb_value.append(read_param(path, i, num))
pb_value.append(test_opt_value)
pb_value = np.array(pb_value)
    
for i in kex_foldername:
    kex_value.append(read_param(path, i, num))
kex_value.append(test_opt_value)
kex_value = np.array(kex_value)

for i in dw_foldername:
    dw_value.append(read_param(path, i, num))
dw_value.append(test_opt_value)
dw_value = np.array(dw_value)

print (pb_value[:,0])

fig, ax = plt.subplots(3, 1, figsize = (8 , 12))

ax[0].scatter(pb_value[:,0], pb_value[:,-1], s=40)
ax[1].scatter(dw_value[:,1], dw_value[:,-1], s=40)
ax[2].scatter(kex_value[:,2], kex_value[:,-1], s=40)

ax[0].scatter(test_opt_value[0], test_opt_value[-1], s=10, color = 'r', label = 'best fit')
ax[1].scatter(test_opt_value[1], test_opt_value[-1], s=10, color = 'r', label = 'best fit')
ax[2].scatter(test_opt_value[2], test_opt_value[-1], s=10, color = 'r', label = 'best fit')

for i in range(3):
    ax[i].set_xlabel(['$p$$_{ES}$ (%)', '$\Delta$$\omega$ (ppm)', '$k$$_{ex}$ (s$^{-1}$)'][i])
    ax[i].set_ylabel('r$\chi$$^2$')

name = sys.argv[1]
ax[0].set_title(name)
plt.tight_layout()
plt.savefig('{}.pdf'.format(name), dpi = 300, transparent = True)




