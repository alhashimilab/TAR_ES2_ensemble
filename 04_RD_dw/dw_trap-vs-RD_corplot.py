import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as  sns
import matplotlib
from matplotlib.offsetbox import AnchoredText

from uncertainties import ufloat, umath
import uncertainties.unumpy as unumpy
from uncertainties.umath import *
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['legend.title_fontsize'] = 20
sns.set_context('poster', font_scale=1.5)
matplotlib.rcParams['axes.linewidth']= 4


df = pd.read_excel('/Users/gengainan/Desktop/00_ES2_ensemble/Data/01_RD/04_RD_dw/CS_uucgES2vsRD.xlsx', sheet_name='Sheet1')
uucg = unumpy.uarray((df.uucg, len(df.uucg)*[0]))
RD = unumpy.uarray((df.RD, df.RD_error))

palette ={"Lower helix": 'black',
          "Bulge": sns.color_palette("tab10")[1],
          "Upper helix":sns.color_palette("tab10")[2],
          'Apical loop' : sns.color_palette()[9]}

df['renum'] =df['resi'].str[-2:]
df['renum']= df['renum'].astype(int)
df.reset_index(drop=True,inplace=True)


c1p = df.loc[df.nuclei == "c1p"]
aro = df.loc[df.nuclei == "aro"]
nh = df.loc[df.nuclei == "N1/N3"]
c4p = df.loc[df.nuclei == "C4'"]
h13 = df.loc[df.nuclei == "H1/H3"]

c1pRD = unumpy.uarray((c1p.RD, c1p.RD_error))

x = df.uucg
y = df.RD
R2 = np.corrcoef(x,y)[0,1]
rmsd = np.sqrt(np.mean((x-y)**2))
plt.figure(figsize=(16, 12))
ax = sns.scatterplot(data=df, x=unumpy.nominal_values(RD), y=df.uucg, hue = 'Motif',s=1600, edgecolor="none", alpha=0.75, palette=palette,style="nuclei")
ax.errorbar(unumpy.nominal_values(RD), df.uucg, xerr=unumpy.std_devs(RD), fmt=' ', ecolor='k',elinewidth=2,capsize=12,capthick=2)

ax.get_legend().remove()

#xpoints = ypoints = ax.get_xlim()
xpoints = ypoints = (-5.5,4.0)
#
ax.set_xlim(xpoints)
ax.set_ylim(xpoints)
ax.set_xticks(np.arange(-4,4.0,2), lw=4)
ax.set_yticks(np.arange(-4,4.0,2), lw=4)
ax.set_xticklabels(np.arange(-4,4.0,2),fontsize=48)
ax.set_yticklabels(np.arange(-4,4.0,2),fontsize=48)

ax.plot(xpoints, ypoints, linestyle='--', color='k', lw=4, scalex=False, scaley=False)
ax.set_ylabel(r'$\rm{\Delta\omega_{TAR^{UUCG-ES}}}$ (ppm)',fontsize=58, fontweight="bold")
# ax.set_ylabel('ES RD d$\omega(ppm)$')
ax.set_xlabel(r'$\rm{\Delta\omega_{RD}}$ (ppm)',fontsize=58,fontweight="bold")

line = 'N' +'= {:.0f}'.format(df.shape[0])+ '\n'\
    'R$^2$' + '= {:.2f}'.format(R2) + '\n' \
       + 'RMSD' + '= {:.2f}'.format(rmsd) + " ppm"

anchored_text = AnchoredText(line, loc=4, prop=dict(size=48,ha="right"), frameon= False)
ax.add_artist(anchored_text)

plt.savefig("./dw_CS-RDvsUUCGES2.pdf", bbox_inches="tight", dpi=300)
plt.show()
