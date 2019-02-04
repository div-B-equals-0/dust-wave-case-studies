# ---
# jupyter:
#   jupytext_format_version: '1.2'
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
#   language_info:
#     codemirror_mode:
#       name: ipython
#       version: 3
#     file_extension: .py
#     mimetype: text/x-python
#     name: python
#     nbconvert_exporter: python
#     pygments_lexer: ipython3
#     version: 3.6.3
# ---

# ## Proper motions of stars near AE Aur

# +
import numpy as np
from astropy.table import Table 

datafile = '../data/AE-Aur-Gaia/1549211065246O-result.vot'

tab = Table.read(datafile)
# -

tab

import pandas as pd
from matplotlib import pyplot as plt
# %matplotlib inline
import seaborn as sns
sns.set_color_codes()

df = tab.to_pandas()
df

ddf = df[((df.pmdec - df.pmdec.mean()) / df.pmdec.std()).abs() < 3]


ddf.describe()

# +
from scipy import stats

def corrfunc(x, y, **kws):
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(x, y)
#    r, p = stats.pearsonr(x, y)
    ax = plt.gca()
    if pvalue < 0.05:
        fontcolor = 'k' if pvalue < 0.003 else 'k'
        fontalpha = 1.0 if pvalue < 0.003 else 0.5
        ax.annotate(f"$r = {rvalue:.3f}$",
                    xy=(.05, .9), xycoords=ax.transAxes, 
                    fontsize='small', color=fontcolor, alpha=fontalpha)
        ax.annotate(rf"$m = {slope:.4f} \pm {stderr:.4f}$",
                    xy=(.05, .05), ha='left', xycoords=ax.transAxes, 
                    fontsize='small', color=fontcolor, alpha=fontalpha)


# -

columns = ['ra', 'dec', 'parallax', 'pmra', 'pmdec']
g = sns.PairGrid(ddf, 
                 vars=columns, 
                )
g = g.map_diag(plt.hist, edgecolor="w")
g = g.map_upper(plt.scatter, alpha=0.5, edgecolor='w')
g = g.map_upper(corrfunc)
g = g.map_lower(sns.kdeplot, cmap="Blues_d")

datafile = '../data/AE-Aur-Gaia/1549214713051O-result.vot'
df2 = Table.read(datafile).to_pandas()
ddf2 = df2[((df2.pmdec - df2.pmdec.mean()) / df2.pmdec.std()).abs() < 3]
ddf2.describe()


g = sns.PairGrid(ddf2, 
                 vars=columns, 
                )
g = g.map_diag(plt.hist, edgecolor="w")
g = g.map_upper(plt.scatter, alpha=0.5, edgecolor='w')
g = g.map_upper(corrfunc)
g = g.map_lower(sns.kdeplot, cmap="Blues_d")

# So, looks like there might be a clump centered on RA = 79.6 or so that could merit further investigation.  
#
# But the average PM is about -10 in dec and about +3 in RA


