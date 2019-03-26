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

# # Gaia distances to Cygnus X

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from astropy.table import Table

# %matplotlib inline
sns.set_color_codes()

tab = Table.read("../data/Cygnus/Gaia/1553609645651O-result.vot")

tab

df = tab.to_pandas()

df.describe()

columns = ['ra', 'dec', 'parallax', 'pmra', 'pmdec']
g = sns.PairGrid(df, 
                 vars=columns, 
                )
g = g.map_diag(plt.hist, edgecolor="w")
g = g.map_upper(plt.scatter, alpha=0.5, edgecolor='w')
g = g.map_lower(sns.kdeplot, cmap="Blues_d")

# ## Now do a wider area

# This time it is 30 arcmin radius around 20:33:16 +41:18:41 

# Parallax error < 0.05 and proper motion limits: -3.5 to -1.5 in RA, -5.5 to -3.5 in Dec

ddf = Table.read("../data/Cygnus/Gaia/1553612017852O-result.vot").to_pandas()


m = (ddf.parallax > 0.3) & (ddf.parallax < 0.85)
ddf[m].describe()

sns.set_context("talk")
columns = ['ra', 'dec', 'parallax', 'pmra', 'pmdec']
g = sns.PairGrid(ddf[m], 
                 vars=columns, 
                )
g = g.map_diag(plt.hist, color="r", edgecolor="w")
g = g.map_upper(plt.scatter, s=15, color="r", alpha=0.4, edgecolor='none')
g = g.map_lower(sns.kdeplot, cmap="Reds_d")
g.fig.savefig("cygnus-x-gaia-pairplot.pdf")

ddf[m][columns].corr()

g.fig?


