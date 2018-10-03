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

import numpy as np
from astropy.table import Table
import pandas as pd
from matplotlib import pyplot as plt
# %matplotlib inline
import seaborn as sns

# ### Read in bow shock data from Kobulnicky:2017a

tab01 = Table.read("data/Kobulnicky2017/J_AJ_154_201_table1.dat.fits")

# Remove coordinate columns that we do not want. Also the color temperatures. *And we have to remove the Name field*. This is because of inconsitencies between the names in the different tables, which cause problems when we merge.

tab01.remove_columns(
   [
       'RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs',
       'T70/160', 'Name',
   ]
)

tab01[300:310]

# Indicate lower limit fluxes by negative values, as in the published table: 

for band in 'F3.6', 'F4.5', 'F5.8', 'F8.0', 'F70', 'F160':
    lband = 'l_' + band
    m = tab01[lband] == '<'
    tab01[band][m] = -tab01[band][m]
    tab01.remove_column(lband)

tab01

tab02 = Table.read("data/Kobulnicky2017/J_AJ_154_201_table2.dat.fits")

tab02.remove_columns(
   [
       'RAh', 'RAm', 'RAs', 'DE-', 'DEd', 'DEm', 'DEs',
       'T70/160', 'Name',
   ]
)
for band in 'F3.4', 'F4.6', 'F12', 'F70', 'F160':
    lband = 'l_' + band
    m = tab02[lband] == '<'
    tab02[band][m] = -tab02[band][m]
    tab02.remove_column(lband)

tab02[20:30]

# There is something wrong with Table 5.  There is an integer column that contains `'---'` strings, which need dealing with or the reader will crash. This doesn't seem to be possible with the FITS table reader, so we resort to ascii, where we can fix it with the `fill_values` parameter.

tab05 = Table.read("data/Kobulnicky2017/table5.dat", 
                   fill_values=[('---', -1)], 
                   format='ascii.cds', 
                   readme="data/Kobulnicky2017/ReadMe")

tab05.remove_columns(
    [
        'TSS', 'T22/T70', 'T70/160', 
    ]
)

import astropy.units as u

import astropy.constants as const

tab05['LIR'] = tab05['FIR'].to(u.erg/u.cm**2/u.s)*4*np.pi*(tab05['Dist']*1e3*u.parsec).to(u.cm)**2 / const.L_sun.to(u.erg/u.s)

tab05

from astropy.table import join

tab05_01 = join(tab05, tab01, keys=('ID'), join_type='left')
tab05_01.remove_columns(['F3.6', 'F4.5', 'F5.8',])
tab05_01

tab05_01_02 = join(tab05_01, tab02, keys=('ID'), join_type='left')
tab05_01_02.remove_columns(['F3.4', 'F4.6',])
tab05_01_02

# Now merge the WISE and Spitzer photometry, taking (8, 12) and (22, 24) as equivalent.

# Make a mask that is true for rows with Spitzer photometry
m_sst = ~tab05_01_02['Rad_1'].mask
m_wise = ~tab05_01_02['Rad_2'].mask
m_sst, m_wise

# +
groups = [
    ['Rad_1', 'Rad_2', 'Rad'],
    ['Height_1', 'Height_2', 'Height'],
    ['F8.0', 'F12', 'F8 or 12'],
    ['F24', 'F22', 'F24 or 22'],
    ['F70_1', 'F70_2', 'F70'],
    ['F160_1', 'F160_2', 'F160'],
    ['I70_1', 'I70_2', 'I70'],
    ['T24/70', 'T22/70', 'T2x/70'],
]

for sst, wise, merge in groups:
    tab05_01_02[merge] = np.where(m_sst, tab05_01_02[sst], np.where(m_wise, tab05_01_02[wise], np.nan))
    tab05_01_02[merge].mask = ~(m_sst | m_wise)
    tab05_01_02.remove_columns([sst, wise])
tab05_01_02['Observatory'] = np.where(m_sst, 'SST', np.where(m_wise, 'WISE', None))
tab05_01_02
# -

# Now work out my own IR flux by weighted sum of the 8 to 160 bands

t = tab05_01_02
t.remove_row(2)
t['FIR_will'] = 1e-10*(2.55*np.abs(t['F8 or 12']) + 0.77*t['F24 or 22'] + 0.24*t['F70'])
t['ID', 'FIR', 'FIR_will']

import matplotlib.pyplot as plt
import seaborn as sns
# %matplotlib inline
sns.set_context("poster")

fig, ax = plt.subplots(figsize=(10, 8))
c = ax.scatter(t['FIR'], 2*t['FIR_will'], 
               c=t['Dist'], cmap='magma_r', 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label('Distance, kpc')
for id_, x, y in zip(t['ID'], t['FIR'], 2*t['FIR_will']):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
fmin, fmax = 2e-10, 2e-7
ax.plot([fmin, fmax], [fmin, fmax], ls='--')
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'Kobulnicky: $F_\mathrm{IR}$',
    ylabel=r'Will: $2 \times F_\mathrm{IR}$',
)
ax.set_aspect('equal')
None

# So everything looks OK, except:
#
# 1. I had to multiply my fluxes by 2
# 2. Source 67 is over 10 times too bright in the Kobulnicky table

# Now add in the table that I transcribed from the 2018 paper:

tab18 = Table.read('kob18.fits')

join(t, tab18, keys='ID')


