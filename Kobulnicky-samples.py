# -*- coding: utf-8 -*-
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
sns.set_color_codes()

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

m = (tab02['RAh'] == 10)# & (tab01['RAm'] == 5)
tab02[m]

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
t['FIR_will'] = 1e-10*(2.4*np.abs(t['F8 or 12']) + 1.6*t['F24 or 22'] + 0.51*t['F70'])
t['ID', 'FIR', 'FIR_will']

import matplotlib.pyplot as plt
import seaborn as sns
# %matplotlib inline
sns.set_context("poster")

fig, ax = plt.subplots(figsize=(10, 8))
c = ax.scatter(t['FIR'], t['FIR_will'], 
               c=t['Dist'], cmap='magma_r', 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label('Distance, kpc')
for id_, x, y in zip(t['ID'], t['FIR'], t['FIR_will']):
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
    ylabel=r'Will: $F_\mathrm{IR}$',
)
ax.set_aspect('equal')
None

# So everything looks OK, except:
#
# * Source 67 is over 10 times too bright in the Kobulnicky table
#
# So, I will use my fluxes instead.  

# Now add in the table that I transcribed from the 2018 paper:

tab18 = Table.read('kob18.fits')

tt = join(t, tab18, keys='ID')
tt['LIR_will'] = tt['FIR_will']*(u.erg/u.cm**2/u.s)*4*np.pi*(tt['Dist']*1e3*u.parsec).to(u.cm)**2 / const.L_sun.to(u.erg/u.s)
tt['LIR/L* will'] = tt['LIR_will']/(1e4*tt['L4'])
tt

fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = 2.0/tt['L*/LIR_1'], 2*tt['LIR/L* will']
c = ax.scatter(xx, yy, 
               c=4.0 + np.log10(tt['L4']), cmap='viridis_r', 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label(r'$\log_{10}\ \left[L_* / L_\odot \right]$')
for id_, x, y in zip(tt['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
fmin, fmax = 2e-4, 5e-1
ax.plot([fmin, fmax], [fmin, fmax], ls='--')
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'Kobulnicky+ (2017): $\tau = 2 L_\mathrm{IR}/L_*$',
    ylabel=r'This paper: $\tau = 2 L_\mathrm{IR}/L_*$',
)
ax.set_aspect('equal')
None

# In this graph we compare the original luminosity ratio taken directly from the Kobulnicky (2017) table (x axis) with the ratio calculated using my new total IR fluxes, combined with the new luminosities in the Kobulnicky (2018) table (y axis).   Most of the points show reasonable agreement between the two methods, with a few exceptions:
#
# * 67: this had the $F_\text{IR}$ overestimated in K17.  Using a more reasonable value gives a lower $\tau$
# * 341: The spectral class has changed from B2 (K17) to O9 (K18), increasing the assumed $L_*$, which lowers $\tau$
# * 411: The luminosity class has changed from Ib (K17) to V (K18), so $L_*$ has been greatly reduced, which increases the estimated $\tau$
#
# And there doesn't seem to be any significant correlation with stellar luminosity.

# Next job is to estimate the shell pressure from the $\tau$:
#
# 1. Assume UV dust opacity per gas mass gives column density: $\Sigma = \tau/\kappa$.
# 2. Using measured thickness, find $\rho = \Sigma / H$
# 3. Assume sound speed, so $P_\mathrm{shell} = \rho a^2$
#
# Putting it all together gives $P_\mathrm{shell} = \tau a^2 / \kappa H$.

# Then we can compare that with the radiation pressure at the shell:
# $$
# P_\mathrm{rad} = \frac{L_*}{4\pi R^2 c}
# $$
# and define an observational momentum trapping efficiency: $\eta_\mathrm{obs} = P_\mathrm{shell} / P_\mathrm{rad}$, so that:
# $$
# \eta_\mathrm{obs} = \frac{\tau a^2}{\kappa H} \frac{4\pi R^2 c}{L_*}
# $$
# If we expand out the $\tau$, we se that $\eta_\mathrm{obs} \propto L_*^{-2}$, which is quite a steep sensitivity (especially since $L_*$ may be just a guess). 

# Make a new table that just has the columns that we want.  We take the $R_0$ from K18 Table 1.  The thickness $H$ could be taken from K17: "Height" in Tables 1 and 2. *But* I don't trust those values.  For instance, zeta Oph clearly has $H < 60''$ from its image, but the table gives $404''$, which is ridiculous given that $R_0 = 299''$. In fact, when I use these columns and calculate $H/R$, then I get a range from 0.5 to 3, which does not make any sense.

# It would be better to simply assume $H/R = 3 / (4 M^2)$, which is $\approx 0.1$ if $V = 30$ km/s, but more likely the velocities are lower.  Let's assume 0.25 for now.  Assume $\kappa = 600$ and $a = 11.4$ km/s

# +
colnames = tt.colnames[0:5] + ['L4', 'LIR_will', 'R0', 'D_kpc', 'U']
ttt = tt[colnames]
ttt['tau'] = np.round(2*tt['LIR/L* will'], decimals=5)
ttt['H/R'] =np.round(tt['Height'] / tt['R0_as'], decimals=2)

cs = (11.4*u.km/u.s).cgs
kappa = 600.0*u.cm**2/u.g
H_R = 0.25

ttt['P_k_shell'] = np.array(ttt['tau'])*(cs**2 / (const.k_B*kappa * H_R * ttt['R0']*u.parsec)).cgs
ttt['n_shell'] = np.round(ttt['P_k_shell']/u.Kelvin/(2 * 1.0e4), decimals=1)
ttt['P_k_rad'] = (1e4*const.L_sun*ttt['L4'] / (4*np.pi * (ttt['R0']*u.pc)**2 * const.c * const.k_B)).cgs
ttt['eta_obs'] = np.round(ttt['P_k_shell'].data/ttt['P_k_rad'].data, decimals=5)
ttt

# +
fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = ttt['tau'], ttt['eta_obs']
c = ax.scatter(xx, yy, 
               c=4.0 + np.log10(tt['L4']),
#               c=np.log10(tt['Teff']), 
               cmap='magma', vmin=4.0, vmax=6.0,
               edgecolors='k', alpha=1.0, s=300)
fig.colorbar(c, ax=ax).set_label(
    r'$\log_{10}\ \left[L_* / L_\odot \right]$'
#    r'$\log_{10}\ \left[T_\mathrm{eff} \right]$'
)
for id_, x, y in zip(tt['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize=10, color='k', 
        fontweight='black', fontstretch='condensed',
        xytext=(0,0), textcoords='offset points', ha='center', va='center',
    )
    ax.annotate(
        str(id_), (x, y), fontsize=10, color='w',
        xytext=(0,0), textcoords='offset points', ha='center', va='center',
    )



fmin, fmax = 8e-5, 8.0
ax.plot([fmin, fmax], [fmin, fmax], ls='--')
ax.fill_between([fmin, fmax], [25*fmin, 25*fmax], [1.5*fmin, 1.5*fmax], color='k', alpha=0.05)
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'UV optical depth of shell: $\tau$',
    ylabel=r'Shell momentum efficiency: $\eta_\mathrm{shell}$',
)
ax.set_aspect('equal')
fig.savefig('K18-eta-tau.pdf')
None
# -

# Now, before proceeding to the mass loss rates, let's do a factor plot of some selected parameters:

from scipy import stats

# +
df = ttt.to_pandas()
df.set_index(keys='ID', inplace=True)
columns = ['D_kpc', 'Teff', 'R0', 'L4', 'LIR_will', 'tau', 'n_shell', 'eta_obs']
pretty = [
    r'$\log_{10}\ D / \mathrm{kpc}$', 
    r'$\log_{10}\ T_\mathrm{eff} / \mathrm{K}$', 
    r'$\log_{10}\ R_0 / \mathrm{pc}$', 
    r'$\log_{10}\ L_* / L_\odot$', 
    r'$\log_{10}\ L_\mathrm{IR} / L_\odot$', 
#    r'$\log_{10}\ F_\mathrm{IR} / \mathrm{erg\ s^{-1}\ cm^{-2}}$', 
    r'$\log_{10}\ \tau$', 
    r'$\log_{10}\ n_\mathrm{shell} / \mathrm{cm^{-3}}$', 
    r'$\log_{10}\ \eta_\mathrm{obs}$',     
]
minmax = [ 
    [-1.3, 0.7], # D
    [4.1, 4.7], # Teff
    [-2.1, 0.4], # R0
    [3.7, 6.3], # was L4 (but now just L/Lsun)
    [0.3, 3.6], # LIR
#    [-10.3, -6.7], # FIR
    [-3.8, -1.2], # tau
    [-0.5, 3.1], # n_shell
    [-3.5, -0.4], # eta_obs
]
df = df.assign(**{col: np.log10(df[col]) for col in columns})
df = df.assign(
    L4=4.0 + df.L4,
    close=pd.Categorical((df['D_kpc'] < 1.5).astype('S5')),
)

def corrfunc(x, y, **kws):
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(x, y)
#    r, p = stats.pearsonr(x, y)
    ax = plt.gca()
    if pvalue < 0.05:
        fontcolor = 'k' if pvalue < 0.003 else 'k'
        fontalpha = 1.0 if pvalue < 0.003 else 0.5
        ax.annotate(f"$r = {rvalue:.2f}$",
                    xy=(.05, .9), xycoords=ax.transAxes, 
                    fontsize='small', color=fontcolor, alpha=fontalpha)
        ax.annotate(rf"$m = {slope:.2f} \pm {stderr:.2f}$",
                    xy=(.95, .05), ha='right', xycoords=ax.transAxes, 
                    fontsize='small', color=fontcolor, alpha=fontalpha)

    
with sns.plotting_context('talk', font_scale=1.0):
    g = sns.pairplot(df, vars=columns, 
                     diag_kind='hist', 
#             diag_kind='kde',
#             hue='close', 
                     kind='reg',
                     plot_kws=dict(scatter_kws=dict(edgecolor='w')),
                     diag_kws=dict(bins=6),
                    )
    g.map_offdiag(corrfunc)
    for j, [[v1, v2], label] in enumerate(zip(minmax, pretty)):
        g.axes[j, j].set(xlim=[v1, v2], ylim=[v1, v2])
        g.axes[-1, j].set(xlabel=label)
        g.axes[j, 0].set(ylabel=label)
    plt.gcf().savefig('K18-pairplot.pdf')
    
#df
# -

# Some observations about these correlations:
#
# 1. The basic observational parameters are $D$, $T_\mathrm{eff}$, $R_0$, $L_*$, $L_\mathrm{IR}$.  In principle, observational errors in these should all be independent, so any correlations between them are real. 
#     1. The first three are the best determined, and there are no significant corrlations between them.
#     2. The luminosities might potentially suffer from Malmquist bias, but there is no evidence for this in the case of $L_*$ since low luminosities are seen at the largest distance.  In the case of $L_\mathrm{IR}$, there is a slight increasing trend with distance, which might indicate bias, but it is only marginally significant.
#     3. $R_0$ could also be Malmquist biased due to angular resolution, but there is no trend in $R_0$–$D$, so it seems not.  This is probably because the smallest bows in this subsample have angular sizes $> 4''$ and are still just resolved at Carina at 8 micron (Spitzer has $2''$ resolution at 8 micron).
#
# 2. $\tau$ is found from $L_\mathrm{IR}/L_*$, meaning that errors in $L_*$ would give negative $\tau$–$L_*$ correlation, but they seem to be uncorrelated so we are probably OK.  
#
# 3. If $L_*$ were determined from photometry, then any distance errors would give a spurious $R_0 \propto L_*^{1/2}$ correlation.  **So it is a good job that it is found from spectral classification instead**.  Also, the distances should be reliable.
#
# 4. What we actually see is a linear $R_0 \propto L_*$ correlation.  This is the most significant correlation of all ($r = 0.79$), implying that >60% of the variation in $R_0$ is driven by the stellar luminosity, leaving only 40% left for the environment (and other stellar properties). 
#
# 5. We also have a moderate correlation between $L_*$ and $T_\mathrm{eff}$, which is clearly because half the stars are MS, with a tight correlation, and the other half are evolved
#
# 6. Other strong correlations:
#     1. $L_\mathrm{IR}$–$L_*$: approximately linear, explaining 45% in variance of $L_\mathrm{IR}$.  This means that there is no trend in $\tau$ with $L_*$
#     2. $\tau$–$L_\mathrm{IR}$ is very well correlated with $m \approx -1$.  **However**, this is exactly what you would get if the dispersion about the linear $L_\mathrm{IR}$–$L_*$ trend were entirely due to measurement errors in $L_\mathrm{IR}$.  Hopefully, that is not the case, but it needs checking.  The std dev of log $\tau$ is about 0.5 dex (the same is true of nearly all the parameters, except for $T_\mathrm{eff}$).
#

# +
df.describe()

# As an aside, we will compare my $G$ with their $U$
ttt['G'] = 0.074*200*ttt['L4']/ttt['R0']**2
ddf = ttt['ID', 'U', 'G'].to_pandas()
fig, ax = plt.subplots(figsize=(10, 10))
vmin, vmax = 150, 4e5
ax.scatter(x='U', y='G', data=ddf)
ax.plot([vmin, vmax], [vmin, vmax])
for id_, x, y in zip(ttt['ID'], ttt['U'], ttt['G']):
    ax.annotate(
        str(id_), (x, y),
        xytext=(4,4), textcoords='offset points',
               )
ax.set(xscale='log', yscale='log', 
       xlim=[vmin, vmax], ylim=[vmin, vmax],
       xlabel='U 2017', ylabel='U 2018',
      )
ax.set_aspect('equal')
# -

ddf['log U/U'] = np.log10(ddf.U/ddf.G)
ddf.describe()

# So this seems to be a discrepancy between their U column and their new luminosities. So maybe I shouldn't worry about it.

from astropy.table import QTable
ltab = QTable(tt['ID', 'Teff', 'R*', 'L4'], masked=False)
ltab['L2017'] = (4*np.pi*ltab['R*']**2 * const.sigma_sb*ltab['Teff']**4).to(u.Lsun)
ltab['L4'] = ltab['L4']*1e4*u.Lsun
ltab

fig, ax = plt.subplots(figsize=(10, 10))
vmin, vmax = 3e3, 2e6
ax.scatter(ltab['L2017'], ltab['L4'])
for id_, x, y in zip(ltab['ID'], ltab['L2017'].data, ltab['L4'].data):
    ax.annotate(
        str(id_), (x, y),
        xytext=(4,4), textcoords='offset points',
               )
ax.plot([vmin, vmax], [vmin, vmax])
ax.set(xscale='log', yscale='log', 
       xlim=[vmin, vmax], ylim=[vmin, vmax],
       xlabel='L2017', ylabel='L2018',
      )
ax.set_aspect('equal')

# Not the same!  But doesn't fully explain the $U$ scatter.  Must be different radii too


