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

tab01[20:30]





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

# Now work out my own IR flux by weighted sum of the 8 to 160 bands.  Originally, here we removed 329 because it lacks the requisite data, being absent from Tables 1 and 2.  However, it is in Table 5 and also in K18, so we have a flux. 

t = tab05_01_02
#t.remove_row(2)
t['FIR_will'] = 1e-10*(2.4*np.abs(t['F8 or 12']) + 1.6*t['F24 or 22'] + 0.51*t['F70'])
t[2]['FIR_will'] = t[2]['FIR']
t['ID', 'FIR', 'FIR_will']

import matplotlib.pyplot as plt
import seaborn as sns
# %matplotlib inline
sns.set_context("poster")

fig, ax = plt.subplots(figsize=(10, 8))
c = ax.scatter(t['FIR'], t['FIR_will'], 
               c=t['Dist'], cmap='viridis_r', vmin=0.0, vmax=2.5,
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
    xlabel=r'Kobulnicky+ (2018): $F_\mathrm{IR}$, mW/m$^2$',
    ylabel=r'This paper: $F_\mathrm{IR}$, mW/m$^2$',
)
ax.set_aspect('equal')
fig.savefig('K18-flux-comparison.pdf')
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
               c=4.0 + np.log10(tt['L4']), cmap='magma', vmin=3.5, vmax=6.0, 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label(r'$\log_{10}\ \left[L_* / L_\odot \right]$')
for id_, x, y in zip(tt['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
fmin, fmax = 2e-4, 5e-1
ax.fill_between([fmin, fmax], [fmin/2, fmax/2], [fmin*2, fmax*2], color='k', alpha=0.2, zorder=-1)
ax.plot([fmin, fmax], [fmin, fmax], ls='--', zorder=-1)
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'Kobulnicky+ (2017): $\tau = 2 L_\mathrm{IR}/L_*$',
    ylabel=r'This paper: $\tau = 2 L_\mathrm{IR}/L_*$',
)
ax.set_aspect('equal')
fig.tight_layout()
fig.savefig('K17-tau-comparison.pdf')
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
#
# _But maybe "height" is not what I think it is._ What is really wanted is the FWHM of the brightness profile, but H is measured at a low contour (small fraction of the peak brightness), so it is too large. 

# It would be better to simply assume $H/R = 3 / (4 M^2)$, which is $\approx 0.1$ if $V = 30$ km/s, but more likely the velocities are lower.  Let's assume 0.25 for now.  Assume $\kappa = 600$ and $a = 11.4$ km/s

# +
colnames = tt.colnames[0:5] + ['L4', 'LIR_will', 'R0', 'D_kpc', 'U', 'Md_-8', 'V3']
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

# Now do the same, but plot $\eta/\tau$ against $\tau$ to remove the linear trend. 

# +
fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = ttt['tau'], ttt['eta_obs']/ttt['tau']
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

ax.axhline(1.5, ls='--')
#ax.plot([fmin, fmax], [fmin, fmax], ls='--')
#ax.fill_between([fmin, fmax], [25*fmin, 25*fmax], [1.5*fmin, 1.5*fmax], color='k', alpha=0.05)
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[0.3, 50.0],
    xlabel=r'UV optical depth of shell: $\tau$',
    ylabel=r'Shell pressure boost: $P_\mathrm{shell} / \tau P_\mathrm{rad}$',
)
#ax.set_aspect('equal')
fig.savefig('K18-eta-tau-compensated.pdf')
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
        fontcolor = 'purple' if pvalue < 0.003 else 'k'
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

df.describe()

# ## Finally the mass loss rates

# In my approach, $\dot M$ comes from balancing the shell pressure with the ram pressure: $\eta_{\mathrm{sh}} = \eta_{\mathrm{w}}$.  So 
# $$
# \dot M_{-7} = 2.02 L_4 \eta / V_3
# $$

ttt['Md'] = 2.02e-7*ttt['L4']*ttt['eta_obs'] / ttt['V3']
ttt['Md_K18'] = 1e-8*ttt['Md_-8']

fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = ttt['Md_K18'], ttt['Md']
c = ax.scatter(xx, yy, 
               c=4.0 + np.log10(tt['L4']), cmap='magma', vmin=4.0, vmax=6.0, 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label(r'$\log_{10}\ \left[L_* / L_\odot \right]$')
for id_, x, y in zip(tt['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
fmin, fmax = 3e-10, 3e-6
ax.plot([fmin, fmax], [fmin, fmax], ls='--')
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'Kobulnicky+ (2018): $\dot M$, M$_\odot$/yr',
    ylabel=r'This paper: $\dot M$, M$_\odot$/yr',
)
ax.set_aspect('equal')
fig.savefig('K18-mdot-comparison.pdf')
None

# Correct for radiation pressure on shell. And calculate errors, assuming a factor 3 uncertainty either way in $\eta$

ttt['Md_t'] = 2.02e-7*ttt['L4']*(ttt['eta_obs'] - 1.25*ttt['tau'])/ ttt['V3']
efactor = 3.0
ttt['err Md+'] = ttt['Md']*(efactor - 1.0)
ttt['err Md-'] = ttt['Md']*(1 - 1/efactor)

fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = ttt['Md_K18'], ttt['Md_t']
c = ax.scatter(xx, yy, 
               c=4.0 + np.log10(tt['L4']), cmap='magma', vmin=4.0, vmax=6.0, 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label(r'$\log_{10}\ \left[L_* / L_\odot \right]$')
for id_, x, y in zip(tt['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
fmin, fmax = 1e-10, 3e-6
ax.plot([fmin, fmax], [fmin, fmax], ls='--')
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'Kobulnicky+ (2018): $\dot M$, M$_\odot$/yr',
    ylabel=r'This paper: $\dot M$, M$_\odot$/yr',
)
ax.set_aspect('equal')
fig.savefig('K18-mdot-tau-corrected-comparison.pdf')
None

fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = ttt['Md_t'], ttt['Md_t'] - ttt['err Md-']
c = ax.scatter(xx, yy, 
               c=4.0 + np.log10(tt['L4']), cmap='magma', vmin=4.0, vmax=6.0, 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label(r'$\log_{10}\ \left[L_* / L_\odot \right]$')
for id_, x, y in zip(tt['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
fmin, fmax = 1e-10, 3e-6
ax.plot([fmin, fmax], [fmin, fmax], ls='--')
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    ylabel=r'$\dot M - \epsilon$',
    xlabel=r'This paper: $\dot M$, M$_\odot$/yr',
)
ax.set_aspect('equal')
fig.savefig('mdot-vs-negative-error.pdf')
None

mask_upper_limits = ttt['Md_t'] - ttt['err Md-'] < 0.0
ttt['rel err'] = ttt['err Md-']/ttt['Md_t']
ttt[mask_upper_limits]['ID', 'Md_t', 'err Md-', 'rel err']

# Those are the sources where Mdot becomes an upper limit if we use a factor 3 uncertainty in $\eta$.   So, they don't appear on the above graph since Mdot - epsilon is negative. 

ttt['ID', 'Md_t', 'Md_K18']

# Why do we have such a lack of correlation?  Is it my fault or theirs?  Theirs, I hope.

ttt['LIR_K18'] = tt['L4']/tt['L*/LIR_2']
ttt['LIR_K17'] = tt['L4']/tt['L*/LIR_1']
ttt['R/L'] = tt['R0']/tt['L4']
ttt['LIR/V'] = tt['LIR_will']/tt['V3']
ttt['RL/LV'] = tt['R0']*tt['LIR_will']/(tt['L4']*tt['V3'])
ttt['MdMd'] = ttt['Md_K18']/ttt['Md']
ttt['theta'] = ttt['R0']/ttt['D_kpc']/206.264
ttt['FIR'] = tt['FIR_will']
columns = ['D_kpc', 'Teff', 'theta', 'R0', 'L4', 'R/L', 'LIR/V', 'FIR', 'LIR_will', 'LIR_K18', 'LIR_K17', 
           'U', 'tau', 'n_shell', 'eta_obs', 'V3', 'RL/LV', 'Md', 'Md_K18', 'MdMd']
mdf = ttt[columns].to_pandas()
mdf = mdf.assign(**{col: np.log10(mdf[col]) for col in mdf.columns})
mdf.corr()

# So, $r = 0.7 \therefore r^2 = 0.477$ between the two $\dot M$ estimates, meaning only half the variance in one is explained by the other. 
#

# My $\dot M$ comes explicitly from $L_* \eta / V$.  The respective $r^2$ are 0.23, 0.74, 0.06, which sum to 1.03, not 1, due to mutual covariance between $L_*$ and $V$ ($\eta$ is totally uncorrelated with either, interestingly).  So my $\dot M$ is _mostly_ determined by $\eta$, but $\eta$ itself comes from $R \tau /L_*$ and $\tau$ comes from $L_\mathrm{IR} / L_*$, so expanding everything: $\eta$ comes from  $R L_\mathrm{IR}/ L_*^2$ so 
# $$
# \dot M \propto \frac{R L_\mathrm{IR}}{L_* V} \ .
# $$

# *Note* we made a column `RL/LV` to check that it is the same as $\dot M$ -- thankfully _it is!_  

# Rather than the `.corr()`, we will use `.cov()` to get the raw covariance matrix.  This is more useful in some respects, since it more clearly shows how much variance there is in the individual parameters.  For instance, stellar wind velocity has very low dispersion.  Since we have taken log of all quantities, there is no problem with disparate numerical scales. 

mdf[['LIR_will', 'L4', 'R0', 'R/L', 'V3', 'Md']].cov()

# So, from that it is clear that _operationally_ the principal determinant of my $\dot M$ is the infrared luminosity.  Also look at correlation matrix for good measure.

mdf[['LIR_will', 'L4', 'R0', 'R/L', 'V3', 'Md']].corr()

# So, infrared luminosity determines 75% of $\dot M$ variance, while $V$ clearly has no influence on anything (< 6% level). Interestingly neither does $R/L_*$ even though that is part of the operational determinant of $\dot M$ and $R$ and $L$ do individually contribute to $\dot M$ variance at 25% level.  Part of the reason must be simply that $R/L_*$ does not vary much, having a rms dispersion of 0.3 dex (about factor 2). 

# So, where does the extra 25% variance in $\dot M$ come from?  There is a weak negative correlation of $L_{\mathrm{IR}}$ with $R/L_*$ that might have something to do with it.  This means that mass loss rate shows less dispersion (0.55 dex) than the IR luminosity (0.62 dex).  Whatever, it doesn't really matter.
#
# Note that the diagonal elements of the covariance array are the variances, so std is sqrt of that.

# _Alternatively ..._ we can put it in exclusively empirical terms: $F_{\mathrm{IR}}$, $D$, $\theta$.
# $$
# \dot M \propto \frac{\theta D^3 F_{\mathrm{IR}}}{L_* V}  
# $$

mdf[['D_kpc', 'theta', 'FIR', 'L4', 'V3', 'Md']].corr()

# But, as can be seen, $\dot M$ is correlated with none of these very well. So, best stick with $L$ and $R$. 

# ### Now, look at their Mdots

# What do they really depend on?  
#
# 1. They use the LOS column density, instead of the radial column density, which they estimate from the 70 micron surface brightness, together with an estimate of the 70 micron emissivity, $j_\nu$, which depends on the radiation field $U$.  They claim that $j_\nu \propto U^{1/2}$ ($j_\nu$ is emissivity 
#
# 2. They skip out the middle man of the shell and directly balance the ram pressure of the wind with the ram pressure of the ambient stream. _Except that not really, since they use the density of the shell and then say that it is factor of 4 times the ambient density._
#
# 3. They just assume a stream velocity of 30 km/s, so they don't depend on the gas temperature in the shell (although, in reality it should effect the compression ratio). 
#
# $$
# \frac{\dot M V}{4 \pi R^2} = (\rho_s / \delta) v_a^2  
# $$
# Also, apparently 
# $$
# \rho_s = m I_\nu / \ell j_\nu
# $$
# From their equation (8) this gives the following dependency for $\dot M$:
# $$
# \dot M = 4 \pi \left[ \frac{R^2 I_\nu}{\ell V_w}  \right] \left[ \frac{V_a^2 m }{j_\nu(U) \delta } \right]  
# $$
# where first term in square brackets is measured (or at least estimated per star) quantities, while second term is in square brackets is assumed parameters, except that $j_\nu (U)$ is the dust emissivity (Jy.cm2/sr/nucleon = 1e-23 erg/s/sr/nucleon) as a function of radiation field where $U F_0 = L / 4 \pi R^2$ with $F_0 = 0.0217$ erg/cm2/s being the ISRF flux, but integrated across the whole SED. 

# We show below in the next section that (1) they are not using the $j_\nu(U)$ that they claim to be using, and (2) they should be using a somewhat different one anyway (because of SED shape).  However, we have corrected that, so we can proceed with the analysis.  As in K18, approximate the $U$ dependence as $U^{1/2} = L_*^{1/2} R^{-1}$. So, the $\dot M$ breakdown becomes:
# $$
# \dot M \propto \left[ \frac{R^3 I_\nu}{\ell L_*^{1/2} V_w}  \right] \left[ \frac{V_a^2 m }{j_0 \delta } \right]
# $$

# Or
# $$
# \dot M \propto \left[ \frac{R I_\nu}{(\ell/R) U^{1/2} V_w}  \right] \left[ \frac{V_a^2 m }{j_0 \delta } \right]
# $$

# ##  An investigation of radiation field, $U$, dust emissivity, $j_\nu$, and mass loss rate for K17 and K18

# _Revisiting this, now that I have established that the DL07 emissivities are too low for OB stars if one uses bolometric flux for the U (because the ISRF has only a small fraction at UV wavelengths)._


# What is the `U` column in the table?  It comes from Table 5 of K17 and supposedly comes from the $T_{\mathrm{eff}}$, $R_*$ and $R_0$ from the same table.   We have this table, so we can check this.

tab05['UU'] = 0.01329*(tab05['Teff']/1e4)**4 * tab05['R*']**2 / tab05['Dist2']**2
tab05

ddf = tab05['ID', 'U', 'UU'].to_pandas()
fig, ax = plt.subplots(figsize=(10, 10))
vmin, vmax = 150, 4e5
ax.scatter(x='U', y='UU', data=ddf)
ax.plot([vmin, vmax], [vmin, vmax])
for id_, x, y in zip(tab05['ID'], tab05['U'], tab05['UU']):
    ax.annotate(
        str(id_), (x, y),
        xytext=(4,4), textcoords='offset points',
               )
ax.set(xscale='log', yscale='log', 
       xlim=[vmin, vmax], ylim=[vmin, vmax],
       xlabel='U 2017', ylabel='UU 2017',
      )
ax.set_aspect('equal')

# So the previous plot shows the `U` column from K17 Tab 5 on the x axis, and the $U$ that I calculate from the `Teff`, `R*`, and `Dist2` (equals radius) columns of the same table on the y axis. 
#
# It is strange that there is any scatter at all, but this shows that there isn't a serious problem with the way that K17 calculated their $U$. 

# Next we compare with the K18 values, to work out where they got their $U$ from.  We have got the original tables from the ApJ website and exported them to FITS (see `dust-wave-case-studies.org`).
#
# _Note that we have to explicitly change some columns from string to float._

k18tab1 = Table.read('data/Kobulnicky2018/k18tab1.fits')
k18tab2 = Table.read('data/Kobulnicky2018/k18tab2.fits')
k18tab = join(k18tab1, k18tab2, keys=('ID', 'Name', 'Alt. name'), join_type='left')
for col in 'U', 'j_nu', 'Mdot':
    k18tab[col] = k18tab[col].astype('float')
k18tab

# First, check the U values against what they should be: $U = 14.7 L_4 R_{\mathrm{pc}}^{-2}$. 

k18tab['UU'] = 14.7*k18tab['Lum.']/k18tab['R_0']**2
ddf = k18tab['ID', 'U', 'UU'].to_pandas()
fig, ax = plt.subplots(figsize=(10, 10))
vmin, vmax = 150, 4e5
ax.scatter(x='U', y='UU', data=ddf)
ax.plot([vmin, vmax], [vmin, vmax])
for id_, x, y in zip(k18tab['ID'], k18tab['U'], k18tab['UU']):
    ax.annotate(
        str(id_), (x, y),
        xytext=(4,4), textcoords='offset points',
               )
ax.set(xscale='log', yscale='log', 
       xlim=[vmin, vmax], ylim=[vmin, vmax],
       xlabel='U 2018', ylabel='UU 2018',
      )
ax.set_aspect('equal')
(ddf['U']/ddf['UU']).describe()

# So, it looks like they are using a factor of about 16.3 instead of 14.7, which is odd.  This makes their values about 1.1 times higher than mine.  But apart from that, the U values look fine.

# Now, look at the emissivity as a function of $U$, and compare it with the DL07 values.

DL07tab = Table.read('../cloudy-dust-charging/DL07-data/emissivities.fits')

ddf = k18tab['ID', 'U', 'j_nu'].to_pandas()
fig, ax = plt.subplots(figsize=(8, 7))
umin, umax = 110, 5e5
jmin, jmax = 2e-13, 5e-11
ax.plot(k18tab['U'], k18tab['j_nu'], 'x', alpha=1.0, markeredgewidth=2, color='r', label='K18 sources')
ax.plot(DL07tab['U'], DL07tab['70'], '-', color='k', alpha=0.3, lw=10, label='DL07 models')
ax.plot(DL07tab['U']/8.0, DL07tab['70'], ':', color='k', alpha=0.3, lw=10, label=r'DL07($U \times 8$)')
ax.legend(loc='lower right')
ax.set(xscale='log', yscale='log', 
       xlim=[umin, umax], 
       ylim=[jmin, jmax],
       xlabel=r'Radiation field: $U = F / F_{\mathrm{MMP83}}$', 
       ylabel=r'70 $\mu$m emissivity: $j_\nu$, Jy cm$^{2}$ sr$^{-1}$ H$^{-1}$',
      )
sns.despine()
fig.tight_layout()
fig.savefig('K18-emissivity-vs-U.pdf')
None

# So, they are not even using the emissivities that they say they are using.  But we need to check if it is just a problem with the table, or if they are actuslly using these values to calculate the $\dot M$.  _Yes, they are – see below._

# So, we will calculate the $\dot M$ from their table quantities, using their equation (8). 

k18tab['Va'] = 30.0
k18tab['Va'][0] = 26.5
k18tab['Mdot2'] = 1.67e-28*k18tab['R_0,as']**2 * k18tab['D'] * k18tab['Va']**2 * 1e7*k18tab['Peak_70'] / (k18tab['V_inf_{}'] * k18tab['ell'] * k18tab['j_nu'])

fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = k18tab['Mdot'], k18tab['Mdot2']
c = ax.scatter(xx, yy, 
               c=4.0 + np.log10(k18tab['Lum.']), cmap='magma', vmin=4.0, vmax=6.0, 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label(r'$\log_{10}\ \left[L_* / L_\odot \right]$')
for id_, x, y in zip(k18tab['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
fmin, fmax = 1e-9, 3e-6
ax.plot([fmin, fmax], [fmin, fmax], ls='--')
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'From K18 Table 2, $\dot M$',
    ylabel=r'From K18 eq. (8), $\dot M$',
)
ax.set_aspect('equal')
fig.savefig('K18-mdot-internal-comparison.pdf')
None

k18tab['Md/Md'] = k18tab['Mdot2']/k18tab['Mdot']
print(f"Ratio of Mdot = {k18tab['Md/Md'].mean():.4f} +/- {k18tab['Md/Md'].std():.4f}")


# So, they do seem to be using the emissivity values from their Table.  Although there is a strange factor of 0.83 appearing from somwhere.

# ### What the mass loss rates _should_ have been in K18

# First, we will use the $j_\nu$ values directly from DL07, which is what they claim they were doing (although they were not, see above).  

k18tab['jnu_dl07'] = np.interp(k18tab['U'], DL07tab['U'], DL07tab['70'])

# Check that that worked:

fig, ax = plt.subplots(figsize=(10, 10))
umin, umax = 150, 5e5
jmin, jmax = 5e-13, 5e-11
ax.plot(k18tab['U'], k18tab['jnu_dl07'], 'o', alpha=0.5, label='K18 sources')
ax.plot(DL07tab['U'], DL07tab['70'], '-', color='k', alpha=0.3, lw=10, label='DL07 models')
ax.plot(DL07tab['U']/8.0, DL07tab['70'], ':', color='k', alpha=0.3, lw=10, label=r'DL07($U \times 8$)')
for i, [id_, x, y] in enumerate(zip(k18tab['ID'], k18tab['U'], k18tab['jnu_dl07'])):
    ax.annotate(
        str(id_), (x, y),
        xytext=(np.random.randint(-10, 10), 35 - 10*(5*i % 7) + np.random.randint(0, 15)), 
        textcoords='offset points', 
        fontsize=10, ha='center', alpha=0.8,
               )
ax.legend()
ax.set(xscale='log', yscale='log', 
       xlim=[umin, umax], 
       ylim=[jmin, jmax],
       xlabel=r'K18 radiation field: $U$', 
       ylabel=r'70 micron emissivity: $j_\nu$',
      )
None

k18tab['Mdot3'] = 1.67e-28*k18tab['R_0,as']**2 * k18tab['D'] * k18tab['Va']**2 * 1e7*k18tab['Peak_70'] / (k18tab['V_inf_{}'] * k18tab['ell'] * k18tab['jnu_dl07'])

fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = k18tab['Mdot'], k18tab['Mdot3']
c = ax.scatter(xx, yy, 
               c=4.0 + np.log10(k18tab['Lum.']), cmap='magma', vmin=4.0, vmax=6.0, 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label(r'$\log_{10}\ \left[L_* / L_\odot \right]$')
for id_, x, y in zip(k18tab['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
fmin, fmax = 1e-10, 3e-6
ax.plot([fmin, fmax], [fmin, fmax], ls='--')
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'From K18 table, $\dot M$',
    ylabel=r'Corrected to DL07, $\dot M$',
)
ax.set_aspect('equal')
fig.savefig('K18-mdot-DL07-comparison.pdf')
None

# Next, we will correct the emissivities to account for the different SED between O stars and ISRF.  We do this by multiplying $U$ by eight.

k18tab['jnu_Ux8'] = np.interp(8*k18tab['U'], DL07tab['U'], DL07tab['70'])
k18tab['Mdot4'] = 1.67e-28*k18tab['R_0,as']**2 * k18tab['D'] * k18tab['Va']**2 * 1e7*k18tab['Peak_70'] / (k18tab['V_inf_{}'] * k18tab['ell'] * k18tab['jnu_Ux8'])

fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = k18tab['Mdot'], k18tab['Mdot4']
c = ax.scatter(xx, yy, 
               c=np.log10(k18tab['U']), cmap='Blues_r', vmin=2.3, vmax=5.3, 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label(r'$\log_{10}\ U$')
for id_, x, y in zip(k18tab['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
fmin, fmax = 1e-10, 3e-6
ax.plot([fmin, fmax], [fmin, fmax], ls='--', zorder=-1)
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'K18 published $\dot M$',
    ylabel=r'K18 corrected to DL07 with $U \times 8$, $\dot M$',
)
ax.set_aspect('equal')
fig.tight_layout()
fig.savefig('K18-mdot-Ux8-comparison.pdf')
None

# ### Comparison between my mass loss and the K18 corrected values

fig, ax = plt.subplots(figsize=(10, 8))
xx, yy, ye = k18tab['Mdot4'], ttt['Md_t'], ttt['err Md']
ax.errorbar(xx, yy, yerr=ye, fmt='none')
c = ax.scatter(xx, yy, 
               c=4.0 + np.log10(tt['L4']), cmap='magma', vmin=4.0, vmax=6.0, 
               edgecolors='k', alpha=1.0)
fig.colorbar(c, ax=ax).set_label(r'$\log_{10}\ \left[L_* / L_\odot \right]$')
for id_, x, y in zip(tt['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
fmin, fmax = 1e-10, 3e-6
ax.plot([fmin, fmax], [fmin, fmax], ls='--')
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'Kobulnicky+ (2018) corrected: $\dot M$, M$_\odot$/yr',
    ylabel=r'This paper: $\dot M$, M$_\odot$/yr',
)
ax.set_aspect('equal')
fig.tight_layout()
fig.savefig('K18-mdot-corrected-comparison.pdf')
None

# +
fig, ax = plt.subplots(figsize=(10, 8))


for m, alpha in [mask_upper_limits, 0.3], [~mask_upper_limits, 1.0]:
    xx, yy, ye = (k18tab[m]['Mdot4'], ttt[m]['Md_t'], 
                  np.stack((ttt[m]['err Md-'], ttt[m]['err Md+'])))
    ax.errorbar(xx, yy, yerr=1.0*ye, alpha=alpha,
                fmt='none', zorder=-1, color='b')
    c = ax.scatter(xx, yy, 
                   c=np.log10(tt[m]['R0']), cmap='YlOrRd_r', 
                   vmin=-1.5, vmax=0.0, 
                   edgecolors='k', alpha=alpha)
    for id_, x, y in zip(tt[m]['ID'], xx, yy):
        ax.annotate(
            str(id_), (x, y), fontsize='xx-small', 
            xytext=(4,4), textcoords='offset points', alpha=alpha)
        
fig.colorbar(c, ax=ax).set_label(r'$\log_{10}\ \left[R_0 / \mathrm{pc} \right]$')    
    
fmin, fmax = 1e-10, 3e-6
ax.plot([fmin, fmax], [fmin, fmax], 
        ls='--', zorder=-1, color=(0.9, 0.6, 0.1))
ax.set(
    xscale='log', yscale='log', 
    xlim=[fmin, fmax], ylim=[fmin, fmax],
    xlabel=r'Kobulnicky+ (2018) corrected: $\dot M$, M$_\odot$/yr',
    ylabel=r'This paper: $\dot M$, M$_\odot$/yr',
)
ax.set_aspect('equal')
fig.tight_layout()
fig.savefig('K18-mdot-corrected-comparison-R0.pdf')
None
# -

k18tab['Mdot_will'] = ttt['Md_t']
k18tab['Md_Md'] = k18tab['Mdot_will'] / k18tab['Mdot4']
k18tab['LIR'] = ttt['LIR_will']

k18tab['Mdot4', 'Mdot_will', 'Md_Md', 'Lum.', 'R_0', 'U', 'Peak_70', 'LIR'].to_pandas().applymap(np.log10).corr()

k18tab['Mdot', 'Mdot4', 'Mdot_will'].to_pandas().applymap(np.log10).describe()

k18tab['ID''Mdot4', 'Mdot_will', 'Md_Md']

k18tab[mask_upper_limits]

['Sub-Giant' if 'IV' in s else 'Dwarf' if 'V' in s else 'Giant' for s in k18tab['Sp.T._1']]

k18tab['Lum Class'] = ['Sub-Giant' if 'IV' in s else 'Dwarf' if 'V' in s else 'Giant' for s in k18tab['Sp.T._1']]

k18tab['ell,pc'] = k18tab['ell']*k18tab['R_0']/k18tab['R_0,as']
k18tab['ell/R0'] = k18tab['ell']/k18tab['R_0,as']

ldf = k18tab['Mdot4', 'Mdot_will', 'Md_Md', 'Lum.', 'T_eff', 'V_inf_{}',
             'D', 'R_0,as', 'R_0', 'ell,pc', 'ell/R0', 
             'U', 'Peak_70', 'LIR'].to_pandas().applymap(np.log10)
ldf['Lum Class'] = k18tab['Lum Class']
sns.pairplot(ldf, hue='Lum Class', hue_order=['Dwarf', 'Sub-Giant', 'Giant'], 
             palette='husl', markers=["o", "s", "D"],
             diag_kind='kde', diag_kws=dict(bw='silverman'))

ldf.corr()

sns.pairplot(ldf, hue='Lum Class', hue_order=['Dwarf', 'Sub-Giant', 'Giant'], 
             palette='husl', markers=["o", "s", "D"],
             x_vars=['Lum.', 'T_eff', 'D', 'V_inf_{}', 'R_0', 'ell,pc', 'ell/R0', 'U', 'Peak_70', 'LIR'],
             y_vars=['Mdot4', 'Mdot_will', 'Md_Md'],
             plot_kws=dict(alpha=0.7))

# This shows that 
#
# 1. The K18 (corrected) mass loss rates are very well correlated with the bow physical size $R_0$ ($r = 0.86$) and negatively correlated with radiation field $U$ ($r = -0.83$).  It is not clear which of these is the driving factor, since they are very well correlated between themselves ($r = -0.93$).  It would be good if they were correlated with $U$, since there is the prospect of estimating $U$ in a distance-independant manner, using 70/24 flux ratio. They are also equally well negatively correlated with the relative path length $\ell/R_0$ ($r = -0.84$), even though they are weakly positively correlated with $\ell$ itself. 
#
# 2. My mass loss rates are very well correlated with the shell luminosity $L_{\mathrm{IR}}$ ($r = 0.88$), as we already knew.
#
# 3. The ratio Me/K18 is best correlated with angular size of bow ($r = -0.70$), so for poorly resolved objects, I get a higher $\dot M$ than K18, whereas for best resolved objects I get a lower $\dot M$ than K18.  This is probably not due to underestimating the peak brightness for small objects, as one might think, since `Peak_70` has a slight tendency to _fall_ with `R_0,as` ($r = -0.51$). 

# ## Mass loss versus luminosity

# Finally, we do the plots for my values and the corrected K18 ones.

# O star models from Krticka & Kubat (2017)

def Mdot_Krticka18(L):
    return 10**(-5.69 + 1.63*np.log10(L/1e6))

# B star models from Krticka (2014)

Lum_Krt14 = 10**np.array([4.39, 4.20, 3.99, 3.77, 3.54, 3.30, 3.03, 2.75, 2.43])
Mdot_Krt14 = 1.e-11*np.array([210, 160, 96, 39, 7.9, 3.4, 0.91, 0.072, 0.0])

# Weak wind observations from Marcolino (2009)

Lum_M09 = 10**np.array([4.73, 4.74, 4.96, 4.86, 4.79])
Mdot_M09 = 10**np.array([-9.35, -9.22, -8.92, -8.80, -9.22])

MarVink_I = Table.read('data/Mdot_tables/martins-vink-O-I.fits')
MarVink_III = Table.read('data/Mdot_tables/martins-vink-O-III.fits')
MarVink_V = Table.read('data/Mdot_tables/martins-vink-O-V.fits')

# +
fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = 1e4*k18tab['Lum.'], k18tab['Mdot_will']
c = ax.scatter(xx, yy, 
               c=1.e-3*k18tab['T_eff'], cmap='viridis', #vmin=4.0, vmax=6.0, 
               edgecolors='k', alpha=1.0, label='_nolabel_', zorder=100)
fig.colorbar(c, ax=ax).set_label(r'$T_{\mathrm{eff}}$, kK')
for id_, x, y in zip(tt['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
xmin, xmax = 8e3, 2e6
ymin, ymax = 1e-10, 3e-6
xgrid = np.logspace(4.3, np.log10(xmax))
ax.plot(xgrid, Mdot_Krticka18(xgrid), ls='-', color='k', alpha=0.5, lw=5.0, label='Krtička & Kubat (2018)')
ax.plot(Lum_Krt14, Mdot_Krt14, ls=':', color='k', alpha=0.5, lw=5.0, label='Krtička (2014)')

ax.plot(10**MarVink_V['log L'], 10**MarVink_V['Mdot'], ls='-', color='r', label='Vink (2000)')
ax.plot(10**MarVink_III['log L'], 10**MarVink_III['Mdot'], ls='--', color='r', label='_nolabel_')
ax.plot(10**MarVink_I['log L'], 10**MarVink_I['Mdot'], ls=':', color='r', label='_nolabel_')

ax.plot(Lum_M09, Mdot_M09, '+', alpha=1.0, label='Marcolino (2009)')

ax.legend(fontsize='xx-small')
ax.set(
    xscale='log', yscale='log', 
    xlim=[xmin, xmax], ylim=[ymin, ymax],
    xlabel=r'$\log_{10}\ \left[L_* / L_\odot \right]$',
    ylabel=r'This paper: $\dot M$, M$_\odot$/yr',
)
fig.tight_layout()
fig.savefig('Mdot-from-eta-vs-luminosity.pdf')
None

# +
fig, ax = plt.subplots(figsize=(10, 8))
xx, yy = 1e4*k18tab['Lum.'], k18tab['Mdot4']
c = ax.scatter(xx, yy, 
               c=1.e-3*k18tab['T_eff'], cmap='viridis', #vmin=4.0, vmax=6.0, 
               edgecolors='k', alpha=1.0, label='_nolabel_', zorder=100)
fig.colorbar(c, ax=ax).set_label(r'$T_{\mathrm{eff}}$, kK')
for id_, x, y in zip(tt['ID'], xx, yy):
    ax.annotate(
        str(id_), (x, y), fontsize='xx-small',
        xytext=(4,4), textcoords='offset points',
               )
xmin, xmax = 8e3, 2e6
ymin, ymax = 1e-10, 3e-6
xgrid = np.logspace(4.3, np.log10(xmax))
ax.plot(xgrid, Mdot_Krticka18(xgrid), ls='-', color='k', alpha=0.5, lw=5.0, label='Krtička & Kubat (2018)')
ax.plot(Lum_Krt14, Mdot_Krt14, ls=':', color='k', alpha=0.5, lw=5.0, label='Krtička (2014)')

ax.plot(10**MarVink_V['log L'], 10**MarVink_V['Mdot'], ls='-', color='r', label='Vink (2000)')
ax.plot(10**MarVink_III['log L'], 10**MarVink_III['Mdot'], ls='--', color='r', label='_nolabel_')
ax.plot(10**MarVink_I['log L'], 10**MarVink_I['Mdot'], ls=':', color='r', label='_nolabel_')

ax.plot(Lum_M09, Mdot_M09, '+', alpha=1.0, label='Marcolino (2009)')

ax.legend(fontsize='xx-small', loc='lower right')
ax.set(
    xscale='log', yscale='log', 
    xlim=[xmin, xmax], ylim=[ymin, ymax],
    xlabel=r'$\log_{10}\ \left[L_* / L_\odot \right]$',
    ylabel=r'K18 (corrected): $\dot M$, M$_\odot$/yr',
)
fig.tight_layout()
fig.savefig('Mdot-K18-corrected-vs-luminosity.pdf')
None
# -



# ### More details on K18 table

# Next, get a general overview of their tables.

pd.options.display.max_columns = None
k18df = k18tab.to_pandas()
k18df.describe()

k18df.corr()

# A few intersting nuggets from the correlation matrix:
#
# 1. `ell` and `R_0,as` are very well correlated ($r = 0.982$).  See correlation plot below.  This can maybe be used to estimate H/R
# 2. But we should have logged the quantities first, which I have now done - see plot.  *Anyway, leave this for later*

(k18df['ell']/k18df['R_0,as']).describe()

ldf = k18df[['R_0,as', 'ell']].applymap(np.log10)
sns.pairplot(ldf, kind='reg')

g = sns.regplot('R_0,as', 'ell', data=ldf)
g.set(xlim=[0.0, 3.0], ylim=[0.0, 3.0])
g.set_aspect('equal')

# ## Looking at 24/70 ratio

tab05_01_02['F70/F2x'] = tab05_01_02['F70']/tab05_01_02['F24 or 22']
tab05_01_02['U18'] = k18tab['U']

def loggify(x):
    try:
        return np.log10(x)
    except:
        return x


df050102 = tab05_01_02['U', 'U18', 
                       'F8 or 12', 'F24 or 22', 'F70', 'F160', 
                       'F70/F2x', 'T2x/70', 'FIR', 'FIR_will',
                       'Observatory'
                      ].to_pandas().applymap(loggify)

df050102

df050102.corr()

g = sns.pairplot(df050102.drop(2).drop('F160', axis=1).fillna(-1.5))
g.map_offdiag(corrfunc)

g = sns.pairplot(df050102.drop(2).fillna(-1.5), 
                 x_vars=['F24 or 22', 'F70/F2x'],
                 y_vars=['U18', 'FIR_will'],
                 hue='Observatory',
                 kind='reg'
                )
plt.gcf().set_size_inches(10, 10)

g = sns.pairplot(df050102.drop(2).fillna(-1.5), 
                 x_vars=['F24 or 22', 'F70/F2x'],
                 y_vars=['U18', 'FIR_will'],
                )
g.map_offdiag(corrfunc)
plt.gcf().set_size_inches(10, 10)

# **Conclusions**: 
#
# 1. The radiation field $U$ is not as well correlated with `F70/F2x` as we might wish. We have $r = -0.54$, meaning only 30% of the variance in $U$ can be predicted by knowing the flux ratio.  On the other hand, the mean relation is exactly what we expect from the Cloudy models - see `grain-jratios-vs-U.pdf`.  Strangely, the correlation is much better for just the WISE sources, but there are too few points to draw any conclusions.
#
# 2. The total FIR flux is very well correlated with the 24 micron value ($r = 0.97$), better than at 70 ($r = 0.87$) or 8 micron ($r = 0.61$ from pairplot, but this is affected by missing values - really should be $0.86$, see `.corr()` matrix above).





# ### Earlier stuff

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

a = np.linspace(0.0, 1.0, 11)
b = a**2
np.savetxt("test4jane.dat", np.stack((a, b), axis=1), fmt='%.4e')
!cat test4jane.dat
