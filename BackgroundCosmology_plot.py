# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 14:54:18 2024
BackgroundCosmology_plot
1.0 To plot data generated from cosmology _cpp
@author: csyang
"""

import re
import os
import pandas as pd
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

os.chdir(r'C:\Users\csyang\source\repos\cosmology_cpp\cosmology_cpp')
mpl.rcParams['figure.dpi'] = 600

#%% Background Cosmology 
file_path = 'cosmology.txt'
column_names = ["x","eta","t","Hp","dHpdx","ddHpddx","OmegaB","OmegaCDM","OmegaLambda","OmegaR","OmegaNu","OmegaK"]

# Read the file into a DataFrame
cosmo_df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=column_names)

# Display the first few rows of the DataFrame
print(cosmo_df)

# Unit conversion
# eta: m -> Mpc
cosmo_df["eta_Mpc"] = cosmo_df["eta"]/3.08567758e22
cosmo_df["t_Gy"] = cosmo_df["t"]/((1e9*365*24*3600))
cosmo_df["Hp_h"] = cosmo_df["Hp"]/(100e3/3.08567758e22)

# Calculating eta*Hp/c
cosmo_df["etaHpc"] = cosmo_df["eta"]*cosmo_df["Hp"]/2.99792458e8

# Calculating combined density parameter
cosmo_df["Omega_rel"] = cosmo_df["OmegaR"] + cosmo_df["OmegaNu"]
cosmo_df["Omega_mat"] = cosmo_df["OmegaB"] + cosmo_df["OmegaCDM"]

# Calculating (1/Hp)*(dH/dx)
cosmo_df["der1divHp"] = cosmo_df["dHpdx"] / cosmo_df["Hp"]

# Calculating (1/Hp)*(d2H/dx2)
cosmo_df["der2divHp"] = cosmo_df["ddHpddx"] / cosmo_df["Hp"]

#%% Plots
# Eta against x
plt.plot(cosmo_df["x"],cosmo_df["eta_Mpc"])
plt.yscale('log')
plt.xlabel(r"$x$")
plt.ylabel(r"$\eta(x)$ (Mpc)")
plt.title(r"Plot of Conformal time, $\eta(x)$ against $x$")

# t against x
plt.plot(cosmo_df["x"],cosmo_df["t_Gy"])
plt.yscale('log')
plt.xlabel(r"$x$")
plt.ylabel(r"$t(x)$ (Gy)")
plt.title(r"Plot of Time, $t(x)$ (Gy) against $x$")

# Hp against x
plt.plot(cosmo_df["x"],cosmo_df["Hp_h"])
plt.yscale('log')
plt.xlabel(r"$x$")
plt.ylabel(r"$\mathcal{H}(x)$ (Mpc)")
plt.title(r"Plot of Scaled hubble parameter, $\mathcal{H}(x)$ against $x$")

# Eta*Hp/c against x
plt.plot(cosmo_df["x"],cosmo_df["etaHpc"])
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{\eta(x)\mathcal{H}(x)}{c}$")
plt.xlim([-15.0,0.0])
plt.ylim([0.5,3.0])
plt.title(r"Plot of $\frac{\eta(x)\mathcal{H}(x)}{c}$ against $x$")

# Hp and (dHpdx/2) against x
plt.plot(cosmo_df["x"],cosmo_df["Hp"],label=r"$\mathcal{H}$")
plt.plot(cosmo_df["x"],cosmo_df["dHpdx"]/2,label=r"$\frac{d\mathcal{H}}{dx}$")
plt.xlabel(r"$x$")
plt.legend()
plt.title(r"Plot of $\mathcal{H}$ and $\frac{1}{2}\frac{d\mathcal{H}}{dx}$ against $x$")

# Component density parameter against x
plt.plot(cosmo_df["x"],cosmo_df[["Omega_rel"]],label=r"$\Omega_{\rm relativistic} = \Omega_r + \Omega_{\nu}$")
plt.plot(cosmo_df["x"],cosmo_df[["Omega_mat"]],label=r"$\Omega_{\rm matter} = \Omega_b + \Omega_{CDM}$")
plt.plot(cosmo_df["x"],cosmo_df[["OmegaLambda"]],label=r"$\Omega_{\Lambda}$")
plt.xlabel(r"$x$")
plt.ylabel(r"$\Omega_i(x)$")
plt.legend()
plt.title(r"Plot of density parameter, $\Omega_i(x)$ against $x$")

# Checking H and its derivatives
# First derivative
plt.plot(cosmo_df["x"],cosmo_df["der1divHp"])
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$")
plt.title(r"Plot of $\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$ against $x$")

# Second derivative
plt.plot(cosmo_df["x"],cosmo_df["der2divHp"])
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$")
plt.title(r"Plot of $\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$ against $x$")

#%% Supernova fitting
filepathSpfit = "results_supernovafitting.txt"
columnNamesSpfit = ["Chi2","h","OmegaM","OmegaK"]

# Read the file into a DataFrame
spFit_df = pd.read_csv(filepathSpfit, delim_whitespace=True,skiprows=200,names=columnNamesSpfit)
spFit_df['OmegaLambda'] = 1 - spFit_df['OmegaM'] -spFit_df['OmegaK'] 
minRow = spFit_df.iloc[spFit_df.Chi2.idxmin()]
minChi2 = minRow.Chi2
spFit_1sigma = spFit_df.loc[spFit_df["Chi2"] < minChi2 + 3.53]
spFit_1sigma.plot.scatter("OmegaM","OmegaLambda")




