# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import casefoam
from casefoam import postFunctions
import seaborn as sns
from sinwave_prosperetti import sinwave_prosperetti

import os
import sys

cases = [['mycCartesian'], 
        ['mlpCurvature'],
        ['coarse','coarser','coarsest','mid','fine']] 

baseCase = 'Cases' 
solutionDir = 'surfaces'
file = 'alpha.water_constantIso.raw'

rho1 = 1  # Density liquid 1
rho2 = 1  # Density liquid 2
wavelength = 0.003  # Wavelength
H0 = 3e-5  # Initial height
nu = 0.001  # Kinematic viscosity
sigma = 1       # Surface tension
ana = sinwave_prosperetti(rho1, rho2, wavelength, H0, nu, sigma)
#analytical = ana.tabulatedData(25, 51)
analytical = ana.tabulatedData(25, 51)
analytical.columns = ['time', "analytical"]
analytical["analytical"] = abs(analytical["analytical"])


postFunction = postFunctions.getFreeSurfaceWallAndCentre

sol = casefoam.posField_to_timeSeries(
    solutionDir, file, postFunction, cases, baseCase, axis=2)
sol = sol.reset_index()
sol.columns = ['time', 'min', 'mean', 'max',
               'interfaceType', 'Method', 'nCells']
sol = sol.replace('coarsest', 8)
sol = sol.replace('coarser', 16)
sol = sol.replace('coarse', 32)
sol = sol.replace('mid', 64)
sol = sol.replace('fine', 128)
sol['max'] /= H0
sol['time'] /= 1.475e-5

for ind in cases[0]:
    print("converting csv for %s"%ind)
    results = sol[sol['interfaceType'] == ind] 
    results.to_csv("sinWaveHex_%s.csv"%ind, index=False)
    
if 0:
    ax = analytical.plot(style='.', x='time', y='analytical',
                     c='black', marker='+', ms=7)
    sns.set_style("ticks")
    plicRDF = sol[sol['interfaceType'] == 'plicRDF'] #sol[sol['interfaceType'] == 'plicRDF']
    ax = sns.lineplot(x='time', y='max', hue='Method', style="nCells", data=plicRDF,
                      hue_order=['gradAlpha', 'heightFunctionStructured2D', 'RDF'], ax=ax)
    plt.ylabel('Relative amplitude')
    plt.xlabel('Non-dimensional time')
    plt.savefig("sinWaveHex_plicRDF.pdf")
    #plt.savefig("sinWaveHex_youngCartesian.pdf")
    #plicRDF.to_csv("sinWaveHex_plicRDF.csv", index=False)
    #plicRDF.to_csv("sinWaveHex_youngCartesian.csv", index=False)
    plt.show()

