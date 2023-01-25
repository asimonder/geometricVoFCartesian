import pandas as pd
import os
import sys
import numpy as np
from sinwave_prosperetti import sinwave_prosperetti

model=sys.argv[1]
curv=sys.argv[2]
norm=sys.argv[3]
Co=sys.argv[4]
#curvScheme=sys.argv[2]

print("model=",model)
print("curv=",curv)
print("norm=",norm)

rho1 = 1  # Density liquid 1
rho2 = 1  # Density liquid 2
wavelength = 0.003  # Wavelength
H0 = 3e-5  # Initial height
nu = 0.001  # Kinematic viscosity
sigma = 1       # Surface tension
ana = sinwave_prosperetti(rho1, rho2, wavelength, H0, nu, sigma)
#analytical = ana.tabulatedData(25, 51)
NAna=5001
analytical = ana.tabulatedData(25, NAna)
analytical.columns = ['time', "analytical"]
analytical["analytical"] = abs(analytical["analytical"])
print(analytical["time"].shape,analytical["analytical"].shape)

ana=np.zeros((NAna,2))
ana[:,0]=analytical["time"]
ana[:,1]=analytical["analytical"]

data=pd.read_csv("%s_%s_%s_%s.csv"%(model,curv,norm,Co))
#data=pd.read_csv("%s"%())
#data=pd.read_csv("sinWaveHex_%s_%s.csv"%(norm,model))
data=data.sort_values('time')
NL=np.array([8,16,32,64,128])
errL2=np.zeros(NL.shape[0])
errLMax=np.zeros(NL.shape[0])
iC=0
for nCells in NL:
    #if data[data['interfaceType'] == 'mlpNormal'].values.shape[0]>0:
    #    norm='mlpNormal'
    #if data[data['Method'] == 'mlpCurvature'].values.shape[0]>0:
    #    curv='mlpCurvature'
    #if data[data['Method'] == 'mlpCurvatureNormal'].values.shape[0]>0:
    #    curv='mlpCurvatureNormal'
    #else:
    #    curv=model

    times=data["time"][data['interfaceType'] == norm]
    times=times[data['Method'] == curv]
    times=times[data['nCells'] == nCells]
    etaMax=data["max"][data['interfaceType'] == norm]
    etaMax=etaMax[data['Method'] == curv]
    etaMax=etaMax[data['nCells'] == nCells]

    #print(nCells,etaMax)

    if data[data['interfaceType'] == norm].values.shape[0]<0:
       sys.exit("No data available for %s!"%norm);
    if data[data['Method'] == curv].values.shape[0]<0:
       sys.exit("No data available for %s!"%curv);
    
    #etaRef=np.interp(times, ana[:,0], ana[:,1])
    #etaRef/=etaRef[0]
    etaSim=np.interp(ana[:,0],times, etaMax.values)
    etaAna=ana[:,1]/ana[0,1]
    N=np.float(etaSim.shape[0])
    errL2[iC]=np.linalg.norm(etaAna-etaSim,ord=2)/np.sqrt(N)
    errLMax[iC]=np.linalg.norm(etaAna-etaSim,ord=np.inf)
    iC+=1

df = pd.DataFrame({'name': ["Co%s_%s_%s"%(Co,curv,model)],
                   'recondScheme': [norm],
                   "%d"%NL[0]: [0.],
                   "%d"%NL[1]: [0.],
                   "%d"%NL[2]: [0.],
                   "%d"%NL[3]: [0.],
                   "%d"%NL[4]: [0.]})

df.iloc[0,2:]=errL2[0:]
fName="resultsL2.dat"
if os.path.exists(fName):
    df2=pd.read_csv(fName)
    df2=df2.append(df,ignore_index=True)
    df2=df2.sort_values(by=['%s'%NL[-2]])
    df2.to_csv(path_or_buf=fName, index=False)    
else:
    df.to_csv(path_or_buf=fName, index=False)
df.iloc[0,2:]=errLMax[0:]
fName="resultsLMax.dat"
if os.path.exists(fName):
    df2=pd.read_csv(fName)
    df2=df2.append(df,ignore_index=True)
    df2=df2.sort_values(by=['%s'%NL[-2]])
    df2.to_csv(path_or_buf=fName, index=False)    
else:
    df.to_csv(path_or_buf=fName, index=False)

    
