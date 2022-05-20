import shapes as shp
from shapes import circle
from shapes import stencil
import os
import sys
import pickle
import numpy as np
import time

NData=int(sys.argv[1])
RMin=float(sys.argv[2])
expR=float(sys.argv[3])
NInt=int(sys.argv[4])
dist=sys.argv[5]
symModel=int(sys.argv[6])
evalPoint=sys.argv[7]
rndOrigin=int(sys.argv[8])
invertCurvature=int(sys.argv[9])
start=time.time()

RMax=RMin*2**expR
data={}
data["theta"]=np.zeros(NData)
data["alpha"]=np.zeros((NData,25))
data["kappa"]=np.zeros(NData)
data["nx"]=np.zeros(NData)
data["ny"]=np.zeros(NData)
for i in range(0,NData):
    if dist=="Exp":
        R=RMin*2**np.random.uniform(0,expR)
    elif dist=="UniR":
        R=np.random.uniform(RMin,RMax)
    elif dist=="UniK":
        K=np.random.uniform(1./RMin,1./RMax)
        R=1./K
    else:
        exit("distribution is not known!")
    theta=np.random.uniform(0,2*np.pi)
    if rndOrigin:
        ox=np.random.uniform(-1,1)
        oy=np.random.uniform(-1,1)
        origin=np.array([ox,oy])
    else:
        origin=np.zeros(2)
    circ=circle(R,origin)
    sten=stencil(theta,circ,lenX=5,lenY=5,N=NInt,evalPoint=evalPoint)
    sgnA=np.random.choice([-1, 1.])
    if sgnA<0:
        sten.alpha=1.-sten.alpha
        sten.K*=-1.
        sten.n*=-1.
    if symModel:
        sten.convertZone(invertCurvature)
    #else:
            
    data["theta"][i]=theta
    data["alpha"][i,:]=sten.alpha.flatten(order='F')
    data["nx"][i]=sten.n[0]
    data["ny"][i]=sten.n[1]
    data["kappa"][i]=sten.K
    print(i)
    
wDir="datasets_NInt%d"%NInt
if rndOrigin:
    name="circleOrg"
else:
    name="circle"

name="%s%s"%(name,evalPoint)
if symModel:
    if invertCurvature:
        fileName="%s%s_R_%.1f_%.1f"%(name,dist,RMin,RMax)
    else:
        fileName="%sAllK%s_R_%.1f_%.1f"%(name,dist,RMin,RMax)
else:
    fileName="%sRaw%s_R_%.1f_%.1f"%(name,dist,RMin,RMax)
shp.writeDataset(wDir,fileName,data)

end=time.time()
print("elapsed time %f seconds."%(end-start))
