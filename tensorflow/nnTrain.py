# mlp for regression
from numpy import sqrt
#import pandas as pd
#from pandas import read_csv
import numpy as np
import tensorflow as tf
from tensorflow import keras
#from sklearn.model_selection import train_test_split
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense
#from tensorflow.keras.layers import Dropout
#from tensorflow.keras.constraints import maxnorm
from tensorflow.keras.models import load_model
#from keras.optimizers import SGD
import os
import sys
import pickle
from sklearn.utils import shuffle

dataFile=sys.argv[1]
NData=int(sys.argv[3]) 
doNormal=int(sys.argv[4]) 
doCurv=int(sys.argv[5]) 
NEpocs=int(sys.argv[6]) 
actFun=sys.argv[7]
actFun2=sys.argv[8]
Nbatch=int(sys.argv[9])
mom=float(sys.argv[10])
xStencil=int(sys.argv[11])
yStencil=int(sys.argv[12])
doSGD=int(sys.argv[13])
continueModel=int(sys.argv[14])
NInt=int(sys.argv[15])
#scaleInput=int(sys.argv[16])
lr=float(sys.argv[16])

scaleInput=False

NL=sys.argv[2].split(',')
NL=[int(i) for i in NL]
NL.append(1)
NLayers=len(NL)
funcs=[actFun]
for i in range(1,NLayers-1):
    funcs.append(actFun)
funcs.append(actFun2)

initial_learning_rate = lr
lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
    initial_learning_rate,
    decay_steps=2000,
    decay_rate=0.99,
    staircase=False)


optim=tf.keras.optimizers.Adam(
    learning_rate=lr_schedule,
    beta_1=0.9,
    beta_2=0.999,
    epsilon=1e-8,
    amsgrad=False,
    name="adam")
if doSGD or continueModel>-1:
    print("Training of the model will be done with the SGD method.")
    optim=tf.keras.optimizers.SGD(
        learning_rate=0.00001, momentum=0.90, nesterov=False, name='SGD')

#mldir='/home/projects/11001025/machineLearning/trainElementaryShapes'
mldir='.'

if doSGD:
    fName="%s_%dx%d_sgd%d_mom%.1f"%(dataFile,xStencil,yStencil,Nbatch,mom)
else:
    fName="%s_%dx%d_adam%d"%(dataFile,xStencil,yStencil,Nbatch)
if scaleInput:
    fName="scaledInput_%s"%fName
mDir="%s_NData%dK_epocs%d"%(fName,int(NData/1000),NEpocs)
wDir='%s/models_NInt%d/%s'%(mldir,NInt,mDir)
    
iC=1
for i in range(0,NLayers):    
    if funcs[i]==funcs[i-1] and NL[i]==NL[i-1]:
        iC+=1
    else:
        if i>0:
            wDir=wDir+'x%d'%(iC)
        iC=1
        wDir=wDir+'_%s%d'%(funcs[i],NL[i])
wDir=wDir+'x%d'%(iC)

nMax=2;NS=5;
iMax=int((xStencil-1)/2)
jMax=int((yStencil-1)/2)
indices=np.zeros(xStencil*yStencil,dtype=int)
iC=0
for i in range(-iMax,iMax+1,1):
    for j in range(-jMax,jMax+1,1):
        indices[iC]=(i+nMax)+NS*(j+nMax)
        iC+=1

use_bias=True
        
# load the dataset
path = "datasets_NInt%d/%s"%(NInt,dataFile)
file = open(path,'rb')
data=pickle.load(file)
print("Size of the dataset: ",data['kappa'].shape[0])
file.close()
x_train = data['alpha'][:NData,indices]
if scaleInput:
    x_train = 2.0*(x_train-0.5)
n_features = x_train.shape[1]
print("x_train.shape=",x_train.shape)

if continueModel>-1:
    Nbatch=32
    NEpocs=1000
    print("Training continuation is done with batch size Nbatch=%d and NEpocs=%d"%(Nbatch,NEpocs))

#first-order derivatives
if doNormal:
    y_train = data['nx']/data['ny']
    NL[-1]=1
    if continueModel>-1:
        wDir="%s_%d"%(wDir,continueModel)
        modelPath="%s/model_Curv.h5"%(wDir)
        if os.path.exists(modelPath):
            modelG=load_model(modelPath)
        else:
            sys.exit('Missing model at %s'%modelPath)
    else:
        modelG = Sequential()
        if funcs[0]=='lrelu':
            modelG.add(Dense(NL[0],activation=tf.keras.layers.LeakyReLU(alpha=0.01),input_shape=(n_features,),use_bias=use_bias))
        else:
            modelG.add(Dense(NL[0], activation=funcs[0], input_shape=(n_features,),use_bias=use_bias))
        
            for i in range(1,NLayers):
                if funcs[i]=='lrelu':
                    modelG.add(Dense(NL[i],activation=tf.keras.layers.LeakyReLU(alpha=0.01),use_bias=use_bias))
                else:
                    modelG.add(Dense(NL[i], activation=funcs[i],use_bias=use_bias))

    modelG.summary()
    modelG.compile(optimizer=optim,loss='mse', metrics=['mae', 'mape'])
    #history =modelG.fit(x_train, y_train, epochs=NEpocs, batch_size=Nbatch, verbose=1, validation_data=(X_val,y_val))
    history =modelG.fit(x_train, y_train, epochs=NEpocs, batch_size=Nbatch, verbose=2, validation_split=0.3)
    for iF in range(100):
        wDir2="%s_%d"%(wDir,iF)
        fName="%s/model_Der1.h5"%(wDir2)
        if not os.path.exists(fName):
            wDir=wDir2
            try:
                os.makedirs(wDir)
            except:
                "Directory %s exists."%wDir
            break

    modelG.save('%s/model_Der1.h5'%wDir)
    f = open("%s/history_loss_Der1"%wDir, "a")
    np.savetxt(f,history.history['loss'])
    f.close()
    f = open("%s/history_val_loss_Der1"%wDir, "a")
    np.savetxt(f,history.history['val_loss'])
    f.close()

    #error = modelG.evaluate(X_test, y_testG, verbose=1)
    #print('test of first derivatives: MSE= %.12f, RMSE= %.12f' % (error, sqrt(error)))

#curvature
if doCurv:
    y_train = data['kappa'][:NData]
    NL[-1]=1
    modelC = Sequential()
    if continueModel>-1:
        wDir="%s_%d"%(wDir,continueModel)
        modelPath="%s/model_Curv.h5"%(wDir)
        if os.path.exists(modelPath):
            modelC=load_model(modelPath)
        else:
            sys.exit('Missing model at %s'%modelPath)
    else:
        if funcs[0]=='lrelu':
            modelC.add(Dense(NL[0],activation=tf.keras.layers.LeakyReLU(alpha=0.01),use_bias=use_bias,input_shape=(n_features,)))
        else:
            modelC.add(Dense(NL[0], activation=funcs[0], input_shape=(n_features,),use_bias=use_bias))
        
        for i in range(1,NLayers):
            if funcs[i]=='lrelu':
                modelC.add(Dense(NL[i],activation=tf.keras.layers.LeakyReLU(alpha=0.01),use_bias=use_bias))
            else:
                #if dropOut>0 and i<NLayers-1:
                #modelC.add(Dropout(dropOut))
                if i==NLayers-1:
                    modelC.add(Dense(NL[i], activation=funcs[i],use_bias=True))
                else:
                    modelC.add(Dense(NL[i], activation=funcs[i],use_bias=use_bias))
                           
    modelC.summary()
    modelC.compile(optimizer=optim,loss='mse', metrics=['mae', 'mape'])
    X, y = shuffle(x_train, y_train, random_state=0)
    history =modelC.fit(X, y, epochs=NEpocs,batch_size=Nbatch, verbose=2, validation_split=0.3)
    #history =modelC.fit(x_train, y_train, epochs=NEpocs,batch_size=Nbatch, verbose=2, validation_split=0.3)
    for iF in range(1000):
        wDir2="%s_%d"%(wDir,iF)
        fName="%s/model_Curv.h5"%(wDir2)
        if not os.path.exists(fName):
            wDir=wDir2
            try:
                os.makedirs(wDir)
            except:
                "Directory %s exists."%wDir
            break

    modelC.save('%s/model_Curv.h5'%wDir)
    f = open("%s/history_loss_Curv"%wDir, "a")
    np.savetxt(f,history.history['loss'])
    f.close()
    f = open("%s/history_val_loss_Curv"%wDir, "a")
    np.savetxt(f,history.history['val_loss'])
    f.close()
    f = open("%s/history_mae"%wDir, "a")
    np.savetxt(f,history.history['mae'])
    f.close()
    f = open("%s/history_val_mae"%wDir, "a")
    np.savetxt(f,history.history['val_mae'])
    f.close()
    f = open("%s/history_mape"%wDir, "a")
    np.savetxt(f,history.history['mape'])
    f.close()
    f = open("%s/history_val_mape"%wDir, "a")
    np.savetxt(f,history.history['val_mape'])
    f.close()

    #evaluate the model
    #error = modelC.evaluate(X_test, y_testC, verbose=1)
    #print('test of curvatures: MSE= %.12f, RMSE= %.12f' % (error, sqrt(error)))

    #NIntM=200
    #dataFile="waveLong_lambMin_8.0_lambMax_2048.0_HsMin_0.00_HsMax_0.26_psi_-0.000.pkl"
    #dataFile="circleExp_R_4.0_4096.0_0.pkl"

