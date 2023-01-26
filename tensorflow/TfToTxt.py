# mlp for regression
import numpy as np
import os
import sys
import ast
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
modelDir=sys.argv[1]       
modelName=sys.argv[2]

#if curv:
fName="%s/%s/curv_activations.txt"%(modelDir,modelName)
if os.path.exists(fName):
    print("Txt files for the curvature model exists already.")
    exit();
modelPath="%s/%s/model_Curv.h5"%(modelDir,modelName)
if not os.path.exists(modelPath):
    print("TF curvature model is not available.")
    exit();
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import load_model
modelCurv=load_model(modelPath)
file1 = open(fName, "w")
model_config = modelCurv.get_config()
iL=0
for layer in model_config['layers']:
    name = layer['class_name']
    if name!="InputLayer":
        iL+=1
        file1.write("%s\n"%layer['config']['activation'])
NLayers=iL
file1.close()
NWeights=len(modelCurv.weights)
for i in range(NLayers):
    if NLayers==NWeights:
        data=modelCurv.weights[i].numpy()
        np.savetxt("%s/%s/curv_weights%d.txt"%(modelDir,modelName,i),data.T,delimiter=',')
        data=np.zeros(modelCurv.weights[i].numpy().shape[1])
        np.savetxt("%s/%s/curv_biases%d.txt"%(modelDir,modelName,i),data,delimiter=',')
    elif int(2*NLayers)==NWeights:
        data=modelCurv.weights[2*i].numpy()
        np.savetxt("%s/%s/curv_weights%d.txt"%(modelDir,modelName,i),data.T,delimiter=',')
        data=modelCurv.weights[2*i+1].numpy()
        np.savetxt("%s/%s/curv_biases%d.txt"%(modelDir,modelName,i),data,delimiter=',')            
    elif int(NLayers+1)==NWeights:
        if i<NLayers-1:
            data=modelCurv.weights[i].numpy()
            np.savetxt("%s/%s/curv_weights%d.txt"%(modelDir,modelName,i),data.T,delimiter=',')
            data=np.zeros(modelCurv.weights[i].numpy().shape[1])
            np.savetxt("%s/%s/curv_biases%d.txt"%(modelDir,modelName,i),data,delimiter=',')
        else:
            data=modelCurv.weights[i].numpy()
            np.savetxt("%s/%s/curv_weights%d.txt"%(modelDir,modelName,i),data.T,delimiter=',')
            data=modelCurv.weights[i+1].numpy()
            np.savetxt("%s/%s/curv_biases%d.txt"%(modelDir,modelName,i),data,delimiter=',')

