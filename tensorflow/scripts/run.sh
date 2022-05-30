#!/bin/bash -l   

conda activate tf


dataFile=circleOrgCentroidAllNormalExp_R_2.0_4096.0.pkl
#dataFile=circleOrgCentroidRawExp_R_2.0_4096.0.pkl


NData=500000
doNormal=0
doCurv=1
#NEpocs=1500
#NEpocs=5000
NEpocs=2000
#actFun=relu
actFun=tanh
actFun2=linear
NBatch=512
mom=0.9
#NLayers=100,100
#NLayers=15,15,15,15
#NLayers=30,30,30,30
NLayers=30,30,30,30
#NLayers=60,60,60,60
#NLayers=320,320
#NLayers=200,200,200
#NLayers=120,120,120,120
#NLayers=8,8,8,8
#NLayers=15,15
#NLayers=15,15,15
#NLayers=120
doSGD=0
xStencil=3;yStencil=3;
#xStencil=5;yStencil=3;
#xStencil=3;yStencil=5;
#xStencil=5;yStencil=5;
#continueModel=0 #6
continueModel=-1
NInt=200

lr=0.002
use_bias=0

python nnTrain.py $dataFile $NLayers $NData $doNormal $doCurv $NEpocs $actFun $actFun2 $NBatch $mom $xStencil $yStencil $doSGD $continueModel $NInt $lr $use_bias
