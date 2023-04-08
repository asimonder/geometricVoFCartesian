#!/bin/bash -l   

#RMin=2;expR=9;
#RMin=4;expR=9;
#RMin=4;expR=8;
#RMin=2;expR=11;
RMin=2;expR=12;
#RMin=4;expR=12;
NInt=200
symModel=1
dist=Exp
#dist=UniR
#dist=UniK
evalPoint=Centroid
rndOrg=1
invK=1

for t in {0..10}
do
    echo $t
    python genCircles.py 10000 $RMin $expR $NInt $dist $symModel $evalPoint $rndOrg $invK
done
