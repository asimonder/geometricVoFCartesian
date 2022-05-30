#!/bin/bash -l   
#PBS -N circle
#PBS -r n
#PBS -j oe
#PBS -m abe
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1:mem=96gb
#PBS -P 11001025
# PBS -q largemem
 
#RMin=2;expR=9;
#RMin=4;expR=9;
#RMin=4;expR=8;
RMin=2;expR=11;
#RMin=4;expR=12;
NInt=200
symModel=1
dist=Exp
#dist=UniR
#dist=UniK
evalPoint=Centroid
rndOrg=1
invK=0

for t in {0..100}
do
    echo $t
    python genCircles.py 1000 $RMin $expR $NInt $dist $symModel $evalPoint $rndOrg $invK
done

