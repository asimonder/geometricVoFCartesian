#!/bin/bash -l   
#PBS -N hf
#PBS -r n
#PBS -j oe
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=12:mem=48gb
#PBS -P 11001025
# PBS -q largemem
 
cd $PBS_O_WORKDIR


modelType=conv

curv=heightFunction #mlpCurvature
normal=mycCartesian
Co=0.0005

curvModel=$curv #$(basename $d)
echo $curvModel
runDir="$curvModel"_"$curv"_"$normal"_"$Co"
fName="$runDir".csv

sed "s/\['mycCartesian'\]/\['"$normal"'\]/g" genCases.py>genCases_"$runDir".py
sed -i "s/\['mlpCurvature'\]/\['"$curv"'\]/g" genCases_"$runDir".py
#sed -i "s/sinWave/sinWave/g" genCases_"$runDir".py
sed "s/\['mycCartesian'\]/\['"$normal"'\]/g" getResults.py>getResults_"$runDir".py
sed -i "s/\['mlpCurvature'\]/\['"$curv"'\]/g" getResults_"$runDir".py

rm -r $runDir
mkdir $runDir
cp -r sinWave $runDir/
cd $runDir
#cp $models/$normalModel/*der1* sinWave/machineLearningModels/
sed -i 's/interfaceTol 1e-4/interfaceTol 1e-3/g' sinWave/constant/transportProperties
sed -i 's/fillNeighbours 0/fillNeighbours 4/g' sinWave/constant/transportProperties
sed -i "s/0.0005/"$Co"/g" sinWave/system/controlDict

python ../genCases_"$runDir".py
sed -i 's/&/ /g' Allrun
sed -i 's/fine\/Allrun/fine\/AllrunParallel/g' Allrun
./Allrun
python ../getResults_"$runDir".py
	
mv sinWaveHex_"$normal".csv ../$fName
cd ..
python calcErrors.py $curvModel $curv $normal $Co
rm -r $runDir
#rm genCases_"$runDir".py
#rm getResults_"$runDir".py
