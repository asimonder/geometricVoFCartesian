#!/bin/bash -l   
#PBS -N ascw
#PBS -r n
#PBS -j oe
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=12:mem=48gb
#PBS -P 11001025
# PBS -q largemem
 
cd $PBS_O_WORKDIR


modelType=normal_circleOrgCentroidAllNormalExp_R_2.0_4096.0.pkl_3x3_adam512_NData500K_epocs2000_tanh30x4_linear1x1_99_v3_1e3_ext
modelType=normal_circleOrgCentroidAllNormalExp_R_2.0_4096.0.pkl_3x3_adam512_NData500K_epocs2000_tanh30x4_linear1x1_99_v4_1e3_emp
modelType=circleOrgCentroidRawExp_R_2.0_4096.0.pkl_3x3_adam512_NData500K_epocs2000_tanh30x4_linear1x1_7_v4_1e3_emp
modelType=QiJCP2019_v4_1e3_emp

#normal=mlpNormal
curv=mlpCurvature
normal=mycCartesian
Co=0.0005

[ -f "$FILE" ]
for d in $models/"$modelType"*;
do
    curvModel=$(basename $d)
    echo $curvModel
    if grep -RF "$curvModel," resultsL2.dat 
    then
	continue;
    else
	runDir="$curvModel"_"$curv"_"$normal"_"$Co"
	fName="$runDir".csv

	sed "s/\['mycCartesian'\]/\['"$normal"'\]/g" genCases.py>genCases_"$runDir".py
	sed -i "s/\['mlpCurvature'\]/\['"$curv"'\]/g" genCases_"$runDir".py
	sed "s/\['mycCartesian'\]/\['"$normal"'\]/g" getResults.py>getResults_"$runDir".py
	sed -i "s/\['mlpCurvature'\]/\['"$curv"'\]/g" getResults_"$runDir".py

	rm -r $runDir
	mkdir $runDir
	cp -r sinWave $runDir/
	cd $runDir
	rm -r sinWave/machineLearningModels
	cp -r $models/$curvModel sinWave/machineLearningModels
	#cp $models/$normalModel/*der1* sinWave/machineLearningModels/
	sed -i 's/zonalModel false/zonalModel false/g' sinWave/constant/transportProperties
	#sed -i 's/zonalModel false/zonalModel true/g' sinWave/constant/transportProperties
	sed -i 's/useScaling true/useScaling true/g' sinWave/constant/transportProperties
	sed -i 's/iMax 2/iMax 1/g' sinWave/constant/transportProperties
	sed -i 's/jMax 1/jMax 1/g' sinWave/constant/transportProperties
	sed -i 's/nMax 2/nMax 1/g' sinWave/constant/transportProperties
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
	#rm -r $runDir
	rm genCases_"$runDir".py
	rm getResults_"$runDir".py
    fi
done
