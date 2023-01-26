# geometricVoFCartesian
The *geometricVoFCartesian* is an extension library for simulating two-phase flows in OpenFOAM. It features conventional and machine-learning methods to estimate the normal vector and curvature of the fluid-fluid interfaces in geometric Volume-of-Fluid (VoF) framework. 

## Feature Overview
Available methods are as follows:
- **Interface normal vector**: Youngs' Method, Central-Columns Differences Method, Mixed Youngs-Central Method
- **Interface curvature**: Height Function Method, Multilayer-Perceptron (MLP) models (only in 2D currently)

The methods are designed for uniform (isotropic) Cartesian grids. In cases where the interface motion occurs in a smaller subset of the domain, the uniform grid can also be restricted to a smaller zone containing the interface. This is currently achieved using the cellSet functionality. The rest of the domain can use generic unstructured grids. 

### Fundamental Classes
- **ijkZone**: Creation and manipulation of the regular-grid zone containing the interface motions.
- **multilayerPerceptron**: Basic functionalities and data structures for an MLP. Implementation is done using exclusively the C++ Standard Library. Thus, no additional package is required to use new machine learning models.  
- **uniformStencil**: Parallel operations on local cell blocks.
- **interfaceForces**: Base class for new curvature models.

### Solver
**interIsoCartFoam**: Extension to ```interIsoFoam``` to feature newly implemented schemes.

## Machine Learning
Machine learning models employ deep MLP architectures to estimate the interfacial curvarture. The code to develop these models can be found in ```tensorflow``` directory. Models are developed as follows:
1. A synthetic dataset composed of circular arcs of varying sizes is generated with:

```bash tensorflow/scripts/genCircles.sh```

2. About hundred Models are trained with mini-batch optimization using the script:

```bash tensorflow/scripts/mlpTrain.sh```

3. The best performing models are selected, and stored in *mlpCurvatureModels*. Two different types of models are developed: 
- **SymMLP**: the symmetry-preserving MLP model using bias-free neurons
- **StdMLP**: the standard MLP model
4. The weights and biases of the TensorFlow models are converted to raw text format, e.g.,

``` python tensorflow/TfToTxt.py mlpCurvatureModels SymMLP```

The ```multilayerPerceptron``` class in OpenFOAM reads these parameters and constructs the corresponding MLP model.

The details of MLP models can be found in:

Önder, A., & Liu, P. L.-F. (2022). Deep learning of interfacial curvature: a symmetry-preserving approach for the volume of fluid method. arXiv. http://arxiv.org/abs/2206.06041


## Prerequisites
OpenFOAM v2006 must be installed:

https://www.openfoam.com/news/main-news/openfoam-v20-06
## Installation
```
git clone https://github.com/asimonder/geometricVoFCartesian.git
cd geometricVoFCartesian
./Allwclean
./Allwmake
```

## Usage
The library exclusively works on Cartesian zones. Class ```ijkZone``` handles the regular-grid functionality. It requires the cellset where the interface is located to be specified. If cyclic BCs are considered, they need to defined as well. The dictionary is located in ```system/ijkDict```:

```
{
cellSet refinedCells;
cyclicX true;
}
```
If the grid is uniform everywhere in the domain then simply use "```cellSet domain;```".

## Examples 
Several 2D benchmark cases are provided.

## Author
Asim Önder (asim.onder[at]gmail.com)
