
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
Machine learning models employ deep MLP architectures to estimate the interfacial curvarture. The code to develop these models can be found in tensorflow directory. Models are developed in \
three steps:
1. A synthetic dataset composed of circular arcs of varying sizes is generated with *tensorflow/scripts/genCircles.sh*.
2. Models are trained with mini-batch optimization using the script *tensorflow/scripts/mlpTrain.sh*.
3. The best performing models are selected, and the weights and biases of these Tensorflows models are converted to standard ascii format. The *multilayerPerceptron* class in OpenFOAM reads these parameters and constructs the corresponding MLP model.

The details of MLP models can be found in:

Önder, A., & Liu, P. L.-F. (2022). Deep learning of interfacial curvature: a symmetry-preserving approach for the volume of fluid method. arXiv. http://arxiv.org/abs/2206.06041


## Getting Started
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

## Examples 
Several 2D benchmark cases are provided.

## Author
Asim Önder



