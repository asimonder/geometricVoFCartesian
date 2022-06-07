
# geometricVoFCartesian
The *geometricVoFCartesian* is an extension library for simulating two-phase flows in OpenFOAM. It features conventional and machine-learning methods to estimate the normal vector and curvature of the fluid-fluid interfaces in geometric Volume-of-Fluid (VOF) framework. 

## General Description
The methods are designed for uniform (isotropic) Cartesian grids. In cases where the interface motion occurs in a smaller subset of the domain, the uniform grid can also be restricted to a smaller zone containing the interface. This is currently achieved using the cellSet functionality. The rest of the domain can use generic unstructured grids. 

Available methods are as follows:
- **Interface normal vector**: Youngs' Method, Central-Columns Differences Method, Mixed Youngs-Central (MYC) Method
- **Interface curvature**: Height Function Method, Multilayer-Perceptron (MLP) models (only in 2D currently)

### Fundamental Classes
- **ijkZone**: Creation and manipulation of the regular-grid zone containing the interface motions.
- **multilayerPerceptron**: Basic functionalities and data structures for an MLP. Implementation is done using exclusively the C++ Standard Library. Thus, no additional package is required to use new machine learning models.  
- **uniformStencil**: Parallel operations on local cell blocks.
- **interfaceForces**: Base class for new curvature models.


### Solver
**interIsoCartFoam**: Extension to ```interIsoFoam``` to feature newly implemented schemes.

### Machine Learning
Machine learning models employ deep MLP architectures to estimate the interfacial curvarture. Models are developed in three steps. First, a synthetic dataset composed of circular arcs of varying sizes is generated. Subsequently, about hundred MLP models are trained in TensorFlow for each hyperparameters configuration. There is a scatter in model performance due to inherent stochasticity involved in training process. Finally, the best performing models are selected, and the weights and biases of these Tensorflows models are converted to standard ascii format. The *multilayerPerceptron* class in OpenFOAM reads these parameters and constructs the corresponding MLP model. 

## Installation
Make sure you have installed OpenFOAM v2006 and run the following scripts in the project directory:
```
./Allwclean
./Allwmake
```

## Examples 
Several 2D benchmark cases are provided.

## Prerequisites
OpenFOAM v2006.

## Author
Asim Önder



