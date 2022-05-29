
# geometricVoFCartesian
The *geometricVoFCartesian* is an extension library for simulating two-phase flows in OpenFOAM. It features conventional and machine-learning models to estimate the normal vector and curvature of the fluid-fluid interfaces in geometric Volume-of-Fluid (VOF) framework. 

### General Description
The methods are designed for uniform (isotropic) Cartesian grids. In cases where the interface motion occurs in a smaller subset of the domain, the uniform grid can also be restricted to a smaller zone containing the interface. This is currently achieved using the cellSet functionality. The rest of the domain can use generic unstructured grids. 

Available methods are as follows:
- **Interface normal vector**: Youngs' Method, Central-Columns Differences Method, Mixed Youngs-Central (MYC) Method
- **Interface curvature**: Height Function Method, Multilayer-Perceptron (MLP) models

MLP models currently support 2D geometries only. 



Machine learning models employ multilayer perceptron (MLP) architectures to estimate interface normal vector and curvarture. 
Models are developed in three steps. First, MLPs are trained in Tensorflow using synthetic datasets composed of either circles or waves of varying sizes. Subsequently, the weights and biases of Tensorflows models are converted to ascii format. Finally, a newly implemented MLP class in OpenFOAM reads these parameters and constructs the corresponding MLP model. The basic functionalities of MLPs are implemented into OpenFOAM using the C++ Standard Library. Thus, no additional package is required to use new machine learning models.  

In addition to neural-network models, several popular conentional schemes are also implemented: (i) interface normal: Youngs' Method, Central-Columns Differences Method, Mixed Youngs-Central Method; (ii) interface curvature: the Height Function Method.

## Example 
Benchmark cases are provided in soon.

## Prerequisites
OpenFOAM v2006.

## Author
Asim Önder



