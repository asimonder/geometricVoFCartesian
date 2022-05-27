
# geometricVoFCartesian
geometricVoFCartesian is an extension library for two-phase modelling in OpenFOAM. It features conventional and machine-learning models to estimate the normal vector and curvature of the fluid-fluid interface in geometric Volume-of-Fluid (VOF) framework. 

## Features
The methods are designed for uniform (isotropic) Cartesian grids. In cases where the interface is restricted to a subset of the domain, e.g., deep water waves, the uniform grid can also be restricted to a smaller zone containing the interface. This is currently achieved using the cellSet functionality. The rest of the domain can use unstructured grids. 


The models employ multilayer perceptron (MLP) architectures to estimate interface normal vector and curvarture. 
Models are developed in three steps. First, MLPs are trained in Tensorflow using synthetic datasets composed of either circles or waves of varying sizes. Subsequently, the weights and biases of Tensorflows models are converted to ascii format. Finally, a newly implemented MLP class in OpenFOAM reads these parameters and constructs the corresponding MLP model. The basic functionalities of MLPs are implemented into OpenFOAM using the C++ Standard Library. Thus, no additional package is required to use new machine learning models.  

The methods are designed for uniform Cartesian grids. In cases where the interface is restricted to a subset of the domain, e.g., deep water waves, the uniform grid can also be restricted to a smaller zone containing the interface. This is currently achieved using the cellSet functionality. The rest of the domain can use unstructured grids. 

In addition to neural-network models, several popular conentional schemes are also implemented: (i) interface normal: Youngs' Method, Central-Columns Differences Method, Mixed Youngs-Central Method; (ii) interface curvature: the Height Function Method.

## Example 
Benchmark cases are provided in soon.

## Prerequisites
OpenFOAM v2006.

## Author
Asim Önder



