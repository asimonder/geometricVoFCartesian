
# machineLearningVOF
An extension library containing new neural-network based VOF models for two-phase flow simulations in OpenFOAM. This is an alpha pre-release containing only a portion of the package. The whole package along with example cases will be released soon.  

## Description
New interface reconstruction schemes and surface tension models are provided. The models employ multilayer perceptron (MLP) architectures to estimate interface normal vector and curvarture. Models are developed in three steps. First, MLPs are trained in Tensorflow using synthetic datasets composed of either circles or waves of varying sizes. Subsequently, the weights and biases of Tensorflows models are converted to ascii format. Finally, a newly implemented MLP class in OpenFOAM reads these parameters and constructs the corresponding MLP model. The basic functionalities of MLPs are implemented into OpenFOAM using the C++ Standard Library. Thus, no additional package is required to use new machine learning models.  

The methods are designed for uniform Cartesian grids. In cases where the interface is restricted to a subset of the domain, e.g., deep water waves, the uniform grid can also be restricted to a smaller zone containing the interface. This is currently achieved using the cellSet functionality. The rest of the domain can use unstructured grids. Such an approach offers the best of both worlds: the accuracy and speed of Cartesian methods around the interface where they are needed, and the flexibility of the unstructured grids elsewhere. 

In addition to neural-network models, several popular conentional schemes are also implemented: (i) interface normal: Youngs' Method, Central-Columns Differences Method, Mixed Youngs-Central Method; (ii) interface curvature: the Height Function Method.

## Example 
Examples and benchmark cases will be provided soon.

## Prerequisites
OpenFOAM v2006.

## Author
Asim Önder



