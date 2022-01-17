
# machineLearningVOF
An extension library containing new neural-network based VOF models for two-phase flow simulations in OpenFOAM. This is an alpha pre-release containing only a portion of the package. The whole package along with example cases will be released once testing is complete.  

## Description
New interface reconstruction schemes and surface tension models are provided. The models employ multilayer perceptron architectures to estimate interface normal vector and curvarture. Multiplayer perceptrons are implemented using C++ Standard Library and does not require any additional package. 

The methods are designed for uniform Cartesian grids. In cases where the interface is restricted to a subset of the domain, e.g., deep water waves, the uniform grid can also be restricted to a smaller zone containing the interface. This is currently achieved using the cellSet functionality. The rest of the domain can use unstructured grids. Such an approach offers the best of both worlds: the accuracy and speed of Cartesian methods around the interface where they are needed, and the flexibility of the unstructured grids elsewhere. 

In addition to neural-network models, several popular conentional schemes are also implemented: (i) interface normal: Youngs' Method, Central-Columns Differences Method, Mixed Youngs-Central Method; (ii) interface curvature: the Height Function Method.

## Example 
Examples and benchmark cases will be provided soon.

## Prerequisites
OpenFOAM v2006.

## Author
Asim Önder



