Library of new reconstruction schemes and surface tension models for two-phase flow simulations in OpenFOAM.  Compatible with OpenFOAM v2006.

The methods are designed for uniform Cartesian grids. In cases where the interface is restriced to a subset of a domain, e.g.,problems involving deep water waves, the uniform grid can also be restricted to a smaller zone containing the interface. This is currently achieved using the cellSet functionality. The rest of the domain can use unstructured grids. Such an approach allows best of both worlds: the accuracy and speed of Cartesian methods around the interface where they are needed, and the flexibility of the unstructured grids elsewhere. 

The library is still in development. At the moment, the available methods are:
Interface normal: Young's method, central-columns differences, mixed Young-central method
Interface curvature: height function method

An example for gravity-capillary waves is provided. More benchmark cases will follow.

A paper documenting the details of the library will be added when ready.

Author: Asim Onder



