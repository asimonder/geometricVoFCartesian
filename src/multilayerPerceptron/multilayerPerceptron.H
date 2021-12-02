/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2021 Asim Onder
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::multilayerPerceptron

Description
    Developed by Asim Onder (2021)


\*---------------------------------------------------------------------------*/

#ifndef multilayerPerceptron_H
#define multilayerPerceptron_H


#include "fvMesh.H"
/*#include "ijkZone.H"
#include "volFields.H"
#include "surfaceMesh.H"
#include "volMesh.H"
#include "IOobject.H"
#include "MeshObject.H"
#include "IFstream.H"
#include "cellSet.H"*/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include<complex>
#include<cstdlib>
#include <vector>
#include <numeric>
#include <math.h>



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class multilayerPerceptron
{
  std::string fileName_;
  int NLayers_;
  std::vector<std::vector<std::vector<double>>> weights_;
  std::vector<std::vector<double>> biases_;
  std::vector<std::string> activations_;

  //- No copy construct
  multilayerPerceptron(const multilayerPerceptron&) = delete;
  
  //- No bitwise assignment
  //void operator=(const multilayerPerceptron&) = delete;


public:
  multilayerPerceptron()=default;
  
  multilayerPerceptron(std::string filename);
  //- Destructor
  virtual ~multilayerPerceptron() = default;

  std::vector<std::string>  activations();
  std::vector<std::vector<std::vector<double>>>  weights();
  std::vector<std::vector<double>> biases();
  std::vector<double> predict(const std::vector<double>& data);
  int stencilSize()
  {
    return weights_[0][0].size();
  }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam


#endif

// ************************************************************************* //
