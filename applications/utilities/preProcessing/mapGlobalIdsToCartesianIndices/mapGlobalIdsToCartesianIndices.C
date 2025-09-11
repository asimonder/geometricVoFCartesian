/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Asim Önder
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

Application
    interIsoFoamCart


Description
    IJK indexing of the global ids of the cells in a Cartesian zone. This offline mapping  avoids MPI problems in massively-parallel large cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "ijkZone.H" // Your ijkZone header

int main(int argc, char *argv[])
{
    // Set up the case to run in serial, reading the master mesh
    #include "setRootCase.H"
    #include "createTime.H"
    // The 'true' flag tells createMesh to ignore processor directories
    // and read the master constant/polyMesh
    #include "createMesh.H"

    Info<< "Running in serial on a single node to generate global map..." << endl;

    // Construct your ijkZone object. This runs the original, slow logic
    // but it's fine because we are in serial with lots of memory.
    ijkZone ijkMesh(mesh);

    // Now that the object is constructed, get the globalIds list
    // and write it to a file in the 'constant' directory.
    Info<< "Writing globalIds map to disk..." << endl;

    IOobject header
    (
        "globalIds.dat",
        runTime.constant(), // The path is constant/
        mesh,
        IOobject::NO_READ,
        IOobject::WRITE_IF_PRESENT,
        false
    );

    // Write the list to disk in ASCII format for inspection
    Ostream& os = header.writeStream(IOstream::ASCII);
    os << myZone.globalIds(); // You need a getter: const labelList& globalIds() const;
    os.flush();

    Info<< "Finished writing constant/globalIds.dat" << endl;
    Info<< "You can now recompile your library and run the parallel solver." << endl;

    return 0;
}
