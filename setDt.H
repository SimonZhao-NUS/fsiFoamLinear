/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar minDtLoc = 0.0;
double cfl = maxCo;

if (mesh.nInternalFaces())
{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    scalarField contraVSumX(mesh.nCells(), 0);
    scalarField contraVSumY(mesh.nCells(), 0);
    scalarField contraVSumZ(mesh.nCells(), 0);
    forAll(owner, fi) {
	double contraVX = .5*std::fabs(mesh.Sf()[fi].x());
	double contraVY = .5*std::fabs(mesh.Sf()[fi].y());
	double contraVZ = .5*std::fabs(mesh.Sf()[fi].z());
	contraVSumX[owner[fi]] += contraVX;
	contraVSumY[owner[fi]] += contraVY;
	contraVSumZ[owner[fi]] += contraVZ;
	contraVSumX[neighbour[fi]] += contraVX;
	contraVSumY[neighbour[fi]] += contraVY;
	contraVSumZ[neighbour[fi]] += contraVZ;
    }
//    volScalarField c(sqrt(thermo.Cp()/thermo.Cv()*rPsi));
    volScalarField c(sqrt(gamma*rPsi));
    scalarField dTLoc
    (
	cfl*mesh.V()/
	(
	    (mag(U.internalField().component(0))+c.internalField())*contraVSumX 
	    + (mag(U.internalField().component(1))+c.internalField())*contraVSumY 
	    + (mag(U.internalField().component(2))+c.internalField())*contraVSumZ 
	)
    );
    minDtLoc = min(dTLoc);
//    Info<< minDtLoc << endl;
//    cin.get();
    runTime.setDeltaT(minDtLoc);
}

Info<< "Courant Number:" << cfl
    << ", min locate dt: " << minDtLoc << endl;

// ************************************************************************* //
