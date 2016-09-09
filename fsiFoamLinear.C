/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    rhoCentralDyMFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "myHeader.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
//    #include "createDynamicFvMesh.H"
    #include "createMesh.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    #include "readIfLinear.H"
//    #include "readTimeScheme.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "readTimeControls.H"
    #include "readFarFieldValue.H"


    Info<< "\nStarting time loop\n" << endl;

    //////bd corr//////////
    thermo.correct();
    volScalarField rPsi(1.0/psi);
    volScalarField c(sqrt(thermo.gamma()*rPsi));
    e = rhoE/rho - 0.5*magSqr(U);
    p = rho*e*(thermo.gamma()-1);
    U.dimensionedInternalField() =
	rhoU.dimensionedInternalField()
       /rho.dimensionedInternalField();
    farFieldBoundaryCondition(rho, U, p, c, U_inf, rho_inf, p_inf);
    T = p/rho/(thermo.gamma()-1)/thermo.Cv();
    thermo.correct();
    rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
    rhoE.boundaryField() =
	rho.boundaryField()*
	(
	    e.boundaryField() + 0.5*magSqr(U.boundaryField())
	);
    //////bd corr//////

    while ( runTime.run() || (ifLinear == "Yes"))
    {
	#include "setDt.H"
        runTime++;
	
	//p boundary condition
        e = rhoE/rho - 0.5*magSqr(U);
        p = rho*e*(thermo.gamma()-1);
/*        p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();*/
//        p.correctBoundaryConditions();
	//rho boundary condition ***
//	rho.boundaryField() = psi.boundaryField()*p.boundaryField();
	//U boundary condition
	U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
//        U.correctBoundaryConditions();
	farFieldBoundaryCondition(rho, U, p, c, U_inf, rho_inf, p_inf);

	//e boundary condition ***
//        e.correctBoundaryConditions();
//        e.boundaryField() = p.boundaryField()/rho.boundaryField()/(thermo.gamma()().boundaryField()-1);
//	T.boundaryField() = p.boundaryField()/rho.boundaryField()/(thermo.gamma()().boundaryField()-1)/thermo.Cv()().boundaryField();
	T = p/rho/(thermo.gamma()-1)/thermo.Cv();
        thermo.correct();
        //
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
        rhoE.boundaryField() =
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );
/////////boundary condtition end ////////////////

        // --- upwind interpolation of primitive fields on faces
        surfaceScalarField gamma
        (
	    "gamma",
            fvc::interpolate(thermo.gamma(), "reconstruct(gamma)")
        );
        surfaceScalarField rho_pos
        (
	    "rho_pos",
            myInterpolatePos(rho)
        );

        surfaceScalarField rho_neg
        (
	    "rho_neg",
            myInterpolateNeg(rho)
        );
     
        surfaceVectorField rhoU_pos
        (
	    "rhoU_pos",
            myInterpolatePos(rhoU)
        );
        surfaceVectorField rhoU_neg
        (
	    "rhoU_neg",
            myInterpolateNeg(rhoU)
        );

        surfaceScalarField rhoE_pos
        (
	    "rhoE_pos",
            myInterpolatePos(rhoE)
        );
        surfaceScalarField rhoE_neg
        (
	    "rhoE_neg",
            myInterpolateNeg(rhoE)
        );

/*	myMat matQ;
	for(unsigned int ni = 0; ni < mesh.nCells(); ni++) {
	    matQ.push_back(std::vector<double>(1,rho.internalField()[ni]));
	    matQ.push_back(std::vector<double>(1,rhoU.internalField()[ni].x()));
	    matQ.push_back(std::vector<double>(1,rhoU.internalField()[ni].y()));
	    matQ.push_back(std::vector<double>(1,rhoU.internalField()[ni].z()));
	    matQ.push_back(std::vector<double>(1,rhoE.internalField()[ni]));
	}
	myMatWrite(matQ, "matQ");
	
	SpMat matQf_pos(mesh.owner().size()*5, 1);
	SpMat matQf_neg(mesh.owner().size()*5, 1);
	TripletList QfList_pos, QfList_neg;
	for (int fi=0; fi<mesh.owner().size(); fi++) {
	    QfList_pos.push_back(Triplet(fi*5+0, 0, rho_pos[fi]));
	    QfList_neg.push_back(Triplet(fi*5+0, 0, rho_neg[fi]));
	}
	matQf_pos.setFromTriplets(QfList_pos.begin(), QfList_pos.end());
	matQf_neg.setFromTriplets(QfList_neg.begin(), QfList_neg.end());

	myMatWriteDenSp(matQf_pos, "matQf_pos");	
	myMatWriteDenSp(matQf_neg, "matQf_neg");
	Info << "Qf" << endl;
	cin.get();*/
	//
        //
	volScalarField rPsi(1.0/psi);
	volScalarField c(sqrt(thermo.gamma()*rPsi));
        volScalarField H(rhoE/rho+rPsi);
	//
        surfaceVectorField U_pos(rhoU_pos/rho_pos); 
	surfaceVectorField U_neg(rhoU_neg/rho_neg);
	
        surfaceScalarField rPsi_pos((gamma-1)*(rhoE_pos/rho_pos-.5*magSqr(U_pos))); 
	surfaceScalarField rPsi_neg((gamma-1)*(rhoE_neg/rho_neg-.5*magSqr(U_neg)));

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos); 
	surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);
	
        surfaceScalarField E_pos
        (
	    rhoE_pos/rho_pos
        );

        surfaceScalarField E_neg
        (
	    rhoE_neg/rho_neg
        );

        surfaceScalarField H_pos
        (
	    E_pos + rPsi_pos
        );
        surfaceScalarField H_neg
        (
	    E_neg + rPsi_neg
        );

	surfaceScalarField dU1(rho_pos*(U_pos&mesh.Sf()));
	surfaceVectorField dU234(p_pos*mesh.Sf());
	surfaceScalarField dU5(rho_pos*H_pos*(U_pos&mesh.Sf()));
	
	surfaceScalarField rho_tilt(rho_pos);
	surfaceVectorField U_tilt(U_pos);
	surfaceScalarField H_tilt(H_pos);

	roeAverage
	(
	     dU1,
	     dU234,
	     dU5,
	     rho_tilt,
	     U_tilt,
	     H_tilt,
	     rho_pos,
	     rho_neg,
	     U_pos,
	     U_neg,
	     H_pos,
	     H_neg,
	     p_pos,
	     p_neg,
	     gamma,
	     U
	);

//	#include "setDt.H"
//        runTime++;
	Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
        Info<< "Time = " << runTime.timeName() << nl << endl;
	
	surfaceScalarField USf_pos(U_pos&mesh.Sf());
	surfaceScalarField USf_neg(U_neg&mesh.Sf());

	// -mesh.phi() == + xi_t	
	fvc::makeRelative(USf_pos, U);
	fvc::makeRelative(USf_neg, U);

	phi = 0.5*(rho_pos*(USf_pos) + rho_neg*(USf_neg)
		  - dU1);

	surfaceVectorField phiUp 
	(
	    "phiUp",
	    0.5*
	    (rho_pos*(USf_pos)*U_pos
	    + rho_neg*(USf_neg)*U_neg
            + p_pos*mesh.Sf()
	    + p_neg*mesh.Sf()
	    - dU234
	    )
	);

        surfaceScalarField phiEp
        (
	    "phiEp",
	    0.5*(rho_pos*H_pos*(USf_pos)
	    + rho_neg*H_neg*(USf_neg)
//	    + mesh.phi()*p_pos
//	    + mesh.phi()*p_neg
	    - dU5)
	);

	if (ifLinear == "Yes") {
	    #include "linearizationSp.H"
	}

	// --- Solve density          
/*        solve(fvm::ddt(rho) + fvc::div(phi));

    // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));
      
      // --- Solve energy
	solve(fvm::ddt(rhoE) + fvc::div(phiEp));*/


	scalar dT = runTime.deltaTValue();
	rho.internalField() = rho.internalField() - dT*fvc::div(phi)().internalField();
	rhoU.internalField() = rhoU.internalField() - dT*fvc::div(phiUp)().internalField();
	for(unsigned int ni = 0; ni < rho.internalField().size(); ni++) {
	    rhoU.internalField()[ni].z() = 0;
	}
	rhoE.internalField() = rhoE.internalField() - dT*fvc::div(phiEp)().internalField();


	//p boundary condition
        e = rhoE/rho - 0.5*magSqr(U);
        p = rho*e*(thermo.gamma()-1);
//        p.correctBoundaryConditions();
	//rho boundary condition ***
//	rho.boundaryField() = psi.boundaryField()*p.boundaryField();
	//U boundary condition
	U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
//        U.correctBoundaryConditions();
	farFieldBoundaryCondition(rho, U, p, c, U_inf, rho_inf, p_inf);

	//e boundary condition ***
	T = p/rho/(thermo.gamma()-1)/thermo.Cv();
        thermo.correct();
        //
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
        rhoE.boundaryField() =
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
//
//
/*	if (ifLinear == "Yes") {
	    #include "linearizationSp.H"
	}

	// --- Solve density          
        solve(fvm::ddt(rho) + fvc::div(phi));

    // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));
      
      // --- Solve energy
	solve(fvm::ddt(rhoE) + fvc::div(phiEp));

	//p boundary condition
        e = rhoE/rho - 0.5*magSqr(U);
        p = rho*e*(thermo.gamma()-1);
//        p.correctBoundaryConditions();
	//rho boundary condition ***
//	rho.boundaryField() = psi.boundaryField()*p.boundaryField();
	//U boundary condition
	U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
//        U.correctBoundaryConditions();
	farFieldBoundaryCondition(rho, U, p, c, U_inf, rho_inf, p_inf);

	//e boundary condition ***
	T = p/rho/(thermo.gamma()-1)/thermo.Cv();
        thermo.correct();
        //
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
        rhoE.boundaryField() =
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );

        runTime.write();
	*/


/*    forAll(mesh.boundary(), patchi) {
	if (mesh.boundary()[patchi].name() != "inlet" &&
	    mesh.boundary()[patchi].name() != "outlet"
	) {
	    for (int n = 0; n < mesh.boundary()[patchi].size(); n++) {
		int fn = mesh.boundary()[patchi].start() + n;
		forAll (mesh.faces()[fn], p) {
		    int pn = mesh.faces()[fn][p];
		    scalar Lx = mesh.points()[pn].x() - .25;
		    scalar Ly = mesh.points()[pn].y();
		    scalar dy = Lx*alpha;
		    scalar dx = -Ly*alpha;
		    newPoints[pn].x() = mesh.points()[pn].x()+dx;
		    newPoints[pn].y() = mesh.points()[pn].y()+dy;
		}
	    }
	}
    }*/
