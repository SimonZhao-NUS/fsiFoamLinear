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
    fsiFoamLinear

Description
    Roe based fvm solver. Dynamic linearization implemented.
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
//    #include "readFarFieldValue.H"


    Info<< "\nStarting time loop\n" << endl;


    while ( runTime.run() || (ifLinear == "Yes"))
    {
	//reconstruct conservative flow variables on cell interfaces
	#include "conVarRecon.H"

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

	#include "setDt.H"
        runTime++;
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
        solve(fvm::ddt(rho) + fvc::div(phi));

	// --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));
        U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
        U.correctBoundaryConditions();
        rhoU.boundaryField() == rho.boundaryField()*U.boundaryField();
      
        // --- Solve energy
	solve(fvm::ddt(rhoE) + fvc::div(phiEp));
        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryField() ==
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );
        p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();
        p.correctBoundaryConditions();
        rho.boundaryField() == psi.boundaryField()*p.boundaryField();

        turbulence->correct();

	// update values
	// correct boundary conditions
	// correct conservative variables
	// correct intermediate variables


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
/*
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


	//p boundary condition
        e = rhoE/rho - 0.5*magSqr(U);
        p = rho*e*(thermo.gamma()-1);
//        p.dimensionedInternalField() =
//            rho.dimensionedInternalField()
//           /psi.dimensionedInternalField();
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
*/
