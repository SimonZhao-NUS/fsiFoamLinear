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


    Info<< "\nStarting time loop\n" << endl;

/*    std::vector<std::vector<double> > map
    (
	mesh.nCells(),
	std::vector<double>
	(
	    7,
	    0
	)
    );

    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    forAll(P, fi) {
	map[P[fi]][map[P[fi]][0]+1] = fi;
	map[N[fi]][map[N[fi]][0]+1] = fi;
	
	++map[P[fi]][0];
	++map[N[fi]][0];	
    }

    int kkkk= 7449;
    Info << P[kkkk] << endl;
    Info << map[P[kkkk]][0] << endl;
    Info << map[P[kkkk]][1] << endl;
    Info << map[P[kkkk]][2] << endl;
    Info << map[P[kkkk]][3] << endl;
    Info << map[P[kkkk]][4] << endl;
    cin.get();
*/
//    #include "mode_a.H"
    while ( runTime.run() || (ifLinear == "Yes"))
    {
//	#include "setDt.H"
	runTime.setDeltaT(0.0009);

        runTime++;

	Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
        Info<< "Time = " << runTime.timeName() << nl << endl;
	
//	mesh.movePoints(newPoints);

	#include "varRecon.H"
/*    Info << rho[P[kkkk]] << endl;
    Info << map[P[kkkk]][0] << endl;
    Info << rho_pos[map[P[kkkk]][1]] << endl;
    Info << rho_pos[map[P[kkkk]][2]] << endl;
    Info << rho_pos[map[P[kkkk]][3]] << endl;
    Info << map[P[kkkk]][4] << endl;
    cin.get();
*/	
        surfaceVectorField U_pos(rhoU_pos/rho_pos); 
	surfaceVectorField U_neg(rhoU_neg/rho_neg);
	
        surfaceScalarField rPsi_pos((gammaSurface-1)*(rhoE_pos/rho_pos-.5*magSqr(U_pos))); 
	surfaceScalarField rPsi_neg((gammaSurface-1)*(rhoE_neg/rho_neg-.5*magSqr(U_neg)));

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos); 
	surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);
	
        surfaceScalarField E_pos(rhoE_pos/rho_pos);
        surfaceScalarField E_neg(rhoE_neg/rho_neg);

        surfaceScalarField H_pos(E_pos + rPsi_pos);
        surfaceScalarField H_neg(E_neg + rPsi_neg);

	surfaceScalarField dU1
	(
	    IOobject
	    (
		"dU1",
		runTime.timeName(),
		mesh
	    ),
	    mesh,
	    dimDensity*dimVelocity*dimArea
	);

	surfaceVectorField dU234
	(
	    IOobject
	    (
		"dU234",
		runTime.timeName(),
		mesh
	    ),
	    mesh,
	    dimPressure*dimArea
	);

	surfaceScalarField dU5
	(
	    IOobject
	    (
		"dU5",
		runTime.timeName(),
		mesh
	    ),
	    mesh,
	    dimEnergy/dimTime
	);

	surfaceScalarField rho_tilt
	(
	    IOobject
	    (
		"rho_tilt",
		runTime.timeName(),
		mesh
	    ),
	    mesh,
	    dimDensity
	);

	surfaceVectorField U_tilt
	(
	    IOobject
	    (
		"U_tilt",
		runTime.timeName(),
		mesh
	    ),
	    mesh,
	    dimVelocity
	);

	surfaceScalarField H_tilt
	(
	    IOobject
	    (
		"H_tilt",
		runTime.timeName(),
		mesh
	    ),
	    mesh,
	    dimEnergy/dimMass
	);

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
	     gammaSurface,
	     U
	);

	// begin to solve the governing equations
	
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

	if (mesh.moving()) {
	    phiEp += mesh.phi()*p_pos + mesh.phi()*p_neg;
	}

	if (ifLinear == "Yes") {
	    #include "linearizationSp.H"
	}

	word ifDiff("No");
	mesh.schemesDict().readIfPresent("ifDiff", ifDiff);
	if (ifDiff == "Yes")
	{
	    scalarField divPhi(fvc::div(phi)().internalField());
	    vectorField divPhiUp(fvc::div(phiUp)().internalField());
	    scalarField divPhiEp(fvc::div(phiEp)().internalField());
	    SpMat matBeq_diff(rho.internalField().size()*5,1);
	    SpMat matdQ(rho.internalField().size()*5,1);
	    SpMat matcellV(rho.internalField().size()*5,1);
	    SpMat matW(rho.internalField().size()*5,1);
	    TripletList BList, dQList, cellVList, WList;
	    BList.reserve(rho.internalField().size()*5);
	    dQList.reserve(rho.internalField().size()*5);
	    cellVList.reserve(rho.internalField().size()*5);

	    for(unsigned int ni = 0; ni < rho.internalField().size(); ni++) {
		BList.push_back(Triplet(ni*5, 0, divPhi[ni]));
		BList.push_back(Triplet(ni*5+1, 0, divPhiUp[ni].x()));
		BList.push_back(Triplet(ni*5+2, 0, divPhiUp[ni].y()));
		BList.push_back(Triplet(ni*5+3, 0, divPhiUp[ni].z()));
		BList.push_back(Triplet(ni*5+4, 0, divPhiEp[ni]));
		dQList.push_back(Triplet(ni*5, 0, rho[ni]));
		dQList.push_back(Triplet(ni*5+1, 0, rhoU[ni].x()));
		dQList.push_back(Triplet(ni*5+2, 0, rhoU[ni].y()));
		dQList.push_back(Triplet(ni*5+3, 0, rhoU[ni].z()));
		dQList.push_back(Triplet(ni*5+4, 0, rhoE[ni]));
		WList.push_back(Triplet(ni*5, 0, rho[ni]));
		WList.push_back(Triplet(ni*5+1, 0, U[ni].x()));
		WList.push_back(Triplet(ni*5+2, 0, U[ni].y()));
		WList.push_back(Triplet(ni*5+3, 0, U[ni].z()));
		WList.push_back(Triplet(ni*5+4, 0, p[ni]));
		cellVList.push_back(Triplet(ni*5, 0, mesh.V()[ni]));
		cellVList.push_back(Triplet(ni*5+1, 0, mesh.V()[ni]));
		cellVList.push_back(Triplet(ni*5+2, 0, mesh.V()[ni]));
		cellVList.push_back(Triplet(ni*5+3, 0, mesh.V()[ni]));
		cellVList.push_back(Triplet(ni*5+4, 0, mesh.V()[ni]));
	    }
	    matdQ.setFromTriplets(dQList.begin(), dQList.end());
	    matBeq_diff.setFromTriplets(BList.begin(), BList.end());
	    matcellV.setFromTriplets(cellVList.begin(), cellVList.end());
	    matW.setFromTriplets(WList.begin(), WList.end());
	    myMatWriteDenSp(matBeq_diff, "matBeq_diff");
	    myMatWriteDenSp(matdQ, "matdQ0");
	    myMatWriteDenSp(matW, "matW");
	
	    SpMat matQf_pos(mesh.owner().size()*5, 1);
	    SpMat matQf_neg(mesh.owner().size()*5, 1);
	    SpMat matf(mesh.owner().size()*5, 1);
	    TripletList QfList_pos, QfList_neg, fList;
	    for (int fi=0; fi<mesh.owner().size(); fi++) {
		QfList_pos.push_back(Triplet(fi*5+0, 0, rho_pos[fi]));
		QfList_pos.push_back(Triplet(fi*5+1, 0, rhoU_pos[fi].x()));
		QfList_pos.push_back(Triplet(fi*5+2, 0, rhoU_pos[fi].y()));
		QfList_pos.push_back(Triplet(fi*5+3, 0, rhoU_pos[fi].z()));
		QfList_pos.push_back(Triplet(fi*5+4, 0, rhoE_pos[fi]));
		QfList_neg.push_back(Triplet(fi*5+0, 0, rho_neg[fi]));
		QfList_neg.push_back(Triplet(fi*5+1, 0, rhoU_neg[fi].x()));
		QfList_neg.push_back(Triplet(fi*5+2, 0, rhoU_neg[fi].y()));
		QfList_neg.push_back(Triplet(fi*5+3, 0, rhoU_neg[fi].z()));
		QfList_neg.push_back(Triplet(fi*5+4, 0, rhoE_neg[fi]));
		fList.push_back(Triplet(fi*5, 0, phi[fi]));
		fList.push_back(Triplet(fi*5+1, 0, phiUp[fi].x()));
		fList.push_back(Triplet(fi*5+2, 0, phiUp[fi].y()));
		fList.push_back(Triplet(fi*5+3, 0, phiUp[fi].z()));
		fList.push_back(Triplet(fi*5+4, 0, phiEp[fi]));
	    }
	    matQf_pos.setFromTriplets(QfList_pos.begin(), QfList_pos.end());
	    matQf_neg.setFromTriplets(QfList_neg.begin(), QfList_neg.end());
	    matf.setFromTriplets(fList.begin(), fList.end());

	    myMatWriteDenSp(matQf_pos, "matQf_pos");	
	    myMatWriteDenSp(matQf_neg, "matQf_neg");
	    myMatWriteDenSp(matf, "matf");

	    unsigned int nBd = 0;
	    forAll(mesh.boundary(), patchi) {
		forAll(mesh.boundary()[patchi], facei) {
		    ++nBd;
		}
	    }
	    SpMat matfb(nBd*5, 1);
	    SpMat matQfb(nBd*5, 1);
	    SpMat matWfb(nBd*5, 1);
	    TripletList QfbList, WfbList, fbList;
	    unsigned int nBdIndex = 0;
	    forAll(mesh.boundary(), patchi) {
		forAll(mesh.boundary()[patchi], facei) {
		    QfbList.push_back(Triplet(nBdIndex*5+0, 0, rho.boundaryField()[patchi][facei]));
		    QfbList.push_back(Triplet(nBdIndex*5+1, 0, rhoU.boundaryField()[patchi][facei].x()));
		    QfbList.push_back(Triplet(nBdIndex*5+2, 0, rhoU.boundaryField()[patchi][facei].y()));
		    QfbList.push_back(Triplet(nBdIndex*5+3, 0, rhoU.boundaryField()[patchi][facei].z()));
		    QfbList.push_back(Triplet(nBdIndex*5+4, 0, rhoE.boundaryField()[patchi][facei]));
		    WfbList.push_back(Triplet(nBdIndex*5+0, 0, rho.boundaryField()[patchi][facei]));
		    WfbList.push_back(Triplet(nBdIndex*5+1, 0, U.boundaryField()[patchi][facei].x()));
		    WfbList.push_back(Triplet(nBdIndex*5+2, 0, U.boundaryField()[patchi][facei].y()));
		    WfbList.push_back(Triplet(nBdIndex*5+3, 0, U.boundaryField()[patchi][facei].z()));
		    WfbList.push_back(Triplet(nBdIndex*5+4, 0, p.boundaryField()[patchi][facei]));
		    fbList.push_back(Triplet(nBdIndex*5+0, 0, phi.boundaryField()[patchi][facei]));
		    fbList.push_back(Triplet(nBdIndex*5+1, 0, phiUp.boundaryField()[patchi][facei].x()));
		    fbList.push_back(Triplet(nBdIndex*5+2, 0, phiUp.boundaryField()[patchi][facei].y()));
		    fbList.push_back(Triplet(nBdIndex*5+3, 0, phiUp.boundaryField()[patchi][facei].z()));
		    fbList.push_back(Triplet(nBdIndex*5+4, 0, phiEp.boundaryField()[patchi][facei]));
		    ++nBdIndex;
		}
	    }
	    matQfb.setFromTriplets(QfbList.begin(), QfbList.end());
	    matWfb.setFromTriplets(WfbList.begin(), WfbList.end());
	    matfb.setFromTriplets(fbList.begin(), fbList.end());
	    myMatWriteDenSp(matQfb, "matQfb");
	    myMatWriteDenSp(matWfb, "matWfb");
	    myMatWriteDenSp(matfb, "matfb");
	}

	// --- Solve density          
        solve(fvm::ddt(rho) + fvc::div(phi));

	// --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));
      
        // --- Solve energy
	solve(fvm::ddt(rhoE) + fvc::div(phiEp));

	if (ifDiff == "Yes")
	{
	    SpMat matdQ(rho.internalField().size()*5,1);
	    TripletList dQList;
	    dQList.reserve(rho.internalField().size()*5);

	    for(unsigned int ni = 0; ni < rho.internalField().size(); ni++) {
		dQList.push_back(Triplet(ni*5, 0, rho[ni]));
		dQList.push_back(Triplet(ni*5+1, 0, rhoU[ni].x()));
		dQList.push_back(Triplet(ni*5+2, 0, rhoU[ni].y()));
		dQList.push_back(Triplet(ni*5+3, 0, rhoU[ni].z()));
		dQList.push_back(Triplet(ni*5+4, 0, rhoE[ni]));
	    }
	    matdQ.setFromTriplets(dQList.begin(), dQList.end());
	    myMatWriteDenSp(matdQ, "matdQ");

	    Info << "matBeq_diff written" << endl;
	    cin.get();
	}
	//correct boundary condtions
	//p, T, U, rho, psi, rPsi, e, h, rhoU, rhoE, 10 vars to correct
	U = rhoU/rho;
/*Info << Cv << endl;
Info <<"Cv" << endl;
cin.get();
Info << Cp << endl;
Info <<"Cp" << endl;
cin.get();
Info << e << endl;
Info <<"e" << endl;
cin.get();
Info << p << endl;
Info <<"p" << endl;
cin.get();
Info << gamma << endl;
Info <<"gamma" << endl;
cin.get();*/
	p = (rhoE/rho - .5*magSqr(U))*(gamma-1)*rho;
	farFieldBoundaryCondition(rho, U, p);
	T = p/(Rspec*rho);
	psi = rho/p;
	rPsi = p/rho;
	e = Cv*T;
	h = e+rPsi;
	rhoU = rho*U;
	rhoE = rho*(e + 0.5*magSqr(U));
/*Info << rho << endl;
Info <<"rho" << endl;
cin.get();
Info << rhoU << endl;
Info <<"rhoU" << endl;
cin.get();
Info << rhoE << endl;
Info <<"rhoE" << endl;
cin.get();*/

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
	scalar dT = runTime.deltaTValue();
	rho.internalField() = rho.internalField() - dT*fvc::div(phi)().internalField();
	rhoU.internalField() = rhoU.internalField() - dT*fvc::div(phiUp)().internalField();
	for(unsigned int ni = 0; ni < rho.internalField().size(); ni++) {
	    rhoU.internalField()[ni].z() = 0;
	}
	rhoE.internalField() = rhoE.internalField() - dT*fvc::div(phiEp)().internalField();

	//p boundary condition
        e = rhoE/rho - 0.5*magSqr(U);
        p = rho*e*(gamma()-1);
//        p.correctBoundaryConditions();
	//rho boundary condition ***
//	rho.boundaryField() = psi.boundaryField()*p.boundaryField();
	//U boundary condition
	U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
//        U.correctBoundaryConditions();
	farFieldBoundaryCondition(rho, U, p);

	//e boundary condition ***
	T = p/rho/(gamma()-1)/thermo.Cv();
        thermo.correct();
        //
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
        rhoE.boundaryField() =
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );
*/
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
