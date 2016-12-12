void farFieldBoundaryCondition
(
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p
)
{
    const volScalarField& gamma = p.db().lookupObject<volScalarField>("mygamma");   
    const volScalarField c(sqrt(gamma*p/rho));
    const fvMesh& mesh = p.mesh();
    #include "readFarFieldValue.H"

    forAll(mesh.boundary(), patchi) {
	//Us: + for outflow, - for inflow
//	const scalarField Us = U.boundaryField()[patchi]&mesh.Sf().boundaryField()[patchi];
	scalarField Us = U.boundaryField()[patchi]&mesh.Sf().boundaryField()[patchi];
	if (mesh.moving()) {
	    Us -= mesh.phi().boundaryField()[patchi];
	}
	const scalarField Cs = c.boundaryField()[patchi]*mesh.magSf().boundaryField()[patchi];
	//flowRegime + for supersonic, - for subsonic
	const scalarField flowRegime = mag(Us) - Cs;

	fvPatchField<scalar>& pPatch = p.boundaryField()[patchi];
	fvPatchField<scalar>& rhoPatch = rho.boundaryField()[patchi];
	fvPatchField<vector>& UPatch = U.boundaryField()[patchi];
	const fvsPatchField<vector>& SfPatch = mesh.Sf().boundaryField()[patchi];
	const fvsPatchField<scalar>& magSfPatch = mesh.magSf().boundaryField()[patchi];
	const scalarField pPif = p.boundaryField()[patchi].patchInternalField();
	const scalarField rhoPif = rho.boundaryField()[patchi].patchInternalField();
	const vectorField UPif = U.boundaryField()[patchi].patchInternalField();
	const scalarField cPif = c.boundaryField()[patchi].patchInternalField();
	const labelUList& pFaceCells =
	mesh.boundary()[patchi].faceCells();

	if ((mesh.boundary()[patchi].name() == "inlet") || (mesh.boundary()[patchi].name() == "outlet")) {
	    forAll(mesh.boundary()[patchi], facei) {
		if(Us[facei] > 0) { // outflow
		    if(flowRegime[facei] > 0) { //supersonic
			pPatch[facei] = pPif[facei];
			rhoPatch[facei] = rhoPif[facei];
			UPatch[facei] = UPif[facei];
			UPatch[facei].z() = 0;
		    } else { //subsonic
			vector nf = SfPatch[facei]/magSfPatch[facei];
			pPatch[facei] = p_inf;
			rhoPatch[facei] = rhoPif[facei] + (pPatch[facei]-pPif[facei])/cPif[facei]/cPif[facei];
			UPatch[facei] = UPif[facei]-nf*(pPatch[facei]-pPif[facei])/rhoPif[facei]/cPif[facei];
			UPatch[facei].z() = 0;
		    }
		} else { // inflow
		    if(flowRegime[facei] > 0) { //supersonic
			pPatch[facei] = p_inf;
			rhoPatch[facei] = rho_inf;
			UPatch[facei] = U_inf;
			UPatch[facei].z() = 0;
		    } else { //subsonic
			vector nf = SfPatch[facei]/magSfPatch[facei];
			pPatch[facei] = .5*(p_inf+pPif[facei]-rhoPif[facei]*cPif[facei]*((U_inf-UPif[facei])&nf));
			rhoPatch[facei] = rho_inf + (pPatch[facei]-p_inf)/cPif[facei]/cPif[facei];
			UPatch[facei] = U_inf-nf*(p_inf-pPatch[facei])/rhoPif[facei]/cPif[facei];
			UPatch[facei].z() = 0;
		    }
		}
	    }
	}
	//make this boundary confirmed to OpenFOAM rules, note z direction
       	if (mesh.boundary()[patchi].name() == "wing") {
	    forAll(mesh.boundary()[patchi], facei) {
		vector nf = SfPatch[facei]/magSfPatch[facei];
		pPatch[facei] = pPif[facei];
		rhoPatch[facei] = rhoPif[facei];
//		UPatch[facei] = UPif[facei] - (UPif[facei]&nf)*nf;
		if (mesh.moving()) {
//		    UPatch[facei] = UPif[facei] - (UPif[facei]&nf)*nf;
		    UPatch[facei] = UPif[facei] - ((UPif[facei]&nf) - mesh.phi().boundaryField()[patchi][facei]/magSfPatch[facei])*nf ;
		} else {
		    UPatch[facei] = UPif[facei] - (UPif[facei]&nf)*nf;
		}
		UPatch[facei].z() = 0;

/*		Info<< nf <<endl;
		Info<< UPif[facei] <<endl;
		Info<< UPatch[facei] <<endl;
		cin.get();*/
	    }
	}
    }
/*	Info<< U <<endl;
	cin.get();
	Info<< p <<endl;
	cin.get();
	Info<< rho <<endl;
	cin.get();*/
}

/*
		Info<< nf <<endl;
		Info<< UPif[facei] <<endl;
		Info<< U.boundaryField()[patchi].patchInternalField()()[facei] <<endl;
		Info<< U.internalField()[pFaceCells[facei]] <<endl;
		Info<< UPatch[facei] <<endl;
		cin.get();
 * */
