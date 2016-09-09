void cenDiff
(
    surfaceScalarField& cen1,
    surfaceVectorField& cen234,
    surfaceScalarField& cen5,
    const surfaceScalarField& rho_pos,
    const surfaceScalarField& rho_neg,
    const surfaceVectorField& U_pos,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& H_pos,
    const surfaceScalarField& H_neg,
    const surfaceScalarField& p_pos,
    const surfaceScalarField& p_neg,
    const surfaceVectorField& rhoU_pos,
    const surfaceVectorField& rhoU_neg,
    const surfaceScalarField& rhoE_pos,
    const surfaceScalarField& rhoE_neg,
    const surfaceScalarField& gamma
)
{
    const fvMesh& mesh = rho_pos.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    const labelUList& owner = mesh.owner();

    forAll(owner, facei) {
	const vector sf = Sf.internalField()[facei];
	const vector U = U_pos.internalField()[facei];
	const scalar H = H_pos.internalField()[facei];
	const scalar gamma_scalar = gamma.internalField()[facei];

	myMat dQ
	(
	    5,
	    std::vector<double>(1,0)
	);
	dQ[0][0] = rho_pos.internalField()[facei];
	dQ[1][0] = rhoU_pos.internalField()[facei].x();
	dQ[2][0] = rhoU_pos.internalField()[facei].y(); 
	dQ[3][0] = rhoU_pos.internalField()[facei].z(); 
	dQ[4][0] = rhoE_pos.internalField()[facei]; 

	myMat eigREigVal
	(
	    innProd 
	    (
		calcEigVecR
		(
		    sf,
		    U,
		    H,
		    gamma_scalar
		),
		calcEigVal
		(
		    sf,
		    U,
		    H,
		    gamma_scalar
		)
	    )
	);
	myMat eigLdQ
	(
	    innProd 
	    (
		calcEigVecL
		(
		    sf,
		    U,
		    H,
		    gamma_scalar
		),
		dQ
	    )
	);
	
	myMat matdU
	(
	    innProd 
	    (
		eigREigVal,
		eigLdQ
	    )
	);

	cen1.internalField()[facei] = matdU[0][0];
	cen234.internalField()[facei].x() = matdU[1][0];
	cen234.internalField()[facei].y() = matdU[2][0]; 
	cen234.internalField()[facei].z() = matdU[3][0]; 
	cen5.internalField()[facei] = matdU[4][0];
    }
    forAll(owner, facei) {
	const vector sf = Sf.internalField()[facei];
	const vector U = U_neg.internalField()[facei];
	const scalar H = H_neg.internalField()[facei];
	const scalar gamma_scalar = gamma.internalField()[facei];

	myMat dQ
	(
	    5,
	    std::vector<double>(1,0)
	);
	dQ[0][0] = rho_neg.internalField()[facei];
	dQ[1][0] = rhoU_neg.internalField()[facei].x();
	dQ[2][0] = rhoU_neg.internalField()[facei].y(); 
	dQ[3][0] = rhoU_neg.internalField()[facei].z(); 
	dQ[4][0] = rhoE_neg.internalField()[facei]; 

	myMat eigREigVal
	(
	    innProd 
	    (
		calcEigVecR
		(
		    sf,
		    U,
		    H,
		    gamma_scalar
		),
		calcEigVal
		(
		    sf,
		    U,
		    H,
		    gamma_scalar
		)
	    )
	);
	myMat eigLdQ
	(
	    innProd 
	    (
		calcEigVecL
		(
		    sf,
		    U,
		    H,
		    gamma_scalar
		),
		dQ
	    )
	);
	
	myMat matdU
	(
	    innProd 
	    (
		eigREigVal,
		eigLdQ
	    )
	);

	cen1.internalField()[facei] += matdU[0][0];
	cen234.internalField()[facei].x() += matdU[1][0];
	cen234.internalField()[facei].y() += matdU[2][0]; 
	cen234.internalField()[facei].z() += matdU[3][0]; 
	cen5.internalField()[facei] += matdU[4][0];
    }
}
/*
void roeAverage
(
    surfaceScalarField& dU1,
    surfaceVectorField& dU234,
    surfaceScalarField& dU5,
    surfaceScalarField& rho_tilt,
    surfaceVectorField& U_tilt,
    surfaceScalarField& H_tilt,
    const surfaceScalarField& rho_pos,
    const surfaceScalarField& rho_neg,
    const surfaceVectorField& U_pos,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& H_pos,
    const surfaceScalarField& H_neg,
    const surfaceScalarField& p_pos,
    const surfaceScalarField& p_neg,
    const surfaceVectorField& rhoU_pos,
    const surfaceVectorField& rhoU_neg,
    const surfaceScalarField& rhoE_pos,
    const surfaceScalarField& rhoE_neg,
    const surfaceScalarField& gamma
)
{
	const fvMesh& mesh = rho_pos.mesh();
    
	rho_tilt =
        (
	    sqrt(rho_pos*rho_neg)
	);

        U_tilt =
        (
	    (U_neg+U_pos*sqrt(rho_pos/rho_neg))/(1+sqrt(rho_pos/rho_neg))
	);

	H_tilt =
        (
            (H_neg+H_pos*sqrt(rho_pos/rho_neg))/(1+sqrt(rho_pos/rho_neg))
	);
        
	surfaceScalarField U_bar_tilt
        (
	    (U_tilt&mesh.Sf())/mesh.magSf()
	);

	surfaceScalarField c_tilt
	(
 	    sqrt((gamma-1)*(H_tilt-magSqr(U_tilt)/2))
	);

        surfaceScalarField delta_rho
        (
	    (rho_neg-rho_pos)
	);

        surfaceScalarField delta_p
        (
	    (p_neg-p_pos)
	);

        surfaceVectorField delta_U
        (
	    (U_neg-U_pos)
	);

        surfaceScalarField delta_U_bar
        (
	    ((delta_U&mesh.Sf())/mesh.magSf())
	);

	surfaceScalarField alpha1
	(
	    mesh.magSf()*mag(U_bar_tilt)*(delta_rho-delta_p/magSqr(c_tilt))
	);

	surfaceScalarField alpha2
	(
	    mesh.magSf()*mag(U_bar_tilt+c_tilt)*(delta_p+rho_tilt*c_tilt*delta_U_bar)/(2*magSqr(c_tilt))
	);

	surfaceScalarField alpha3
	(
	    mesh.magSf()*mag(U_bar_tilt-c_tilt)*(delta_p-rho_tilt*c_tilt*delta_U_bar)/(2*magSqr(c_tilt))
	);

	surfaceScalarField alpha4
	(
	    alpha1+alpha2+alpha3
	);
	
	surfaceScalarField alpha5
	(
	    c_tilt*(alpha2-alpha3)
	);
	
	surfaceVectorField alpha678
	(
	    mesh.magSf()*mag(U_bar_tilt)*(rho_tilt*delta_U-mesh.Sf()/mesh.magSf()*rho_tilt*delta_U_bar)
	);


//        Info<< "Time = " << runTime.timeName() << nl << endl;

//	mesh.update();

    const surfaceVectorField& Sf = mesh.Sf();
    const labelUList& owner = mesh.owner();

    surfaceVectorField delta_rhoU
    (
	(rhoU_neg-rhoU_pos)
    );
    
    surfaceScalarField delta_rhoE
    (
	(rhoE_neg-rhoE_pos)
    );

    forAll(owner, facei) {
	const vector sf = Sf.internalField()[facei];
	const vector U = U_tilt.internalField()[facei];
	const scalar H = H_tilt.internalField()[facei];
	const scalar gamma_scalar = gamma.internalField()[facei];

	myMat dQ
	(
	    5,
	    std::vector<double>(1,0)
	);
	dQ[0][0] = delta_rho.internalField()[facei];
	dQ[1][0] = delta_rhoU.internalField()[facei].x();
	dQ[2][0] = delta_rhoU.internalField()[facei].y(); 
	dQ[3][0] = delta_rhoU.internalField()[facei].z(); 
	dQ[4][0] = delta_rhoE.internalField()[facei]; 

	myMat eigREigVal_abs
	(
	    innProd 
	    (
		calcEigVecR
		(
		    sf,
		    U,
		    H,
		    gamma_scalar
		),
		calcEigVal_abs
		(
		    sf,
		    U,
		    H,
		    gamma_scalar
		)
	    )
	);
	myMat eigLdQ
	(
	    innProd 
	    (
		calcEigVecL
		(
		    sf,
		    U,
		    H,
		    gamma_scalar
		),
		dQ
	    )
	);
	
	myMat matdU
	(
	    innProd 
	    (
		eigREigVal_abs,
		eigLdQ
	    )
	);

	dU1.internalField()[facei] = matdU[0][0];
	dU234.internalField()[facei].x() = matdU[1][0];
	dU234.internalField()[facei].y() = matdU[2][0]; 
	dU234.internalField()[facei].z() = matdU[3][0]; 
	dU5.internalField()[facei] = matdU[4][0];
    }
}
*/

void roeAverage
(
    surfaceScalarField& dU1,
    surfaceVectorField& dU234,
    surfaceScalarField& dU5,
    surfaceScalarField& rho_tilt,
    surfaceVectorField& U_tilt,
    surfaceScalarField& H_tilt,
    const surfaceScalarField& rho_pos,
    const surfaceScalarField& rho_neg,
    const surfaceVectorField& U_pos,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& H_pos,
    const surfaceScalarField& H_neg,
    const surfaceScalarField& p_pos,
    const surfaceScalarField& p_neg,
    const surfaceScalarField& gamma,
    const volVectorField& U
)
{
	const fvMesh& mesh = rho_pos.mesh();
    
	rho_tilt =
        (
	    sqrt(rho_pos*rho_neg)
	);

        U_tilt =
        (
	    (U_neg+U_pos*sqrt(rho_pos/rho_neg))/(1+sqrt(rho_pos/rho_neg))
	);

	H_tilt =
        (
            (H_neg+H_pos*sqrt(rho_pos/rho_neg))/(1+sqrt(rho_pos/rho_neg))
	);
        
	surfaceScalarField U_tiltSf
        (
	    (U_tilt&mesh.Sf())
	);
	////for dynamic mesh
	fvc::makeRelative(U_tiltSf, U);
	
	surfaceScalarField U_bar_tilt
        (
	    (U_tiltSf)/mesh.magSf()
	);

	surfaceScalarField c_tilt
	(
 	    sqrt((gamma-1)*(H_tilt-magSqr(U_tilt)/2))
	);

        surfaceScalarField delta_rho
        (
	    (rho_neg-rho_pos)
	);

        surfaceScalarField delta_p
        (
	    (p_neg-p_pos)
	);

        surfaceVectorField delta_U
        (
	    (U_neg-U_pos)
	);

        surfaceScalarField delta_U_bar
        (
	    ((delta_U&mesh.Sf())/mesh.magSf())
	);

	surfaceScalarField alpha1
	(
	    mesh.magSf()*mag(U_bar_tilt)*(delta_rho-delta_p/magSqr(c_tilt))
	);

	surfaceScalarField alpha2
	(
	    mesh.magSf()*mag(U_bar_tilt+c_tilt)*(delta_p+rho_tilt*c_tilt*delta_U_bar)/(2*magSqr(c_tilt))
	);

	surfaceScalarField alpha3
	(
	    mesh.magSf()*mag(U_bar_tilt-c_tilt)*(delta_p-rho_tilt*c_tilt*delta_U_bar)/(2*magSqr(c_tilt))
	);

	surfaceScalarField alpha4
	(
	    alpha1+alpha2+alpha3
	);
	
	surfaceScalarField alpha5
	(
	    c_tilt*(alpha2-alpha3)
	);
	
	surfaceVectorField alpha678
	(
	    mesh.magSf()*mag(U_bar_tilt)*(rho_tilt*delta_U-mesh.Sf()/mesh.magSf()*rho_tilt*delta_U_bar)
	);


//        Info<< "Time = " << runTime.timeName() << nl << endl;

//	mesh.update();

	dU1 =  alpha4;
        dU234 =
	(
	    (U_tilt*alpha4
	    + mesh.Sf()/mesh.magSf()*alpha5
	    + alpha678)
	);
        
	dU5 =
	(
	      (H_tilt*alpha4
//	    + U_bar_tilt*alpha5
	    + (U_tilt&mesh.Sf())/mesh.magSf()*alpha5
	    + (U_tilt&alpha678)
	    - ((magSqr(c_tilt))/(gamma-1))*alpha1)
	);
}
