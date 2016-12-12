        // --- upwind interpolation of primitive fields on faces
        surfaceScalarField gammaSurface
        (
	    "gammaSurface",
            fvc::interpolate(gamma, "reconstruct(gamma)")
        );
	
	scalar K_lmt = 20;
        const volVectorField rhoGrad = fvc::grad(rho, "Gauss linear");
	szhaoFoam::vkSlopeLimiter<scalar> rhoLimiter(rho, rhoGrad, K_lmt);

        surfaceScalarField rho_pos
        (
	    "rho_pos",
            myInterpolatePos(rho, rhoGrad, rhoLimiter)
        );

        surfaceScalarField rho_neg
        (
	    "rho_neg",
            myInterpolateNeg(rho, rhoGrad, rhoLimiter)
        );
     
        const volTensorField rhoUGrad = fvc::grad(rhoU, "Gauss linear");
	szhaoFoam::vkSlopeLimiter<vector> rhoULimiter(rhoU, rhoUGrad, K_lmt);

        surfaceVectorField rhoU_pos
        (
	    "rhoU_pos",
            myInterpolatePos(rhoU, rhoUGrad, rhoULimiter)
        );
        surfaceVectorField rhoU_neg
        (
	    "rhoU_neg",
            myInterpolateNeg(rhoU, rhoUGrad, rhoULimiter)
        );

        const volVectorField rhoEGrad = fvc::grad(rhoE, "Gauss linear");
	szhaoFoam::vkSlopeLimiter<scalar> rhoELimiter(rhoE, rhoEGrad, K_lmt);
/*Info << rhoELimiter.UmaxInd() << endl;
cin.get();
Info << rhoELimiter.UmaxIndBd() << endl;
cin.get();*/
        surfaceScalarField rhoE_pos
        (
	    "rhoE_pos",
            myInterpolatePos(rhoE, rhoEGrad, rhoELimiter)
        );
        surfaceScalarField rhoE_neg
        (
	    "rhoE_neg",
            myInterpolateNeg(rhoE, rhoEGrad, rhoELimiter)
        );
	// end surface interpolation

