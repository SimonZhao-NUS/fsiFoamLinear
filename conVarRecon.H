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
