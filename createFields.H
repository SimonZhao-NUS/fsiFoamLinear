double U_freestream_dim = 1;
double p_freestream_dim = 1;
double rho_freestream_dim = 1;
double T_freestream_dim = 1;
double cSqr_freestream_dim = 1;

#include "readIfNondim.H"

if (ifNondim == "Yes") {
    #include "nondimension.H"
//    #include "nondimension2.H"
}

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();

volScalarField& p = thermo.p();
volScalarField& e = thermo.he();
volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
const volScalarField& mu = thermo.mu();

bool inviscid(true);
if (max(mu.internalField()) > 0.0)
{
    inviscid = false;
}

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
	"U",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

//#include "rhoBoundaryTypes.H"
volScalarField rho
(
    IOobject
    (
	"rho",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    thermo.rho()
//    rhoBoundaryTypes
);

volVectorField rhoU
(
    IOobject
    (
	"rhoU",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
	"rhoE",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);


surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh.Sf() & fvc::interpolate(rhoU)
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
