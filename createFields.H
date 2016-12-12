/*double U_freestream_dim = 1;
double p_freestream_dim = 1;
double rho_freestream_dim = 1;
double T_freestream_dim = 1;
double cSqr_freestream_dim = 1;

#include "readIfNondim.H"

if (ifNondim == "Yes") {
    #include "nondimension.H"
//    #include "nondimension2.H"
}
*/

Info<< "Reading thermophysical properties\n" << endl;
/*
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
*/

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
	"p",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
	"T",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);


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
Info<< "end Reading fields\n" << endl;

#include "readThermo.H"
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
    p/(Rspec*T)
//    rhoBoundaryTypes
);

volScalarField psi
(
    IOobject
    (
	"psi",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    rho/p
);

volScalarField rPsi(1.0/psi);

volScalarField e
(
    IOobject
    (
	"e",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    Cv*T
//    h - p/rho
//    p/((gamma-1)*rho)
);

volScalarField h
(
    IOobject
    (
	"h",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    e + rPsi
//    Cp*T
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

//correct boundary condtions
//p, T, U, rho, psi, rPsi, e, h, rhoU, rhoE, 10 vars to correct
farFieldBoundaryCondition(rho, U, p);
T.boundaryField() = p.boundaryField()/(Rspec.boundaryField()*rho.boundaryField());
psi.boundaryField() = rho.boundaryField()/p.boundaryField();
rPsi.boundaryField() = p.boundaryField()/rho.boundaryField();
e.boundaryField() = Cv.boundaryField()*T.boundaryField();
h.boundaryField() = e.boundaryField()+rPsi.boundaryField();
rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();
rhoE.boundaryField() = rho.boundaryField()*(e.boundaryField() + 0.5*magSqr(U.boundaryField()));

/*Info << p << endl;
Info <<"p" << endl;
cin.get();
Info << T << endl;
Info <<"T" << endl;
cin.get();
Info << U << endl;
Info <<"U" << endl;
cin.get();
Info << rho << endl;
Info <<"rho" << endl;
cin.get();*/

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

/*Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);*/
