#include "constants.H"

volVectorField U_dim
(
    IOobject
    (
	"U",
	"0_dim",
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);

volScalarField p_dim
(
    IOobject
    (
	"p",
	"0_dim",
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);

volScalarField T_dim
(
    IOobject
    (
	"T",
	"0_dim",
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);

IOdictionary dict
(
    IOobject
    (
	"thermophysicalProperties",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);
scalar Cp = readScalar(dict.subDict("mixture").subDict("thermodynamics").lookup("Cp"));
scalar molWeight = readScalar(dict.subDict("mixture").subDict("specie").lookup("molWeight"));
scalar Rspec = constant::physicoChemical::R.value()*1000/molWeight; //(J/(kg*K))
scalar gamma = Cp/(Cp - Rspec);
forAll(mesh.boundary(), patchi) {
    if (mesh.boundary()[patchi].name() == "inlet") {
	U_freestream_dim = mag(U_dim.boundaryField()[patchi][0]);
	p_freestream_dim = p_dim.boundaryField()[patchi][0];
	T_freestream_dim = T_dim.boundaryField()[patchi][0];
	rho_freestream_dim = p_freestream_dim/Rspec/T_freestream_dim;
	cSqr_freestream_dim = gamma*Rspec*T_freestream_dim;
    }
}
volVectorField U
(
    IOobject
    (
	"U",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    U_dim/std::sqrt(cSqr_freestream_dim),
    U_dim.boundaryField().types()
);

volScalarField p
(
    IOobject
    (
	"p",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    p_dim/rho_freestream_dim/cSqr_freestream_dim,
    p_dim.boundaryField().types()
);

volScalarField T
(
    IOobject
    (
	"T",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    T_dim/(cSqr_freestream_dim),
    T_dim.boundaryField().types()
);

U.write();
p.write();
T.write();

Info << "Variables nondimensionized" << endl;
//cin.get();
