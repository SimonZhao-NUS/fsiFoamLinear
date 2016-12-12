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

/*scalar Cp = readScalar(dict.subDict("mixture").subDict("thermodynamics").lookup("Cp"));
scalar molWeight = readScalar(dict.subDict("mixture").subDict("specie").lookup("molWeight"));
scalar Rspec_const = constant::physicoChemical::R.value()*1000/molWeight; //(J/(kg*K))
scalar Cv = Cp - Rspec_const;
scalar gamma_const = Cp/Cv;

volScalarField Rspec
{
    IOobject
    (
        "Rspec",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Rspec", dimless, Rspec_const)
};

volScalarField gamma
{
    IOobject
    (
        "mygamma",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("mygamma", dimless, gamma_const)
};*/

scalar molWeight = readScalar(dict.subDict("mixture").subDict("specie").lookup("molWeight"));

volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar
    (
	"Cp",
	dimensionSet(0, 2, -2, -1, 0, 0, 0),
	readScalar
	(
	    dict.subDict("mixture").subDict("thermodynamics").lookup("Cp")
	)
    )
);

volScalarField Rspec_const 
(
    IOobject
    (
        "molWeight",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar
    (
        "molWeight",
	dimensionSet(0, 2, -2, -1, 0, 0, 0),
	constant::physicoChemical::R.value()*1000/molWeight //(J/(kg*K))
    )
);

volScalarField Cv = Cp - Rspec_const;
volScalarField gamma_const = Cp/Cv;

volScalarField Rspec
(
    IOobject
    (
        "Rspec",
        runTime.timeName(),
        mesh
    ),
//    dimensionedScalar("Rspec", Rspec_const)
    Rspec_const
);

volScalarField gamma
(
    IOobject
    (
        "mygamma",
        runTime.timeName(),
        mesh
    ),
    gamma_const
);
