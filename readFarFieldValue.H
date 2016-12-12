IOdictionary dictIC
(
    IOobject
    (
	"ICs",
	"0", //mesh.time().timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);

vector U_inf = vector(dictIC.subDict("farField").lookup("velocity_inf"));
scalar rho_inf = readScalar(dictIC.subDict("farField").lookup("density_inf"));
scalar p_inf = readScalar(dictIC.subDict("farField").lookup("pressure_inf"));
/*Info << U_inf << endl;
Info << "U_inf" << endl;
Info << rho_inf << endl;
Info << "rho_inf" << endl;
Info << p_inf << endl;
Info << "p_inf" << endl;
cin.get();*/
