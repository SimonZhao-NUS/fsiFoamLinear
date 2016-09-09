word ifNondim("No");
if (mesh.schemesDict().readIfPresent("ifNondim", ifNondim))
{
    if ((ifNondim == "Yes") || (ifNondim == "No"))
    {
        Info<< "If Running Nondimensionized CFD: " << ifNondim << endl;
    }
    else
    {
        FatalErrorIn
        (
            "fsiFoamNondim::readIfNondim"
        )   << "ifNondim: " << ifNondim
            << " is not a valid choice. "
            << "Options are: Yes, No"
            << abort(FatalError);
    }
}
