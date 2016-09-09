word ifLinear("No");
if (mesh.schemesDict().readIfPresent("ifLinear", ifLinear))
{
    if ((ifLinear == "Yes") || (ifLinear == "No"))
    {
        Info<< "If Running Linearized CFD: " << ifLinear << endl;
    }
    else
    {
        FatalErrorIn
        (
            "fsiFoamLinear::readIfLinear"
        )   << "ifLinear: " << ifLinear
            << " is not a valid choice. "
            << "Options are: Yes, No"
            << abort(FatalError);
    }
}
