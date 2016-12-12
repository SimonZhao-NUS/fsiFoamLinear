    pointIOField newPoints
    (
	IOobject
        (
            "newPoints",
	    runTime.timeName(),
	    "polyMesh",
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh.points()
    );

    scalar alpha = -.1/57.3;

    forAll(mesh.points(), pn) {
	scalar Lx = mesh.points()[pn].x() - .25;
	scalar Ly = mesh.points()[pn].y();
	scalar dx = -Ly*alpha;
	scalar dy = Lx*alpha;
	newPoints[pn].x() = mesh.points()[pn].x()+dx;
	newPoints[pn].y() = mesh.points()[pn].y()+dy;
    }
	
    mesh.movePoints(newPoints);
