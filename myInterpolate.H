#include "vkSlopeLimiter.H"

template <class Type>
const GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
myGaussGrad
(
    const GeometricField<Type, fvPatchField, volMesh> vf
)
{
    const fvMesh& mesh = vf.mesh();
    GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh> sfSum
    (
	IOobject
	(
	    "mySurfaceSum",
	    vf.instance(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimensioned<typename outerProduct<vector, Type>::type>
	(
	    "my0",
	    vf.dimensions()/dimLength,
	    pTraits<typename outerProduct<vector, Type>::type>::zero
	),
	zeroGradientFvPatchField<typename outerProduct<vector, Type>::type>::typeName
    );

    Field<typename outerProduct<vector, Type>::type>& sfSumi = sfSum.internalField();
    const Field<Type>& vfi = vf.internalField();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const Field<vector>& sfi = mesh.Sf().internalField();
    const Field<vector>& cfi = mesh.Cf().internalField();
    const Field<vector>& ci = mesh.C().internalField();

    forAll(owner, facei)
    {
        scalar SfdOwn = mag(sfi[facei] & (cfi[facei] - ci[owner[facei]]));
        scalar SfdNei = mag(sfi[facei] & (ci[neighbour[facei]] - cfi[facei]));
        scalar w = SfdNei/(SfdOwn + SfdNei);

	Type vInterp = w*vfi[owner[facei]] + (1-w)*vfi[neighbour[facei]];
	typename outerProduct<vector, Type>::type sfiVInterp = sfi[facei]*vInterp;
	sfSumi[owner[facei]] += sfiVInterp;
	sfSumi[neighbour[facei]] -= sfiVInterp;

    }
    
    forAll(mesh.boundary(), patchi)
    {
	const labelUList& pFaceCells =
	mesh.boundary()[patchi].faceCells();

	const fvsPatchField<vector>& psf = mesh.Sf().boundaryField()[patchi];
	const fvPatchField<Type>& vfb = vf.boundaryField()[patchi];

	forAll(mesh.boundary()[patchi], facei)
	{
//	    Type vInterp = .5*(vfi[pFaceCells[facei]] + vfb[facei]);
	    Type vInterp =  vfb[facei];
	    typename outerProduct<vector, Type>::type sfbVInterp = psf[facei]*vInterp;
	    sfSumi[pFaceCells[facei]] += sfbVInterp;
	}
    }
    sfSumi /= mesh.V();
    sfSum.correctBoundaryConditions();

    return sfSum;
}

scalar
eleWiseProd
(
    const scalar sca1,
    const scalar sca2
)
{ 
    scalar result(sca1*sca2);
    return result;
}

vector
eleWiseProd
(
    const vector vec1,
    const vector vec2
)
{ 
    vector result
    (
	vec1.x()*vec2.x(),
       	vec1.y()*vec2.y(),
       	vec1.z()*vec2.z()
    );

    return result;
}

const std::vector<std::vector<unsigned int> >
createCellNeigList(const fvMesh& mesh)
{
    std::vector<std::vector<unsigned int> > cellNeiList
    (
	mesh.nCells()
    );
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    forAll(P, fi) {
	cellNeiList[P[fi]].push_back(N[fi]);
	cellNeiList[N[fi]].push_back(P[fi]);
    }
    
    return cellNeiList;
}



Field<scalar> epsilonFunc(const fvMesh &mesh, const scalar K)
{
    Field<scalar> epsilon
    (
	mesh.V().size(),
	0
    );

    return epsilon;
}
Field<scalar> vfiMaxFunc
(
    const Field<scalar>& vfi,
    const fvMesh& mesh
)
{
    Field<scalar> vfiMax = vfi;

    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    forAll(P, fi) {
	scalar vf_n = vfi[N[fi]];
	if (vfiMax[P[fi]] < vf_n) {
	    vfiMax[P[fi]] = vf_n;
	}

	vf_n = vfi[P[fi]];
	if (vfiMax[N[fi]] < vf_n) {
	    vfiMax[N[fi]] = vf_n;
	}
    }

    return vfiMax;
}

Field<scalar> vfiMinFunc
(
    const Field<scalar>& vfi,
    const fvMesh& mesh
)
{
    Field<scalar> vfiMin = vfi;

    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    forAll(P, fi) {
	scalar vf_n = vfi[N[fi]];
	if (vfiMin[P[fi]] > vf_n) {
	    vfiMin[P[fi]] = vf_n;
	}

	vf_n = vfi[P[fi]];
	if (vfiMin[N[fi]] > vf_n) {
	    vfiMin[N[fi]] = vf_n;
	}
    }

    return vfiMin;
}

scalar limiterFunc
(
    const scalar delta1Max,
    const scalar delta1Min,
    const scalar delta2,
    const scalar epsilon
)
{
    scalar limiter = 1;

    if(delta2 > 0) {
	limiter = 
	(
	    (
		(
		    magSqr(delta1Max) +
		    magSqr(epsilon)
		) + 
		2*delta2*delta1Max
	    ) /
	    (
		magSqr(delta1Max) +
		2*magSqr(delta2) +
		delta1Max*delta2 +
		magSqr(epsilon)
	    )
	);
	
    } else if (delta2 < 0) {
	limiter = 
	(
	    (
		(
		    magSqr(delta1Min) +
		    magSqr(epsilon)
		) + 
		2*delta2*delta1Min
	    ) /
	    (
		magSqr(delta1Min) +
		2*magSqr(delta2) +
		delta1Min*delta2 +
		magSqr(epsilon)
	    )
	);
    } else {
	limiter = 1;
    }
    
    return limiter;
}

Field<scalar>
myLimiter
(
    const Field<vector>& fgradi,
    const volScalarField& vf
)
{
    Field<scalar> limiter(vf.internalField().size(), 1);
    const fvMesh& mesh = vf.mesh();

    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    const surfaceVectorField &faceCenter = mesh.Cf();
    const volVectorField &cellCenter = mesh.C();

    const Field<scalar>& vfi = vf.internalField();
//    const Field<scalar> epsilonField = epsilonFunc(mesh, 0);
    const Field<scalar> vfiMax = vfiMaxFunc(vfi, mesh);
    const Field<scalar> vfiMin = vfiMinFunc(vfi, mesh);

/*    Info << fgradi<< endl;
    cin.get();
    Info << vfiMax<< endl;
    cin.get();
    Info << vfiMin<< endl;
    cin.get();*/
    forAll(P, fi)
    {
	label p = P[fi];
	label n = N[fi];
	
	vector drt = faceCenter[fi] - cellCenter[p];
	scalar delta1Max = vfiMax[p] - vfi[p];
	scalar delta1Min = vfiMin[p] - vfi[p];
	scalar delta2 = fgradi[p]&drt;
	scalar epsilon = 0;

	scalar limiter_i = 
	limiterFunc
	(
	    delta1Max,
	    delta1Min,
	    delta2,
	    epsilon
	);
	if (limiter[p] > limiter_i) {
	    limiter[p] = limiter_i;
	}

	drt = faceCenter[fi] - cellCenter[n];
	delta1Max = vfiMax[n] - vfi[n];
	delta1Min = vfiMin[n] - vfi[n];
	delta2 = fgradi[n]&drt;
	epsilon = 0;

	limiter_i = 
	limiterFunc
	(
	    delta1Max,
	    delta1Min,
	    delta2,
	    epsilon
	);
	if (limiter[n] > limiter_i) {
	    limiter[n] = limiter_i;
	}
    }

    return limiter;
}
Field<vector>
myLimiter
(
    const Field<typename outerProduct<vector, vector>::type>& fgradi,
    const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
    Field<vector> limiter(vf.internalField().size(), vector(1, 1, 1));
    const fvMesh& mesh = vf.mesh();

    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    const surfaceVectorField &faceCenter = mesh.Cf();
    const volVectorField &cellCenter = mesh.C();
    const Field<vector>& vfi = vf.internalField();
    
    for(direction cmpt = 0; cmpt < 3; ++cmpt) {
	const Field<scalar> vfiMax = vfiMaxFunc(vfi.component(cmpt), mesh);
	const Field<scalar> vfiMin = vfiMinFunc(vfi.component(cmpt), mesh);

/*    Info << fgradi<< endl;
    cin.get();
    Info << vfiMax<< endl;
    cin.get();
    Info << vfiMin<< endl;
    cin.get();*/
	forAll(P, fi)
	{
	    label p = P[fi];
	    label n = N[fi];
	    
	    vector drt = faceCenter[fi] - cellCenter[p];
	    scalar delta1Max = vfiMax[p] - vfi[p].component(cmpt);
	    scalar delta1Min = vfiMin[p] - vfi[p].component(cmpt);
	    scalar delta2 = (fgradi[p]&drt).component(cmpt);
	    scalar epsilon = 0;

	    scalar limiter_i = 
	    limiterFunc
	    (
		delta1Max,
		delta1Min,
		delta2,
		epsilon
	    );
	    if (limiter[p].component(cmpt) > limiter_i) {
		limiter[p].component(cmpt) = limiter_i;
	    }

	    drt = faceCenter[fi] - cellCenter[n];
	    delta1Max = vfiMax[n] - vfi[n].component(cmpt);
	    delta1Min = vfiMin[n] - vfi[n].component(cmpt);
	    delta2 = (fgradi[n]&drt).component(cmpt);
	    epsilon = 0;

	    limiter_i = 
	    limiterFunc
	    (
		delta1Max,
		delta1Min,
		delta2,
		epsilon
	    );
	    if (limiter[n].component(cmpt) > limiter_i) {
		limiter[n].component(cmpt) = limiter_i;
	    }
	}
    }

    return limiter;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
myInterpolatePos
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const Field<Type>& vfi = vf.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    Field<Type>& sfi = sf.internalField();

    const surfaceVectorField &faceCenter = mesh.Cf();
    const volVectorField &cellCenter = mesh.C();
    const GeometricField <typename outerProduct<vector, Type>::type, fvPatchField, volMesh> fgrad = fvc::grad(vf, "Gauss linear");
//    const GeometricField <typename outerProduct<vector, Type>::type, fvPatchField, volMesh> fgrad(myGaussGrad(vf));
    const Field<typename outerProduct<vector, Type>::type> &fgradi = fgrad.internalField();
    
/*    Info << vf << endl;
    Info << vf.name() << endl;
    cin.get();*/
/*    Info << fgrad << endl;
    Info << vf.name() << endl;
    cin.get();*/
//    szhaoFoam::slopeLimiter<Type> limiter(P.size(), Foam::one());
    szhaoFoam::vkSlopeLimiter<Type> limiter(vf, fgrad, 20);
/*    Info << limiter << endl;
    Info << vf.name() << endl;
    cin.get();*/
    for (register label fi=0; fi<P.size(); fi++)
    {
	vector drt = faceCenter[fi] - cellCenter[P[fi]];
//        sfi[fi] = vfi[P[fi]] + (drt&fgradi[P[fi]]); //2nd order reconstruction w/o limiter
//        sfi[fi] = vfi[P[fi]] + eleWiseProd(limiter[P[fi]], drt&fgradi[P[fi]]); //2nd order reconstruction w limiter
        sfi[fi] = vfi[P[fi]]; //1st order approx
    }

    forAll(vf.boundaryField(), pi)
    {
	sf.boundaryField()[pi] = vf.boundaryField()[pi];
    }

    return tsf;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
myInterpolateNeg
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const Field<Type>& vfi = vf.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    Field<Type>& sfi = sf.internalField();

    const surfaceVectorField &faceCenter = mesh.Cf();
    const volVectorField &cellCenter = mesh.C();
    const GeometricField <typename outerProduct<vector, Type>::type, fvPatchField, volMesh> fgrad = fvc::grad(vf, "Gauss linear");
//    const GeometricField <typename outerProduct<vector, Type>::type, fvPatchField, volMesh> fgrad(myGaussGrad(vf));
    const Field<typename outerProduct<vector, Type>::type> &fgradi = fgrad.internalField();
    szhaoFoam::vkSlopeLimiter<Type> limiter(vf, fgrad, 20);
    for (register label fi=0; fi<N.size(); fi++)
    {
	vector drt = faceCenter[fi] - cellCenter[N[fi]];
//        sfi[fi] = vfi[N[fi]] + (drt&fgradi[N[fi]]);
//        sfi[fi] = vfi[N[fi]] + eleWiseProd(limiter[N[fi]], drt&fgradi[N[fi]]); //2nd order reconstruction w limiter
        sfi[fi] = vfi[N[fi]];
    }

    forAll(vf.boundaryField(), pi)
    {
	sf.boundaryField()[pi] = vf.boundaryField()[pi];
    }

    return tsf;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
myInterpolatePos
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField <typename outerProduct<vector, Type>::type, fvPatchField, volMesh> fgrad,
    szhaoFoam::vkSlopeLimiter<Type> limiter
)
{
    const Field<Type>& vfi = vf.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();
    Field<Type>& sfi = sf.internalField();

    const surfaceVectorField &faceCenter = mesh.Cf();
    const volVectorField &cellCenter = mesh.C();
    
    const Field<typename outerProduct<vector, Type>::type> &fgradi = fgrad.internalField();
    
    for (register label fi=0; fi<P.size(); fi++)
    {
	vector drt = faceCenter[fi] - cellCenter[P[fi]];
//        sfi[fi] = vfi[P[fi]] + (drt&fgradi[P[fi]]); //2nd order reconstruction w/o limiter
        sfi[fi] = vfi[P[fi]] + eleWiseProd(limiter[P[fi]], drt&fgradi[P[fi]]); //2nd order reconstruction w limiter
//        sfi[fi] = vfi[P[fi]]; //1st order approx
    }

    forAll(vf.boundaryField(), pi)
    {
	sf.boundaryField()[pi] = vf.boundaryField()[pi];
    }

    return tsf;
}

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
myInterpolateNeg
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField <typename outerProduct<vector, Type>::type, fvPatchField, volMesh> fgrad,
    szhaoFoam::vkSlopeLimiter<Type> limiter
)
{
    const Field<Type>& vfi = vf.internalField();

    const fvMesh& mesh = vf.mesh();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();
    Field<Type>& sfi = sf.internalField();

    const surfaceVectorField &faceCenter = mesh.Cf();
    const volVectorField &cellCenter = mesh.C();
    
    const Field<typename outerProduct<vector, Type>::type> &fgradi = fgrad.internalField();
    for (register label fi=0; fi<N.size(); fi++)
    {
	vector drt = faceCenter[fi] - cellCenter[N[fi]];
//        sfi[fi] = vfi[N[fi]] + (drt&fgradi[N[fi]]);
        sfi[fi] = vfi[N[fi]] + eleWiseProd(limiter[N[fi]], drt&fgradi[N[fi]]); //2nd order reconstruction w limiter
//        sfi[fi] = vfi[N[fi]];
    }

    forAll(vf.boundaryField(), pi)
    {
	sf.boundaryField()[pi] = vf.boundaryField()[pi];
    }

    return tsf;
}
