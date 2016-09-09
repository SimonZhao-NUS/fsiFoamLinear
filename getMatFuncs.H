void createMap
(
    const fvMesh& mesh,
    myMat& map    
)
{
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    forAll(P, fi) {
	map[P[fi]][map[P[fi]][0]+1] = fi;
	map[N[fi]][map[N[fi]][0]+1] = fi;
	
	++map[P[fi]][0];
	++map[N[fi]][0];	
    }
/*    myPrintMat(map);
    Info << map[149][1] << endl;
    Info << map[149][2] << endl;
    Info << map[149][3] << endl;
    Info << map[149][4] << endl;
    cin.get();*/
}
    
void getMatTSp
(
    SpMat& matT_pos,
    SpMat& matT_neg,
    SpMat& matT2_pos,
    SpMat& matT2_neg,
    SpMat& matT3_pos,
    SpMat& matT3_neg,
    Field<scalar> limiterRho_pos,
    Field<scalar> limiterRho_neg,
    Field<vector> limiterRhoU_pos,
    Field<vector> limiterRhoU_neg,
    Field<scalar> limiterRhoE_pos,
    Field<scalar> limiterRhoE_neg,
    const fvMesh& mesh
)
{
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    const Field<scalar>& vi = mesh.V();
    const Field<scalar>& magSfi = mesh.magSf().internalField();
    const Field<vector>& sfi = mesh.Sf().internalField();
    const Field<vector>& cfi = mesh.Cf().internalField();
    const Field<vector>& ci = mesh.C().internalField();

    TripletList tList_pos, tList_neg, t2List_pos, t2List_neg, t3List_pos, t3List_neg;
    tList_pos.reserve(P.size()*5);
    tList_pos.reserve(P.size()*5);
    t2List_neg.reserve(P.size()*10*4);
    t2List_neg.reserve(P.size()*10*4);
    int nBdFaces = 0;
    forAll(mesh.boundary(), patchi) {
	forAll(mesh.boundary()[patchi], fi) {
	    nBdFaces++;
	}
    }
    t3List_pos.reserve(nBdFaces);
    t3List_neg.reserve(nBdFaces);

    std::vector<std::vector<double> > map
    (
	mesh.nCells(),
	std::vector<double>
	(
	    7,
	    0
	)
    );

    createMap(mesh, map);

    for (int fi=0; fi<P.size(); fi++) {
	vector dp = cfi[fi] - ci[P[fi]];
	for (int i = 0; i < map[P[fi]][0] ; ++i) {
	    int sf = map[P[fi]][i+1];
	    scalar pn = 0;
	    if (mesh.owner()[sf] == P[fi]) {
		pn = 1;
	    } else {
		pn = -1;
	    }
	    scalar SfdOwn = mag(sfi[sf] & (cfi[sf] - ci[P[sf]]));
	    scalar SfdNei = mag(sfi[sf] & (ci[N[sf]] - cfi[sf]));
	    scalar w = SfdNei/(SfdOwn + SfdNei);
	    vector ri = w*sfi[sf]*pn;
	    vector li = (1-w)*sfi[sf]*pn;
	    scalar r = (dp&ri)/mesh.V()[P[fi]];
	    scalar l = (dp&li)/mesh.V()[P[fi]];
	    
/*	    if(fi == 300 && sf==224) {
		Info << r << endl;
		Info << l << endl;
		cin.get();
	    }*/
/*	    t2List_pos.push_back(Triplet(fi*5+0, P[sf]*5+0, r*limiterRho_pos[fi]));
	    t2List_pos.push_back(Triplet(fi*5+0, N[sf]*5+0, l*limiterRho_pos[fi]));
	    t2List_pos.push_back(Triplet(fi*5+1, P[sf]*5+1, r*limiterRhoU_pos[fi].x()));
	    t2List_pos.push_back(Triplet(fi*5+1, N[sf]*5+1, l*limiterRhoU_pos[fi].x()));
	    t2List_pos.push_back(Triplet(fi*5+2, P[sf]*5+2, r*limiterRhoU_pos[fi].y()));
	    t2List_pos.push_back(Triplet(fi*5+2, N[sf]*5+2, l*limiterRhoU_pos[fi].y()));
	    t2List_pos.push_back(Triplet(fi*5+3, P[sf]*5+3, r*limiterRhoU_pos[fi].z()));
	    t2List_pos.push_back(Triplet(fi*5+3, N[sf]*5+3, l*limiterRhoU_pos[fi].z()));
	    t2List_pos.push_back(Triplet(fi*5+4, P[sf]*5+4, r*limiterRhoE_pos[fi]));
	    t2List_pos.push_back(Triplet(fi*5+4, N[sf]*5+4, l*limiterRhoE_pos[fi]));*/
	    t2List_pos.push_back(Triplet(fi*5+0, P[sf]*5+0, r));
	    t2List_pos.push_back(Triplet(fi*5+0, N[sf]*5+0, l));
	    t2List_pos.push_back(Triplet(fi*5+1, P[sf]*5+1, r));
	    t2List_pos.push_back(Triplet(fi*5+1, N[sf]*5+1, l));
	    t2List_pos.push_back(Triplet(fi*5+2, P[sf]*5+2, r));
	    t2List_pos.push_back(Triplet(fi*5+2, N[sf]*5+2, l));
	    t2List_pos.push_back(Triplet(fi*5+3, P[sf]*5+3, r));
	    t2List_pos.push_back(Triplet(fi*5+3, N[sf]*5+3, l));
	    t2List_pos.push_back(Triplet(fi*5+4, P[sf]*5+4, r));
	    t2List_pos.push_back(Triplet(fi*5+4, N[sf]*5+4, l));
/*	    for (unsigned int q = 0; q < 5; q++) { 
		t2List_pos.push_back(Triplet(fi*5+q, P[sf]*5+q, r));
		t2List_pos.push_back(Triplet(fi*5+q, N[sf]*5+q, l));
	    }*/
	}

	vector dn = cfi[fi] - ci[N[fi]];
	for (int i = 0; i < map[N[fi]][0] ; ++i) {
	    int sf = map[N[fi]][i+1];
	    scalar pn = 0;
	    if (mesh.owner()[sf] == N[fi]) {
		pn = 1;
	    } else {
		pn = -1;
	    }
	    scalar SfdOwn = mag(sfi[sf] & (cfi[sf] - ci[P[sf]]));
	    scalar SfdNei = mag(sfi[sf] & (ci[N[sf]] - cfi[sf]));
	    scalar w = SfdNei/(SfdOwn + SfdNei);
	    vector ri = w*sfi[sf]*pn;
	    vector li = (1-w)*sfi[sf]*pn;
	    scalar r = (dn&ri)/mesh.V()[N[fi]];
	    scalar l = (dn&li)/mesh.V()[N[fi]];
	    
/*	    t2List_neg.push_back(Triplet(fi*5+0, P[sf]*5+0, r*limiterRho_neg[fi]));
	    t2List_neg.push_back(Triplet(fi*5+0, N[sf]*5+0, l*limiterRho_neg[fi]));
	    t2List_neg.push_back(Triplet(fi*5+1, P[sf]*5+1, r*limiterRhoU_neg[fi].x()));
	    t2List_neg.push_back(Triplet(fi*5+1, N[sf]*5+1, l*limiterRhoU_neg[fi].x()));
	    t2List_neg.push_back(Triplet(fi*5+2, P[sf]*5+2, r*limiterRhoU_neg[fi].y()));
	    t2List_neg.push_back(Triplet(fi*5+2, N[sf]*5+2, l*limiterRhoU_neg[fi].y()));
	    t2List_neg.push_back(Triplet(fi*5+3, P[sf]*5+3, r*limiterRhoU_neg[fi].z()));
	    t2List_neg.push_back(Triplet(fi*5+3, N[sf]*5+3, l*limiterRhoU_neg[fi].z()));
	    t2List_neg.push_back(Triplet(fi*5+4, P[sf]*5+4, r*limiterRhoE_neg[fi]));
	    t2List_neg.push_back(Triplet(fi*5+4, N[sf]*5+4, l*limiterRhoE_neg[fi]));*/
	    t2List_neg.push_back(Triplet(fi*5+0, P[sf]*5+0, r));
	    t2List_neg.push_back(Triplet(fi*5+0, N[sf]*5+0, l));
	    t2List_neg.push_back(Triplet(fi*5+1, P[sf]*5+1, r));
	    t2List_neg.push_back(Triplet(fi*5+1, N[sf]*5+1, l));
	    t2List_neg.push_back(Triplet(fi*5+2, P[sf]*5+2, r));
	    t2List_neg.push_back(Triplet(fi*5+2, N[sf]*5+2, l));
	    t2List_neg.push_back(Triplet(fi*5+3, P[sf]*5+3, r));
	    t2List_neg.push_back(Triplet(fi*5+3, N[sf]*5+3, l));
	    t2List_neg.push_back(Triplet(fi*5+4, P[sf]*5+4, r));
	    t2List_neg.push_back(Triplet(fi*5+4, N[sf]*5+4, l));
/*	    for (unsigned int q = 0; q < 5; q++) { 
		t2List_neg.push_back(Triplet(fi*5+q, P[sf]*5+q, r));
		t2List_neg.push_back(Triplet(fi*5+q, N[sf]*5+q, l));
	    }*/
	}

	for (unsigned int q = 0; q < 5; q++) { 
	    tList_pos.push_back(Triplet(fi*5+q, P[fi]*5+q, 1.0));
	    tList_neg.push_back(Triplet(fi*5+q, N[fi]*5+q, 1.0));
	}
    }
    
    int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	const fvsPatchField<vector>& psf = mesh.Sf().boundaryField()[patchi];
	const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
	forAll(mesh.boundary()[patchi], fi) {
/*		if (mesh.boundary()[patchi].name() == "wing") {
		Info << mesh.Cf().boundaryField()[patchi][fi] << endl;
		Info << mesh.Sf().boundaryField()[patchi][fi] << endl;
		cin.get();
	    }*/
		
	    for (int i = 0; i < map[faceCells[fi]][0] ; ++i) {
		int sf = map[faceCells[fi]][i+1];
		vector d = cfi[sf] - ci[faceCells[fi]];
		vector ri = psf[fi];
		scalar r = (d&ri)/mesh.V()[faceCells[fi]];
	    /*    if (mesh.owner()[sf] == faceCells[fi]) {
		    t3List_pos.push_back(Triplet(sf*5+0, nBdIndex*5+0, r*limiterRho_pos[sf]));
		    t3List_pos.push_back(Triplet(sf*5+1, nBdIndex*5+1, r*limiterRhoU_pos[sf].x()));
		    t3List_pos.push_back(Triplet(sf*5+2, nBdIndex*5+2, r*limiterRhoU_pos[sf].y()));
		    t3List_pos.push_back(Triplet(sf*5+3, nBdIndex*5+3, r*limiterRhoU_pos[sf].z()));
		    t3List_pos.push_back(Triplet(sf*5+4, nBdIndex*5+4, r*limiterRhoE_pos[sf]));
		} else {
		    t3List_neg.push_back(Triplet(sf*5+0, nBdIndex*5+0, r*limiterRho_neg[sf]));
		    t3List_neg.push_back(Triplet(sf*5+1, nBdIndex*5+1, r*limiterRhoU_neg[sf].x()));
		    t3List_neg.push_back(Triplet(sf*5+2, nBdIndex*5+2, r*limiterRhoU_neg[sf].y()));
		    t3List_neg.push_back(Triplet(sf*5+3, nBdIndex*5+3, r*limiterRhoU_neg[sf].z()));
		    t3List_neg.push_back(Triplet(sf*5+4, nBdIndex*5+4, r*limiterRhoE_neg[sf]));
		}*/
		if (mesh.owner()[sf] == faceCells[fi]) {
		    t3List_pos.push_back(Triplet(sf*5+0, nBdIndex*5+0, r));
		    t3List_pos.push_back(Triplet(sf*5+1, nBdIndex*5+1, r));
		    t3List_pos.push_back(Triplet(sf*5+2, nBdIndex*5+2, r));
		    t3List_pos.push_back(Triplet(sf*5+3, nBdIndex*5+3, r));
		    t3List_pos.push_back(Triplet(sf*5+4, nBdIndex*5+4, r));
		} else {
		    t3List_neg.push_back(Triplet(sf*5+0, nBdIndex*5+0, r));
		    t3List_neg.push_back(Triplet(sf*5+1, nBdIndex*5+1, r));
		    t3List_neg.push_back(Triplet(sf*5+2, nBdIndex*5+2, r));
		    t3List_neg.push_back(Triplet(sf*5+3, nBdIndex*5+3, r));
		    t3List_neg.push_back(Triplet(sf*5+4, nBdIndex*5+4, r));
		}
	    }
	    nBdIndex++;
	}
    }

    matT_pos.setFromTriplets(tList_pos.begin(), tList_pos.end());
    matT_neg.setFromTriplets(tList_neg.begin(), tList_neg.end());
    matT2_pos.setFromTriplets(t2List_pos.begin(), t2List_pos.end());
    matT2_neg.setFromTriplets(t2List_neg.begin(), t2List_neg.end());
    matT3_pos.setFromTriplets(t3List_pos.begin(), t3List_pos.end());
    matT3_neg.setFromTriplets(t3List_neg.begin(), t3List_neg.end());
}

void getMatASp
(
    SpMat& matA_pos,
    SpMat& matA_neg,
    SpMat& matA_abs,
    const surfaceVectorField& U_pos,
    const surfaceScalarField& H_pos,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& H_neg,
    const surfaceVectorField& U_tilt,
    const surfaceScalarField& H_tilt,
    const surfaceScalarField& gamma,
    const fvMesh& mesh
)
{
    const labelUList& P = mesh.owner();
    TripletList tList_pos, tList_neg, tList_abs;
    tList_pos.reserve(P.size()*5);
    tList_neg.reserve(P.size()*5);
    tList_abs.reserve(P.size()*5);
    
    for (int fi=0; fi<P.size(); fi++) {
	myMat jcb_pos
	(
	    innProd
	    (
		innProd
		(
		    calcEigVecR
		    (
			mesh.Sf().internalField()[fi],
			U_pos.internalField()[fi],
			H_pos.internalField()[fi],
			gamma.internalField()[fi]
		    ),
		    calcEigVal
		    (
			mesh.Sf().internalField()[fi],
			U_pos.internalField()[fi],
			H_pos.internalField()[fi],
			gamma.internalField()[fi]
		    )
		),
		calcEigVecL
		(
		    mesh.Sf().internalField()[fi],
		    U_pos.internalField()[fi],
		    H_pos.internalField()[fi],
		    gamma.internalField()[fi]
		)
	    )
	);

	myMat jcb_neg
	(
	    innProd
	    (
		innProd
		(
		    calcEigVecR
		    (
			mesh.Sf().internalField()[fi],
			U_neg.internalField()[fi],
			H_neg.internalField()[fi],
			gamma.internalField()[fi]
		    ),
		    calcEigVal
		    (
			mesh.Sf().internalField()[fi],
			U_neg.internalField()[fi],
			H_neg.internalField()[fi],
			gamma.internalField()[fi]
		    )
		),
		calcEigVecL
		(
		    mesh.Sf().internalField()[fi],
		    U_neg.internalField()[fi],
		    H_neg.internalField()[fi],
		    gamma.internalField()[fi]
		)
	    )
	);
	myMat jcb_abs
	(
	    innProd
	    (
		innProd
		(
		    calcEigVecR
		    (
			mesh.Sf().internalField()[fi],
			U_tilt.internalField()[fi],
			H_tilt.internalField()[fi],
			gamma.internalField()[fi]
		    ),
		    calcEigVal_abs
		    (
			mesh.Sf().internalField()[fi],
			U_tilt.internalField()[fi],
			H_tilt.internalField()[fi],
			gamma.internalField()[fi]
		    )
		),
		calcEigVecL
		(
		    mesh.Sf().internalField()[fi],
		    U_tilt.internalField()[fi],
		    H_tilt.internalField()[fi],
		    gamma.internalField()[fi]
		)
	    )
	);
	
	for(int nr = 0; nr < 5; nr++) {
	    for(int nc = 0; nc < 5; nc++) {
		tList_pos.push_back(Triplet(fi*5+nr, fi*5+nc, jcb_pos[nr][nc]));
		tList_neg.push_back(Triplet(fi*5+nr, fi*5+nc, jcb_neg[nr][nc]));
		tList_abs.push_back(Triplet(fi*5+nr, fi*5+nc, jcb_abs[nr][nc]));
	    }
	}
    }
    
    matA_pos.setFromTriplets(tList_pos.begin(), tList_pos.end());
    matA_neg.setFromTriplets(tList_neg.begin(), tList_neg.end());
    matA_abs.setFromTriplets(tList_abs.begin(), tList_abs.end());
}

myMat con2phy
(
    const scalar& rho,
    const vector& U,
    const scalar& T,
    const scalar& p,
    const scalar& gamma
)
{ //this function is to convert rho, rhoU, rhoE to rho, U, p
	myMat subD
	(
	    5,
	    std::vector<double>(5,0)
	);

	double u = U.x();
	double v = U.y();
	double w = U.z();
	double q = u*u + v*v + w*w;
	//row1
	subD[0][0] = 1;
	subD[0][1] = 0;
	subD[0][2] = 0;
	subD[0][3] = 0;
	subD[0][4] = 0;
	//row2
	subD[1][0] = -1.*u/rho;
	subD[1][1] = 1./rho;
	subD[1][2] = 0;
	subD[1][3] = 0;
	subD[1][4] = 0;
	//row3
	subD[2][0] = -1.*v/rho;
	subD[2][1] = 0;
	subD[2][2] = 1./rho;
	subD[2][3] = 0;
	subD[2][4] = 0;
	//row4
	subD[3][0] = -1.*w/rho;
	subD[3][1] = 0;
	subD[3][2] = 0;
	subD[3][3] = 1./rho;
	subD[3][4] = 0;
	//row5
	subD[4][0] = .5*(gamma-1)*q;
	subD[4][1] = -(gamma-1)*u;
	subD[4][2] = -(gamma-1)*v;
	subD[4][3] = -(gamma-1)*w;
	subD[4][4] = (gamma-1);      

	return subD;
}

myMat phy2con
(
    const scalar& rho,
    const vector& U,
    const scalar& T,
    const scalar& p,
    const scalar& gamma
)
{ //this function is to convert rho, U, p to rho, rhoU, rhoE 
	myMat subD
	(
	    5,
	    std::vector<double>(5,0)
	);

	double u = U.x();
	double v = U.y();
	double w = U.z();
	double q = u*u + v*v + w*w;
	//row1
	subD[0][0] = 1;
	subD[0][1] = 0;
	subD[0][2] = 0;
	subD[0][3] = 0;
	subD[0][4] = 0;
	//row2
	subD[1][0] = u;
	subD[1][1] = rho;
	subD[1][2] = 0;
	subD[1][3] = 0;
	subD[1][4] = 0;
	//row3
	subD[2][0] = v;
	subD[2][1] = 0;
	subD[2][2] = rho;
	subD[2][3] = 0;
	subD[2][4] = 0;
	//row4
	subD[3][0] = w;
	subD[3][1] = 0;
	subD[3][2] = 0;
	subD[3][3] = rho;
	subD[3][4] = 0;
	//row5
	subD[4][0] = .5*q;
	subD[4][1] = rho*u;
	subD[4][2] = rho*v;
	subD[4][3] = rho*w;
	subD[4][4] = 1/(gamma-1);

	return subD;
}

void getMatCSp
(
    SpMat& matC,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& p,
    const volScalarField& gamma,
    const unsigned int nBdFaces,
    const fvMesh& mesh
)
{
    TripletList tList;
    tList.reserve(nBdFaces*5);
    
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	forAll(mesh.boundary()[patchi], facei) {
	    myMat mat
	    (
		con2phy
		(
		    rho.boundaryField()[patchi][facei],
		    U.boundaryField()[patchi][facei],
		    T.boundaryField()[patchi][facei],
		    p.boundaryField()[patchi][facei],
		    gamma.boundaryField()[patchi][facei]
		)
	    );

	    myMat mat0
	    (
		phy2con
		(
		    rho.boundaryField()[patchi][facei],
		    U.boundaryField()[patchi][facei],
		    T.boundaryField()[patchi][facei],
		    p.boundaryField()[patchi][facei],
		    gamma.boundaryField()[patchi][facei]
		)
	    );
//	    myPrintMat(mat);
//	    myPrintMat(mat0);
//	    myPrintMat(innProd(mat,mat0));
//	    cin.get();
	    for(int nr = 0; nr < 5; nr++) {
		for(int nc = 0; nc < 5; nc++) {
		    tList.push_back(Triplet(nBdIndex*5+nr, nBdIndex*5+nc, mat[nr][nc]));
		}
	    }
	    ++nBdIndex;
	}
    }
    
    matC.setFromTriplets(tList.begin(), tList.end());
}

void getMatCinvSp
(
    SpMat& matC,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& p,
    const volScalarField& gamma,
    const unsigned int nBdFaces,
    const fvMesh& mesh
)
{
    TripletList tList;
    tList.reserve(nBdFaces*5);
    
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	forAll(mesh.boundary()[patchi], facei) {
	    myMat mat
	    (
		phy2con
		(
		    rho.boundaryField()[patchi][facei],
		    U.boundaryField()[patchi][facei],
		    T.boundaryField()[patchi][facei],
		    p.boundaryField()[patchi][facei],
		    gamma.boundaryField()[patchi][facei]
		)
	    );
	    for(int nr = 0; nr < 5; nr++) {
		for(int nc = 0; nc < 5; nc++) {
		    tList.push_back(Triplet(nBdIndex*5+nr, nBdIndex*5+nc, mat[nr][nc]));
		}
	    }
	    ++nBdIndex;
	}
    }
    
    matC.setFromTriplets(tList.begin(), tList.end());
}

void getMatDSp
(
    SpMat& matD,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& p,
    const volScalarField& gamma,
    const unsigned int nIntCells
)
{
    TripletList tList;
    tList.reserve(nIntCells*5);

    for (unsigned int ni=0; ni<nIntCells; ni++) {
	myMat mat
	(
	    con2phy
	    (
		rho.internalField()[ni],
		U.internalField()[ni],
		T.internalField()[ni],
		p.internalField()[ni],
		gamma.internalField()[ni]
	    )
	);
	for(int nr = 0; nr < 5; nr++) {
	    for(int nc = 0; nc < 5; nc++) {
		tList.push_back(Triplet(ni*5+nr, ni*5+nc, mat[nr][nc]));
	    }
	}
    }
    
    matD.setFromTriplets(tList.begin(), tList.end());
}

void getMatBSp
(
    SpMat& matB,
    const volVectorField& U,
    const volScalarField& H,
    const volScalarField& gamma,
    const unsigned int nBdFaces,
    const fvMesh& mesh
)
{
    TripletList tList;
    tList.reserve(nBdFaces*5);
    
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	forAll(mesh.boundary()[patchi], facei) {
	    myMat mat
	    (
		innProd
		(
		    innProd
		    (
			calcEigVecR
			(
			    mesh.Sf().boundaryField()[patchi][facei],
			    U.boundaryField()[patchi][facei],
			    H.boundaryField()[patchi][facei],
			    gamma.boundaryField()[patchi][facei]
			),
			calcEigVal
			(
			    mesh.Sf().boundaryField()[patchi][facei],
			    U.boundaryField()[patchi][facei],
			    H.boundaryField()[patchi][facei],
			    gamma.boundaryField()[patchi][facei]
			)
		    ),
		    calcEigVecL
		    (
			mesh.Sf().boundaryField()[patchi][facei],
			U.boundaryField()[patchi][facei],
			H.boundaryField()[patchi][facei],
			gamma.boundaryField()[patchi][facei]
		    )
		)
	    );
	    for(int nr = 0; nr < 5; nr++) {
		for(int nc = 0; nc < 5; nc++) {
		    tList.push_back(Triplet(nBdIndex*5+nr, nBdIndex*5+nc, mat[nr][nc]));
		}
	    }
	    ++nBdIndex;
	}
    }
    
    matB.setFromTriplets(tList.begin(), tList.end());
}

void getMatBfSp
(
    SpMat& matBf,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p,
    const unsigned int nBdFaces,
    const fvMesh& mesh
)
{
    TripletList tList;
    tList.reserve(nBdFaces*5);

    scalar x0 = 0.25; 
    
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	const fvsPatchField<vector>& SfPatch = mesh.Sf().boundaryField()[patchi];
	const fvsPatchField<vector>& CfPatch = mesh.Cf().boundaryField()[patchi];
	const fvsPatchField<scalar>& magSfPatch = mesh.magSf().boundaryField()[patchi];
	const fvPatchField<scalar>& pPatch = p.boundaryField()[patchi];
	forAll(mesh.boundary()[patchi], facei) {
	    if (mesh.boundary()[patchi].name() == "wing") {
		// col 1 for h, zero vector
		// col 2 for alpha
		tList.push_back(Triplet(nBdIndex*5+1, 1, SfPatch[facei].y()*p[facei]));
//		tList.push_back(Triplet(nBdIndex*5+2, 1, -SfPatch[facei].x()*p[facei]));
//		tList.push_back(Triplet(nBdIndex*5+4, 1, -SfPatch[facei].y()*p[facei]));
		// col3 for h_dot
		tList.push_back(Triplet(nBdIndex*5+4, 2, -SfPatch[facei].y()*p[facei]));
		// col4 for alpha_dot
		scalar Lx = CfPatch[facei].x() - x0;
		scalar Ly = CfPatch[facei].y();
		tList.push_back(Triplet(nBdIndex*5+4, 3, -(SfPatch[facei].x()*Ly-SfPatch[facei].y()*Lx)*p[facei]));
//		tList.push_back(Triplet(nBdIndex*5+4, 3, -(-SfPatch[facei].y()*Lx)*p[facei]));

	    }
	    ++nBdIndex;
	}
    }
    
    matBf.setFromTriplets(tList.begin(), tList.end());
}

void getMatBdyn1
(
    SpMat& matBdyn1,
    const surfaceScalarField& rho_pos,
    const surfaceVectorField& U_pos,
    const surfaceScalarField& H_pos,
    const surfaceScalarField& p_pos,
    const surfaceScalarField& rho_neg,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& H_neg,
    const surfaceScalarField& p_neg
)
{
    scalar x0 = .25;
    const fvMesh& mesh = rho_pos.mesh();
    const labelUList& P = mesh.owner();
    TripletList tList;
    tList.reserve(P.size()*5);
   
    HashSet<int> surfacePoints;
    forAll(mesh.boundary(), patchi) {
	if (mesh.boundary()[patchi].name() == "wing") {
	    for (int n = 0; n < mesh.boundary()[patchi].size(); n++) {
		int fn = mesh.boundary()[patchi].start() + n;
		forAll (mesh.faces()[fn], p) {
		    surfacePoints.insert(mesh.faces()[fn][p]);
		}
	    }
	}
    }    

    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceVectorField& Cf = mesh.Cf();
    const pointField& pf = mesh.points();
    for (int fi=0; fi<P.size(); fi++) {
	forAll (mesh.faces()[fi], p) {
	    if (surfacePoints[mesh.faces()[fi][p]]) {
		scalar Lx = pf[mesh.faces()[fi][p]].x() - x0;
		scalar dSx = sign(Cf[fi].y())*sign(Sf[fi].x())*Lx*.15;
/*		Info << Cf[fi].x() << endl;
		Info << pf[mesh.faces()[fi][p]].x() << endl;
		Info << Lx << endl;
		Info << Sf[fi].x() << endl;
		Info << dSx << endl;
		cin.get();*/
		scalar phi_pos = rho_pos[fi]*U_pos[fi].x();
		scalar phi_neg = rho_neg[fi]*U_neg[fi].x();
		tList.push_back(Triplet(fi*5+0, 1, .5*(phi_pos+phi_neg)*dSx));
		tList.push_back(Triplet(fi*5+1, 1, .5*(phi_pos*U_pos[fi].x()+phi_neg*U_neg[fi].x()+p_pos[fi]+p_neg[fi])*dSx));
		tList.push_back(Triplet(fi*5+2, 1, .5*(phi_pos*U_pos[fi].y()+phi_neg*U_neg[fi].y())*dSx));
		tList.push_back(Triplet(fi*5+4, 1, .5*(phi_pos*H_pos[fi]+phi_neg*H_neg[fi])*dSx));
		break;
	    }
	}
    }
    
    matBdyn1.setFromTriplets(tList.begin(), tList.end());
}

void getMatFSp
(
    SpMat& matF,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p,
    const volScalarField& c,
    const unsigned int nBdFaces,
    const fvMesh& mesh
)
{
    // relationship of rho, U, p on the boundary.
    TripletList tList;
//    tList.reserve(nBdFaces*5);
    tList.reserve(nBdFaces*6);
    
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	//Us: + for outflow, - for inflow
	const scalarField Us = U.boundaryField()[patchi]&mesh.Sf().boundaryField()[patchi];
	const scalarField Cs = c.boundaryField()[patchi]*mesh.magSf().boundaryField()[patchi];
	//flowRegime + for supersonic, - for subsonic
	const scalarField flowRegime = mag(Us) - Cs;

	const fvsPatchField<vector>& SfPatch = mesh.Sf().boundaryField()[patchi];
	const fvsPatchField<scalar>& magSfPatch = mesh.magSf().boundaryField()[patchi];
	const scalarField pPif = p.boundaryField()[patchi].patchInternalField();
	const scalarField rhoPif = rho.boundaryField()[patchi].patchInternalField();
	const vectorField UPif = U.boundaryField()[patchi].patchInternalField();
	const scalarField cPif = c.boundaryField()[patchi].patchInternalField();
	const labelUList& faceCells =
	mesh.boundary()[patchi].faceCells();

	forAll(mesh.boundary()[patchi], facei) {
	    if ((mesh.boundary()[patchi].name() == "inlet") || (mesh.boundary()[patchi].name() == "outlet")) {
		if(Us[facei] > 0) { // outflow
		    if(flowRegime[facei] > 0) { //supersonic
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5, 1));
//			rhoPatch[facei] = rhoPif[facei];
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, 1));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, 1));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, 1));
//			UPatch[facei] = UPif[facei];
//			UPatch[facei].z() = 0;
			tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+4, 1));
//			pPatch[facei] = pPif[facei];
		    } else { //subsonic
			vector nf = SfPatch[facei]/magSfPatch[facei];
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5, 1));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+4, -1/cPif[facei]/cPif[facei]));
//			rhoPatch[facei] = rhoPif[facei] + (pPatch[facei]-pPif[facei])/cPif[facei]/cPif[facei];
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, 1));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+4, nf.x()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, 1));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+4, nf.y()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, 1));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+4, nf.z()/rhoPif[facei]/cPif[facei]));
//			UPatch[facei] = UPif[facei]-nf*(pPatch[facei]-pPif[facei])/rhoPif[facei]/cPif[facei];
//			UPatch[facei].z() = 0;
//			pPatch[facei] = p_inf;
		    }
		} else { // inflow
		    if(flowRegime[facei] > 0) { //supersonic
/*			pPatch[facei] = p_inf;
			rhoPatch[facei] = rho_inf;
			UPatch[facei] = U_inf;
			UPatch[facei].z() = 0; */
		    } else { //subsonic
			vector nf = SfPatch[facei]/magSfPatch[facei];
			tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+4, .5));
			tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+1, .5*rhoPif[facei]*cPif[facei]*nf.x()));
			tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+2, .5*rhoPif[facei]*cPif[facei]*nf.y()));
			tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+3, .5*rhoPif[facei]*cPif[facei]*nf.z()));
//			pPatch[facei] = .5*(p_inf+pPif[facei]-rhoPif[facei]*cPif[facei]*((U_inf-UPif[facei])&nf));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+4, .5/cPif[facei]/cPif[facei]));
/*			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+1, .5*rhoPif[facei]*cPif[facei]*nf.x()/cPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+2, .5*rhoPif[facei]*cPif[facei]*nf.y()/cPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+3, .5*rhoPif[facei]*cPif[facei]*nf.z()/cPif[facei]/cPif[facei]));*/
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+1, .5*rhoPif[facei]*nf.x()/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+2, .5*rhoPif[facei]*nf.y()/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+3, .5*rhoPif[facei]*nf.z()/cPif[facei]));
//			rhoPatch[facei] = rho_inf + (pPatch[facei]-p_inf)/cPif[facei]/cPif[facei];
/*			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+4, .5*nf.x()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, .5*rhoPif[facei]*cPif[facei]*nf.x()*nf.x()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+2, .5*rhoPif[facei]*cPif[facei]*nf.y()*nf.x()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+3, .5*rhoPif[facei]*cPif[facei]*nf.z()*nf.x()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+4, .5*nf.y()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+1, .5*rhoPif[facei]*cPif[facei]*nf.x()*nf.y()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, .5*rhoPif[facei]*cPif[facei]*nf.y()*nf.y()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+3, .5*rhoPif[facei]*cPif[facei]*nf.z()*nf.y()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+4, .5*nf.z()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+1, .5*rhoPif[facei]*cPif[facei]*nf.x()*nf.z()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+2, .5*rhoPif[facei]*cPif[facei]*nf.y()*nf.z()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, .5*rhoPif[facei]*cPif[facei]*nf.z()*nf.z()/rhoPif[facei]/cPif[facei]));*/
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+4, .5*nf.x()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, .5*nf.x()*nf.x()));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+2, .5*nf.y()*nf.x()));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+3, .5*nf.z()*nf.x()));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+4, .5*nf.y()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+1, .5*nf.x()*nf.y()));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, .5*nf.y()*nf.y()));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+3, .5*nf.z()*nf.y()));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+4, .5*nf.z()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+1, .5*nf.x()*nf.z()));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+2, .5*nf.y()*nf.z()));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, .5*nf.z()*nf.z()));
//			UPatch[facei] = U_inf-nf*(p_inf-pPatch[facei])/rhoPif[facei]/cPif[facei];
//			UPatch[facei].z() = 0;
		    }
		}
	    }

	    if (mesh.boundary()[patchi].name() == "wing") {
		vector nf = SfPatch[facei]/magSfPatch[facei];
		tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+4, 1));
//		pPatch[facei] = pPif[facei];
		tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5, 1));
//		rhoPatch[facei] = rhoPif[facei];
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, 1-nf.x()*nf.x()));
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+2, -nf.y()*nf.x()));
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+3, -nf.z()*nf.x()));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+1, -nf.x()*nf.y()));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, 1-nf.y()*nf.y()));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+3, -nf.z()*nf.y()));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+1, -nf.x()*nf.z()));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+2, -nf.y()*nf.z()));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, 1-nf.z()*nf.z()));
//		UPatch[facei] = UPif[facei] - (UPif[facei]&nf)*nf;
//		UPatch[facei].z() = 0;
	    }
	    ++nBdIndex;
	}
    }
    matF.setFromTriplets(tList.begin(), tList.end());
}

void getMatHSp
(
    SpMat& matH,
    const volScalarField& rho,
    const volVectorField& U,
    const scalar U_freestream,
    const unsigned int nBdFaces,
    const fvMesh& mesh
)
{
    TripletList tList;
    tList.reserve(nBdFaces*5);
    
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	forAll(mesh.boundary()[patchi], facei) {
	    if (mesh.boundary()[patchi].name() == "wing") {
    //	    tList.push_back(Triplet(nBdIndex*5, 0, 0));
		tList.push_back
		(
		    Triplet
		    (
			nBdIndex*5+1,
			0,
			-std::sin(0.0)*U_freestream*rho.boundaryField()[patchi][facei]
		    )
		);
		tList.push_back
		(
		    Triplet
		    (
			nBdIndex*5+2,
			0,
			std::cos(0.0)*U_freestream*rho.boundaryField()[patchi][facei]
		    )
		);
    //	    tList.push_back(Triplet(nBdIndex*5+3, 0, 0));
		tList.push_back
		(
		    Triplet
		    (
			nBdIndex*5+4,
			0,
			rho.boundaryField()[patchi][facei]*
			(
			    -std::sin(0.0)*U_freestream*U.boundaryField()[patchi][facei].x()+
			    std::cos(0.0)*U_freestream*U.boundaryField()[patchi][facei].y()
			)
		    )
		);
	    }
	    ++nBdIndex;
	}
    }
    
    matH.setFromTriplets(tList.begin(), tList.end());
}

void getMatKSp
(
    SpMat& matK,
    const fvMesh& mesh
)
{
    scalar x0 = 0.25;
    forAll(mesh.boundary(), patchi) {
	const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
	const fvsPatchField<vector>& SfPatch = mesh.Sf().boundaryField()[patchi];
	const fvsPatchField<vector>& CfPatch = mesh.Cf().boundaryField()[patchi];
	if (mesh.boundary()[patchi].name() == "wing") {
	    TripletList tList;
	    tList.reserve(mesh.boundary()[patchi].size()*3);
	    forAll(mesh.boundary()[patchi], facei) {
		scalar Lx = CfPatch[facei].x() - x0;
		scalar Ly = CfPatch[facei].y();
		scalar Sx = SfPatch[facei].x();
		scalar Sy = SfPatch[facei].y();
		tList.push_back(Triplet(0, faceCells[facei]*5+4, Sy));
		tList.push_back(Triplet(1, faceCells[facei]*5+4, (Ly*Sx-Lx*Sy)));
	    }
	    matK.setFromTriplets(tList.begin(), tList.end());
	}
    }
}

void getMatKpSp
(
    SpMat& matKp,
    const fvMesh& mesh
)
{
    forAll(mesh.boundary(), patchi) {
	const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
	if (mesh.boundary()[patchi].name() == "wing") {
	    TripletList tList;
	    tList.reserve(mesh.boundary()[patchi].size()*3);
	    forAll(mesh.boundary()[patchi], facei) {
		tList.push_back(Triplet(facei, faceCells[facei]*5+4, 1));
	    }
	    matKp.setFromTriplets(tList.begin(), tList.end());
	}
    }
}

void getMatDivBdSp
(
    SpMat& matDivBd,
    const unsigned int nBdFaces,
    const fvMesh& mesh
)
{
    TripletList tList;
    tList.reserve(nBdFaces*5);
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
	forAll(mesh.boundary()[patchi], facei) {
	    for(unsigned int vn = 0; vn < 5; vn++) {
		tList.push_back
		(
		    Triplet(faceCells[facei]*5+vn,
		    nBdIndex*5+vn,
		    1/mesh.V()[faceCells[facei]])
		);
	    }
	    nBdIndex++;
	}
    }
    matDivBd.setFromTriplets(tList.begin(), tList.end());
}

void getMatDivSp
(
    SpMat& matDiv,
    const fvMesh& mesh
)
{
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    TripletList tList;
    tList.reserve(P.size()*5);

    for (int fi=0; fi<P.size(); fi++) {
	for(unsigned int vn = 0; vn < 5; vn++) {
	    tList.push_back(Triplet(P[fi]*5+vn, fi*5+vn, 1/mesh.V()[P[fi]]));
	    tList.push_back(Triplet(N[fi]*5+vn, fi*5+vn, -1/mesh.V()[N[fi]]));
	}
    }
    
    matDiv.setFromTriplets(tList.begin(), tList.end());
}

void getMatESp
(
    SpMat& matE,
    const volScalarField& rho,
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const unsigned int nBdFaces,
    const fvMesh& mesh,
    double U_freestream
)
{
    TripletList tList;
    tList.reserve(nBdFaces*3);
    
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
	if 
	(
	    p.boundaryField()[patchi].coupled() ||
	    U.boundaryField()[patchi].coupled() ||
	    T.boundaryField()[patchi].coupled() 
	) {
	    Info << "Coupled Patch exists" << endl;
	    cin.get();
	}

	forAll(mesh.boundary()[patchi], facei) {
	    word UType = U.boundaryField().types()[patchi];
	    const fvPatchField<scalar>& rhob = rho.boundaryField()[patchi];
	    if (mesh.boundary()[patchi].name() == "inlet") {
		if (UType == "zeroGradient") {
		} else if (UType == "fixedValue") {
		    tList.push_back(Triplet(nBdIndex*5+1, 0, -U_freestream*std::sin(0/57.296)*rhob[facei]));
		    tList.push_back(Triplet(nBdIndex*5+2, 0, U_freestream*std::cos(0/57.296)*rhob[facei]));
		    tList.push_back(Triplet(nBdIndex*5+3, 0, 0));
		}
	    }
	    ++nBdIndex;
	}
    }
    
    matE.setFromTriplets(tList.begin(), tList.end());
}

void getMatQsSp
(
    SpMat& matQic,
    SpMat& matQif_pos,
    SpMat& matQif_neg,
    SpMat& matQb,
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    const surfaceScalarField& rho_pos,
    const surfaceScalarField& rho_neg,
    const surfaceVectorField& rhoU_pos,
    const surfaceVectorField& rhoU_neg,
    const surfaceScalarField& rhoE_pos,
    const surfaceScalarField& rhoE_neg
)
{
    const fvMesh& mesh = rho.mesh();
    const labelUList& owner = mesh.owner();
    TripletList listQif_pos, listQif_neg;
    listQif_pos.reserve(owner.size()*5);
    listQif_neg.reserve(owner.size()*5);

    int nBdFaces = 0;
    forAll(mesh.boundary(), patchi) {
	forAll(mesh.boundary()[patchi], fi) {
	    nBdFaces++;
	}
    }
    TripletList listQb;
    listQif_pos.reserve(nBdFaces*5);

    const Field<scalar>& rhoi_pos = rho_pos.internalField();
    const Field<vector>& rhoUi_pos = rhoU_pos.internalField();
    const Field<scalar>& rhoEi_pos = rhoE_pos.internalField();
    const Field<scalar>& rhoi_neg = rho_neg.internalField();
    const Field<vector>& rhoUi_neg = rhoU_neg.internalField();
    const Field<scalar>& rhoEi_neg = rhoE_neg.internalField();

    forAll(owner, fi)
    {
	listQif_pos.push_back(Triplet(fi*5+0, 0, rho_pos[fi]));
	listQif_pos.push_back(Triplet(fi*5+1, 0, rhoU_pos[fi].x()));
	listQif_pos.push_back(Triplet(fi*5+2, 0, rhoU_pos[fi].y()));
	listQif_pos.push_back(Triplet(fi*5+3, 0, rhoU_pos[fi].z()));
	listQif_pos.push_back(Triplet(fi*5+4, 0, rhoE_pos[fi]));

	listQif_neg.push_back(Triplet(fi*5+0, 0, rho_neg[fi]));
	listQif_neg.push_back(Triplet(fi*5+1, 0, rhoU_neg[fi].x()));
	listQif_neg.push_back(Triplet(fi*5+2, 0, rhoU_neg[fi].y()));
	listQif_neg.push_back(Triplet(fi*5+3, 0, rhoU_neg[fi].z()));
	listQif_neg.push_back(Triplet(fi*5+4, 0, rhoE_neg[fi]));
    }
    
    int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	const fvPatchField<scalar>& rhob = rho.boundaryField()[patchi];
	const fvPatchField<vector>& rhoUb = rhoU.boundaryField()[patchi];
	const fvPatchField<scalar>& rhoEb = rhoE.boundaryField()[patchi];
	forAll(mesh.boundary()[patchi], fi) {
	    listQb.push_back(Triplet(nBdIndex*5+0, 0, rhob[fi]));
	    listQb.push_back(Triplet(nBdIndex*5+1, 0, rhoUb[fi].x()));
	    listQb.push_back(Triplet(nBdIndex*5+2, 0, rhoUb[fi].y()));
	    listQb.push_back(Triplet(nBdIndex*5+3, 0, rhoUb[fi].z()));
	    listQb.push_back(Triplet(nBdIndex*5+4, 0, rhoEb[fi]));
	    nBdIndex++;
	}
    }

    matQb.setFromTriplets(listQb.begin(), listQb.end());
    matQif_pos.setFromTriplets(listQif_pos.begin(), listQif_pos.end());
    matQif_neg.setFromTriplets(listQif_neg.begin(), listQif_neg.end());

    /*for(unsigned int ni = 0; ni < nIntCells; ni++) {
    matQic.push_back(std::vector<double>(1,rho.internalField()[ni]));
    matQic.push_back(std::vector<double>(1,rhoU.internalField()[ni].x()));
    matQic.push_back(std::vector<double>(1,rhoU.internalField()[ni].y()));
    matQic.push_back(std::vector<double>(1,rhoU.internalField()[ni].z()));
    matQic.push_back(std::vector<double>(1,rhoE.internalField()[ni]));
}*/
}

void getMatTSp
(
    SpMat& matLimiter_pos,
    SpMat& matLimiter_neg,
    Field<scalar>& limiterRho_pos,
    Field<scalar>& limiterRho_neg,
    Field<vector>& limiterRhoU_pos,
    Field<vector>& limiterRhoU_neg,
    Field<scalar>& limiterRhoE_pos,
    Field<scalar>& limiterRhoE_neg
)
{
    TripletList listLmt_pos, listLmt_neg;
    listLmt_pos.reserve(limiterRho_pos.size()*5);
    listLmt_neg.reserve(limiterRho_pos.size()*5);
    
    for(int i = 0; i < limiterRho_pos.size(); ++i) {
	listLmt_pos.push_back(Triplet(i*5+0, 0, limiterRho_pos[i]));
	listLmt_pos.push_back(Triplet(i*5+1, 0, limiterRhoU_pos[i].x()));
	listLmt_pos.push_back(Triplet(i*5+2, 0, limiterRhoU_pos[i].y()));
	listLmt_pos.push_back(Triplet(i*5+3, 0, limiterRhoU_pos[i].z()));
	listLmt_pos.push_back(Triplet(i*5+4, 0, limiterRhoE_pos[i]));
	
	listLmt_neg.push_back(Triplet(i*5+0, 0, limiterRho_neg[i]));
	listLmt_neg.push_back(Triplet(i*5+1, 0, limiterRhoU_neg[i].x()));
	listLmt_neg.push_back(Triplet(i*5+2, 0, limiterRhoU_neg[i].y()));
	listLmt_neg.push_back(Triplet(i*5+3, 0, limiterRhoU_neg[i].z()));
	listLmt_neg.push_back(Triplet(i*5+4, 0, limiterRhoE_neg[i]));
    }

    matLimiter_pos.setFromTriplets(listLmt_pos.begin(), listLmt_pos.end());
    matLimiter_neg.setFromTriplets(listLmt_neg.begin(), listLmt_neg.end());
}

////////////////////////new Order///////////////////////////////
void newOrdering
(
    SpMat& matOdr,
    int nqIntCells
)
{
    TripletList listOdr;
    listOdr.reserve(nqIntCells);
    
    for (int i=0; i < nqIntCells; ++i) { 
	int nIntCells = (nqIntCells/5);
	int a = (i%nIntCells)*5+(i/nIntCells);
	int b = i;
	listOdr.push_back(Triplet(a, b, 1));
    }
    matOdr.setFromTriplets(listOdr.begin(), listOdr.end());
}

void getMatSSp
(
    SpMat& matS,
    const fvMesh& mesh
)
{
    TripletList tList;
    
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	const fvsPatchField<vector>& psf = mesh.Sf().boundaryField()[patchi];
	forAll(mesh.boundary()[patchi], facei) {
	    tList.push_back(Triplet(nBdIndex*5+1, 0, psf[facei].x()));
	    tList.push_back(Triplet(nBdIndex*5+2, 0, psf[facei].y()));
	    tList.push_back(Triplet(nBdIndex*5+3, 0, psf[facei].z()));
	    ++nBdIndex;
	}
    }
    
    matS.setFromTriplets(tList.begin(), tList.end());
}
/*
	listOdr.push_back(Triplet(i*5+0, i, 1));
	listOdr.push_back(Triplet(i*5+1, nIntCells+i, 1));
	listOdr.push_back(Triplet(i*5+2, 2*nIntCells+i, 1));
	listOdr.push_back(Triplet(i*5+3, 3*nIntCells+i, 1));
	listOdr.push_back(Triplet(i*5+4, 4*nIntCells+i, 1));

	if(i != 0) {
	    listOdr.push_back(Triplet(i, i*5+0, 1));
	}
	listOdr.push_back(Triplet(nIntCells+i, i*5+1, 1));
	listOdr.push_back(Triplet(2*nIntCells+i, i*5+2, 1));
	listOdr.push_back(Triplet(3*nIntCells+i, i*5+3, 1));
	if(i != nIntCells) {
	    listOdr.push_back(Triplet(4*nIntCells+i, i*5+4, 1));
	}
 * */

void getMatCfSp
(
    SpMat& matCf,
    const fvMesh& mesh
)
{
    forAll(mesh.boundary(), patchi) {
	const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
	if (mesh.boundary()[patchi].name() == "wing") {
	    TripletList tList;
	    tList.reserve(mesh.boundary()[patchi].size()*3);
	    forAll(mesh.boundary()[patchi], facei) {
//		tList.push_back(Triplet(0, faceCells[facei]*5, mesh.Cf().boundaryField()[patchi][facei].x()));
		tList.push_back(Triplet(0, facei, mesh.Cf().boundaryField()[patchi][facei].x()));
	    }
	    matCf.setFromTriplets(tList.begin(), tList.end());
	}
    }
}

void getMatSfSp
(
    SpMat& matSf,
    const fvMesh& mesh
)
{
    forAll(mesh.boundary(), patchi) {
	const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
	if (mesh.boundary()[patchi].name() == "wing") {
	    TripletList tList;
	    tList.reserve(mesh.boundary()[patchi].size()*3);
	    forAll(mesh.boundary()[patchi], facei) {
		tList.push_back(Triplet(0, facei, mesh.magSf().boundaryField()[patchi][facei]));
	    }
	    matSf.setFromTriplets(tList.begin(), tList.end());
	}
    }
}


/*
myMat con2phy
(
    const scalar& rho,
    const vector& U,
    const scalar& T,
    const scalar& p,
    const scalar& gamma
)
{
	myMat subD
	(
	    5,
	    std::vector<double>(5,0)
	);

	double u = U.x();
	double v = U.y();
	double w = U.z();
	double q = u*u + v*v + w*w;
	//row1
	subD[0][0] = .5*(gamma-1)*q;
	subD[0][1] = -(gamma-1)*u;
	subD[0][2] = -(gamma-1)*v;
	subD[0][3] = -(gamma-1)*w;
	subD[0][4] = (gamma-1);
	//row2
	subD[1][0] = -1.*u/rho;
	subD[1][1] = 1./rho;
	subD[1][2] = 0;
	subD[1][3] = 0;
	subD[1][4] = 0;
	//row3
	subD[2][0] = -1.*v/rho;
	subD[2][1] = 0;
	subD[2][2] = 1./rho;
	subD[2][3] = 0;
	subD[2][4] = 0;
	//row4
	subD[3][0] = -1.*w/rho;
	subD[3][1] = 0;
	subD[3][2] = 0;
	subD[3][3] = 1./rho;
	subD[3][4] = 0;
	//row5
	subD[4][0] = .5*(gamma-1)*T*q/p - T/rho;
	subD[4][1] = -(gamma-1)*u*T/p;
	subD[4][2] = -(gamma-1)*v*T/p;
	subD[4][3] = -(gamma-1)*w*T/p;
	subD[4][4] = (gamma-1)*T/p;

	return subD;
}

myMat phy2con
(
    const scalar& rho,
    const vector& U,
    const scalar& T,
    const scalar& p,
    const scalar& gamma
)
{
	myMat subD
	(
	    5,
	    std::vector<double>(5,0)
	);

	double u = U.x();
	double v = U.y();
	double w = U.z();
	double q = u*u + v*v + w*w;
	//row1
	subD[0][0] = rho/p;
	subD[0][1] = 0;
	subD[0][2] = 0;
	subD[0][3] = 0;
	subD[0][4] = -rho/T;
	//row2
	subD[1][0] = rho/p*u;
	subD[1][1] = rho;
	subD[1][2] = 0;
	subD[1][3] = 0;
	subD[1][4] = -rho*u/T;
	//row3
	subD[2][0] = rho/p*v;
	subD[2][1] = 0;
	subD[2][2] = rho;
	subD[2][3] = 0;
	subD[2][4] = -rho*v/T;
	//row4
	subD[3][0] = rho/p*w;
	subD[3][1] = 0;
	subD[3][2] = 0;
	subD[3][3] = rho;
	subD[3][4] = -rho*w/T;
	//row5
	subD[4][0] = 1/(gamma-1) + .5*rho*q/p;
	subD[4][1] = rho*u;
	subD[4][2] = rho*v;
	subD[4][3] = rho*w;
	subD[4][4] = -.5*rho*q/T;

	return subD;
}

void getMatFSp
(
    SpMat& matF,
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const unsigned int nBdFaces,
    const fvMesh& mesh
)
{
    TripletList tList;
//    tList.reserve(nBdFaces*5);
    tList.reserve(nBdFaces*6);
    
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
	const fvsPatchField<vector>& psf = mesh.Sf().boundaryField()[patchi];
	if 
	(
	    p.boundaryField()[patchi].coupled() ||
	    U.boundaryField()[patchi].coupled() ||
	    T.boundaryField()[patchi].coupled() 
	) {
	    Info << "Coupled Patch exists" << endl;
	    cin.get();
	}

	forAll(mesh.boundary()[patchi], facei) {
	    word pType = p.boundaryField().types()[patchi];
	    if (pType == "zeroGradient") {
		tList.push_back(Triplet(nBdIndex*5+0, faceCells[facei]*5, 1));
	    } else if (pType == "fixedValue") {
		tList.push_back(Triplet(nBdIndex*5+0, faceCells[facei]*5, 0));
	    }
	    word UType = U.boundaryField().types()[patchi];
	    if (UType == "zeroGradient") {
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, 1));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, 1));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, 1));
	    } else if (UType == "fixedValue") {
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, 0));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, 0));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, 0));
	    } else if (UType == "slip") {
		vector unitSf = psf[facei]/mag(psf[facei]);
		typename outerProduct<vector, vector>::type ss(unitSf*unitSf);
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, 1));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, 1));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, 1));
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, -ss.xx()));
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+2, -ss.xy()));
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+3, -ss.xz()));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+1, -ss.yx()));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, -ss.yy()));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+3, -ss.yz()));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+1, -ss.zx()));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+2, -ss.zy()));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, ss.zz()));
	    }
	    word TType = T.boundaryField().types()[patchi];
	    if (TType == "zeroGradient") {
		tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+4, 1));
	    } else if (TType == "fixedValue") {
		tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+4, 0));
	    }

	    ++nBdIndex;
	}
    }
    
    matF.setFromTriplets(tList.begin(), tList.end());
}

void getMatFSp
(
    SpMat& matF,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p,
    const volScalarField& c,
    const unsigned int nBdFaces,
    const fvMesh& mesh
)
{
    // relationship of rho, U, p on the boundary.
    TripletList tList;
//    tList.reserve(nBdFaces*5);
    tList.reserve(nBdFaces*6);
    
    unsigned int nBdIndex = 0;
    forAll(mesh.boundary(), patchi) {
	//Us: + for outflow, - for inflow
	const scalarField Us = U.boundaryField()[patchi]&mesh.Sf().boundaryField()[patchi];
	const scalarField Cs = c.boundaryField()[patchi]*mesh.magSf().boundaryField()[patchi];
	//flowRegime + for supersonic, - for subsonic
	const scalarField flowRegime = mag(Us) - Cs;

//	fvPatchField<scalar>& pPatch = p.boundaryField()[patchi];
//	fvPatchField<scalar>& rhoPatch = rho.boundaryField()[patchi];
//	fvPatchField<vector>& UPatch = U.boundaryField()[patchi];
	const fvsPatchField<vector>& SfPatch = mesh.Sf().boundaryField()[patchi];
	const fvsPatchField<scalar>& magSfPatch = mesh.magSf().boundaryField()[patchi];
	const scalarField pPif = p.boundaryField()[patchi].patchInternalField();
	const scalarField rhoPif = rho.boundaryField()[patchi].patchInternalField();
	const vectorField UPif = U.boundaryField()[patchi].patchInternalField();
	const scalarField cPif = c.boundaryField()[patchi].patchInternalField();
	const labelUList& faceCells =
	mesh.boundary()[patchi].faceCells();

	forAll(mesh.boundary()[patchi], facei) {
	    if ((mesh.boundary()[patchi].name() == "inlet") || (mesh.boundary()[patchi].name() == "outlet")) {
		if(Us[facei] > 0) { // outflow
		    if(flowRegime[facei] > 0) { //supersonic
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5, 1));
//			rhoPatch[facei] = rhoPif[facei];
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, 1));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, 1));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, 1));
//			UPatch[facei] = UPif[facei];
//			UPatch[facei].z() = 0;
			tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+4, 1));
//			pPatch[facei] = pPif[facei];
		    } else { //subsonic
			vector nf = SfPatch[facei]/magSfPatch[facei];
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5, 1));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+4, -1/cPif[facei]/cPif[facei]));
//			rhoPatch[facei] = rhoPif[facei] + (pPatch[facei]-pPif[facei])/cPif[facei]/cPif[facei];
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, 1));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+4, nf.x()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, 1));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+4, nf.y()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, 1));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+4, nf.z()/rhoPif[facei]/cPif[facei]));
//			UPatch[facei] = UPif[facei]-nf*(pPatch[facei]-pPif[facei])/rhoPif[facei]/cPif[facei];
//			UPatch[facei].z() = 0;
//			pPatch[facei] = p_inf;
		    }
		} else { // inflow
		    if(flowRegime[facei] > 0) { //supersonic
//			pPatch[facei] = p_inf;
//			rhoPatch[facei] = rho_inf;
//			UPatch[facei] = U_inf;
//			UPatch[facei].z() = 0; 
		    } else { //subsonic
			vector nf = SfPatch[facei]/magSfPatch[facei];
			tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+4, .5));
			tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+1, .5*rhoPif[facei]*cPif[facei]*nf.x()));
			tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+2, .5*rhoPif[facei]*cPif[facei]*nf.y()));
			tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+3, .5*rhoPif[facei]*cPif[facei]*nf.z()));
//			pPatch[facei] = .5*(p_inf+pPif[facei]-rhoPif[facei]*cPif[facei]*((U_inf-UPif[facei])&nf));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+4, .5/cPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+1, .5*rhoPif[facei]*cPif[facei]*nf.x()/cPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+2, .5*rhoPif[facei]*cPif[facei]*nf.y()/cPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+3, .5*rhoPif[facei]*cPif[facei]*nf.z()/cPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+1, .5*rhoPif[facei]*nf.x()/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+2, .5*rhoPif[facei]*nf.y()/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5+3, .5*rhoPif[facei]*nf.z()/cPif[facei]));
//			rhoPatch[facei] = rho_inf + (pPatch[facei]-p_inf)/cPif[facei]/cPif[facei];
//			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+4, .5*nf.x()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, .5*rhoPif[facei]*cPif[facei]*nf.x()*nf.x()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+2, .5*rhoPif[facei]*cPif[facei]*nf.y()*nf.x()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+3, .5*rhoPif[facei]*cPif[facei]*nf.z()*nf.x()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+4, .5*nf.y()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+1, .5*rhoPif[facei]*cPif[facei]*nf.x()*nf.y()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, .5*rhoPif[facei]*cPif[facei]*nf.y()*nf.y()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+3, .5*rhoPif[facei]*cPif[facei]*nf.z()*nf.y()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+4, .5*nf.z()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+1, .5*rhoPif[facei]*cPif[facei]*nf.x()*nf.z()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+2, .5*rhoPif[facei]*cPif[facei]*nf.y()*nf.z()/rhoPif[facei]/cPif[facei]));
//			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, .5*rhoPif[facei]*cPif[facei]*nf.z()*nf.z()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+4, .5*nf.x()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, .5*nf.x()*nf.x()));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+2, .5*nf.y()*nf.x()));
			tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+3, .5*nf.z()*nf.x()));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+4, .5*nf.y()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+1, .5*nf.x()*nf.y()));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, .5*nf.y()*nf.y()));
			tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+3, .5*nf.z()*nf.y()));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+4, .5*nf.z()/rhoPif[facei]/cPif[facei]));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+1, .5*nf.x()*nf.z()));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+2, .5*nf.y()*nf.z()));
			tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, .5*nf.z()*nf.z()));
//			UPatch[facei] = U_inf-nf*(p_inf-pPatch[facei])/rhoPif[facei]/cPif[facei];
//			UPatch[facei].z() = 0;
		    }
		}
	    }

	    if (mesh.boundary()[patchi].name() == "wing") {
		vector nf = SfPatch[facei]/magSfPatch[facei];
		tList.push_back(Triplet(nBdIndex*5+4, faceCells[facei]*5+4, 1));
//		pPatch[facei] = pPif[facei];
		tList.push_back(Triplet(nBdIndex*5, faceCells[facei]*5, 1));
//		rhoPatch[facei] = rhoPif[facei];
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+1, 1-nf.x()*nf.x()));
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+2, -nf.y()*nf.x()));
		tList.push_back(Triplet(nBdIndex*5+1, faceCells[facei]*5+3, -nf.z()*nf.x()));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+1, -nf.x()*nf.y()));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+2, 1-nf.y()*nf.y()));
		tList.push_back(Triplet(nBdIndex*5+2, faceCells[facei]*5+3, -nf.z()*nf.y()));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+1, -nf.x()*nf.z()));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+2, -nf.y()*nf.z()));
		tList.push_back(Triplet(nBdIndex*5+3, faceCells[facei]*5+3, 1-nf.z()*nf.z()));
//		UPatch[facei] = UPif[facei] - (UPif[facei]&nf)*nf;
//		UPatch[facei].z() = 0;
	    }
	    ++nBdIndex;
	}
    }
    matF.setFromTriplets(tList.begin(), tList.end());
}
 * */
