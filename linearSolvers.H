#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

using Eigen::VectorXd;
typedef Eigen::Map<VectorXd> MapVec;

void myVecWriteSp(const MapVec& mat, const char* fileName)
{
    std::ofstream out(fileName);
    out.precision(17);
    out << mat << std::endl;
}

void myVecWriteSp(const VectorXd& mat, const char* fileName)
{
    std::ofstream out(fileName);
    out.precision(17);
    out << mat << std::endl;
}

int myDdt
(
    SpMat& A,
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const psiThermo& thermo,
    const surfaceVectorField& U_pos,
    const surfaceScalarField& H_pos,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& H_neg,
    const surfaceVectorField& U_tilt,
    const surfaceScalarField& H_tilt,
    const surfaceScalarField& gamma,
    double cfl,
    const word scheme = "EulerExplicit"
);

int mySolve
(
    volScalarField& rho,
    volVectorField& rhoU,
    volScalarField& rhoE,
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const psiThermo& thermo,
    const volScalarField& phi,
    const volVectorField& phiUp,
    const volScalarField& phiEp,
    const surfaceVectorField& U_pos,
    const surfaceScalarField& H_pos,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& H_neg,
    const surfaceVectorField& U_tilt,
    const surfaceScalarField& H_tilt,
    const surfaceScalarField& gamma,
    double cfl,
    const word scheme = "EulerExplicit"
)
{
    //import to Eigen; try to speed up this part
    label n = rho.internalField().size();
    int nq = 5*n;

    scalarField flux(phi);
    flux.append(phiUp.component(0));
    flux.append(phiUp.component(1));
    flux.append(phiUp.component(2));
    flux.append(phiEp);

    MapVec b0(flux.data(), nq);
    SpMat matOdr(nq, nq);
    newOrdering
    (
	matOdr,
	nq
    );
    
    VectorXd b(nq);
    b = matOdr*b0;

    VectorXd x(nq);
    SpMat A(nq, nq);

    myDdt
    (
	A,
       	rho,
       	rhoU,
       	rhoE,
	p,
	U,
	T,
       	thermo,
	U_pos,
	H_pos,
	U_neg,
	H_neg,
	U_tilt,
	H_tilt,
	gamma,
	cfl,
	scheme
//       	"EulerExplicit"
    ); //coeff update to speed up this part
    
    TripletList tList, eyeList;
    tList.reserve(nq);
    eyeList.reserve(nq);

    for(int i = 0; i < nq; ++i) {
	eyeList.push_back(Triplet(i, i, 1));
    }
    const fvMesh& mesh = rho.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    scalarField contraVSumX(mesh.nCells(), 0);
    scalarField contraVSumY(mesh.nCells(), 0);
    scalarField contraVSumZ(mesh.nCells(), 0);
    forAll(owner, fi) {
	double contraVX = .5*std::fabs(mesh.Sf()[fi].x());
	double contraVY = .5*std::fabs(mesh.Sf()[fi].y());
	double contraVZ = .5*std::fabs(mesh.Sf()[fi].z());
	contraVSumX[owner[fi]] += contraVX;
	contraVSumY[owner[fi]] += contraVY;
	contraVSumZ[owner[fi]] += contraVZ;
	contraVSumX[neighbour[fi]] += contraVX;
	contraVSumY[neighbour[fi]] += contraVY;
	contraVSumZ[neighbour[fi]] += contraVZ;
    }
    volScalarField rPsi(1.0/thermo.psi());
    volScalarField c(sqrt(thermo.Cp()/thermo.Cv()*rPsi));
    scalarField dTLoc
    (
	cfl*mesh.V()/
	(
	    (mag(U.internalField().component(0))+c.internalField())*contraVSumX 
	    + (mag(U.internalField().component(1))+c.internalField())*contraVSumY 
	    + (mag(U.internalField().component(2))+c.internalField())*contraVSumZ 
	)
    );
    for(int i = 0; i < n; ++i) {
	for(int j = 0; j < 5; ++j) {
	    tList.push_back(Triplet(i*5+j, i*5+j, dTLoc[i]));
	}
    }
    SpMat diag(nq, nq);
    diag.setFromTriplets(tList.begin(), tList.end());
    SpMat eye(nq, nq);
    eye.setFromTriplets(eyeList.begin(), eyeList.end());
    scalar minDtLoc = min(dTLoc);
    A = (eye + diag*A);
    b = diag*b;
    //solved by Eigen using new Eigen
//    Eigen::GMRES<SpMat> solver(A);
//    Eigen::GMRES<SpMat, Eigen::IncompleteLUT<double> > solver;
    Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double> > solver;
//    Eigen::BiCGSTAB<SpMat> solver;
    solver.setMaxIterations(100);
    solver.compute(A);
//    solver.analyzePattern(A);
//    solver.factorize(A);
    x = solver.solve(b);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
//    std::cout << "maxIterations:   " << solver.maxIterations() << std::endl;
//    std::cout << "tol: " << solver.tolerance()      << std::endl;
//    cin.get();
    

/*//    Eigen::ConjugateGradient<SpMat, (Eigen::Lower)|(Eigen::Lower), Eigen::IncompleteLUT<double> > cg;
    Eigen::ConjugateGradient<SpMat, (Eigen::Lower)|(Eigen::Lower)> cg;
//    Eigen::ConjugateGradient<SpMat> cg;
//    cg.setMaxIterations(1e4);
    cg.compute(A);
//    cg.analyzePattern(A);
//    cg.factorize(A);
    x = cg.solve(b);
    std::cout << "#iterations:     " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
//    std::cout << "maxIterations:     " << cg.maxIterations() << std::endl;
//    std::cout << "tol: " << cg.tolerance()      << std::endl;
//    cin.get();*/

    //solved by Eigen using inv
//    Eigen::MatrixXd Aden(A);
//    SpMat Ainv(Aden.inverse().sparseView());
//    x = Ainv*b;

    //solved by Eigen using old Eigen
//    Eigen::BiCGSTAB<SpMat> solver;
//    solver.compute(A);

//    solver.analyzePattern(A);
//    solver.factorize(A);
//    x = solver.solve(b);

/*    myMat matV, matSumX, matSumY, matSumZ;
    for(int ni = 0; ni < rhoE.internalField().size(); ni++) {
	matV.push_back(std::vector<double>(1,mesh.V()[ni]));
	matSumX.push_back(std::vector<double>(1,contraVSumX[ni]));
	matSumY.push_back(std::vector<double>(1,contraVSumY[ni]));
	matSumZ.push_back(std::vector<double>(1,contraVSumZ[ni]));
    }
    myMatWrite(matV, "matV");
    myMatWrite(matSumX, "matSumX");
    myMatWrite(matSumY, "matSumY");
    myMatWrite(matSumZ, "matSumZ");
    myMatWriteDen(eye, "eye");
    myMatWriteDen(diag, "diag");
    myMatWriteDen(A, "A");
    Info << "done" << endl;
    cin.get();*/
/*    myMatWriteDen(A, "A");
    myMatWriteDen(Ainv, "Ainv");
    myVecWriteSp(x, "matX");
    myVecWriteSp(b, "matb");
    Info << "done" << endl;
    cin.get();*/

    //export to main; try to speed up this part
    VectorXd x0(nq);
    x0 = matOdr.transpose()*x;

/*    myVecWriteSp(x, "matx");
    myVecWriteSp(b, "matb");
    myMatWriteDen(matJcb, "matJcb");
    Info << "matJcb" << endl;
    myMatWriteDen(matOdr, "matOdr");
    Info << "matOdr" << endl;
    Info << "done0" << endl;*/

    const Foam::Vector<Foam::label> validD = rhoU.mesh().solutionD();

    rho.internalField() -= scalarField(UList<double>(x0.data(), n));
    for(unsigned int cmpt = 0; cmpt < 3; ++cmpt) { 
	if (validD[cmpt] == -1) continue;
	rhoU.internalField().replace
	(
	    cmpt,
	    rhoU.internalField().component(cmpt) -
	    scalarField(UList<double>(x0.data()+(cmpt+1)*n, n))
	);
    }
    rhoE.internalField() -= scalarField(UList<double>(x0.data()+4*n, n));
    
/*    Info << rho << endl;
    Info << rhoU << endl;
    Info << rhoE << endl;
    myVecWriteSp(x, "matX");
    myVecWriteSp(b, "matb");*/
/*    myMatWriteDen(A, "A");
    myVecWriteSp(x, "matX");
    myVecWriteSp(b, "matb");
    myMat matdQdt, matQ_new;
    for(int ni = 0; ni < rhoE.internalField().size(); ni++) {
	    matdQdt.push_back(std::vector<double>(1,phi.internalField()[ni]));
	    matdQdt.push_back(std::vector<double>(1,phiUp.internalField()[ni].x()));
	    matdQdt.push_back(std::vector<double>(1,phiUp.internalField()[ni].y()));
	    matdQdt.push_back(std::vector<double>(1,phiUp.internalField()[ni].z()));
	    matdQdt.push_back(std::vector<double>(1,phiEp.internalField()[ni]));
	    matQ_new.push_back(std::vector<double>(1,rho.internalField()[ni]));
	    matQ_new.push_back(std::vector<double>(1,rhoU.internalField()[ni].x()));
	    matQ_new.push_back(std::vector<double>(1,rhoU.internalField()[ni].y()));
	    matQ_new.push_back(std::vector<double>(1,rhoU.internalField()[ni].z()));
	    matQ_new.push_back(std::vector<double>(1,rhoE.internalField()[ni]));
    }
    myMatWrite(matdQdt, "matdQdt");
    myMatWrite(matQ_new, "matQ_new");

    Info << "done" << endl;
    cin.get();*/

    return 1;
}

int myDdtEulerExp
(
    SpMat& A,
    const fvMesh& mesh
)
{
    TripletList tList;
    int n = mesh.nCells();
    int nq = 5*n;
    tList.reserve(nq);
	
//    scalarField A_coeff = 1.0/mesh.time().deltaTValue()*mesh.Vsc();
    double rdt = 1.0/mesh.time().deltaTValue();
    for(int i = 0; i < nq; ++i) {
//	int vi = i%n;
//	tList.push_back(Triplet(i, i, A_coeff[vi]));
	tList.push_back(Triplet(i, i, rdt));
    }

    A.setFromTriplets(tList.begin(), tList.end());

    return 1;
}

int getJcb 
(
    const surfaceVectorField& U_pos,
    const surfaceScalarField& H_pos,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& H_neg,
    const surfaceVectorField& U_tilt,
    const surfaceScalarField& H_tilt,
    const surfaceScalarField& gamma,
    const volVectorField& U,
    const volScalarField& H,
    const volScalarField& p,
    const volScalarField& T,
    const psiThermo& thermo,
    const fvMesh& mesh,
    SpMat& matJcb
)
{
//Aeq = -Div*0.5*(A_pos*(T_pos+SP2_pos)+A_neg*(T_neg+SP2_neg)-A_abs*((T_neg+SP2_neg)-(T_pos+SP2_pos))) - DivBd*B*(C^(-1)*F*D);
/*    int nqIntCells = 1710;
    int nqBdFaces = 380;
    int nqIntFaces = 3230;*/
    unsigned int nBdFaces = 0;

    forAll (mesh.boundary(), patchi)
    {
	nBdFaces += mesh.boundary()[patchi].size();
    }

    int nqBdFaces =  nBdFaces*5;

    int nqIntCells = (p.internalField().size())*5;
    int nqIntFaces = (mesh.owner().size())*5;
    SpMat matDivBd(nqIntCells, nqBdFaces);
    getMatDivBdSp
    (
	matDivBd,
	nqBdFaces/5,
	mesh
    );

    SpMat matDiv(nqIntCells, nqIntFaces);
    getMatDivSp
    (
	matDiv,
	mesh
    );

    SpMat matA_pos(nqIntFaces, nqIntFaces);
    SpMat matA_neg(nqIntFaces, nqIntFaces);
    SpMat matA_abs(nqIntFaces, nqIntFaces);
    getMatASp
    (
	matA_pos,
	matA_neg,
	matA_abs,
	U_pos,
	H_pos,
	U_neg,
	H_neg,
	U_tilt,
	H_tilt,
	gamma,
	mesh
    );

    SpMat matB(nqBdFaces, nqBdFaces);
    getMatBSp
    (
	matB,
	U,
	H,
	thermo.gamma()(),
	nqBdFaces/5,
	mesh
    );

    SpMat matT_pos(nqIntFaces, nqIntCells);
    SpMat matT_neg(nqIntFaces, nqIntCells);
    SpMat matT2_pos(nqIntFaces, nqIntCells);
    SpMat matT2_neg(nqIntFaces, nqIntCells);
    SpMat matT3_pos(nqIntFaces, nqBdFaces);
    SpMat matT3_neg(nqIntFaces, nqBdFaces);
    Field<scalar> limiterRho_pos;
    Field<scalar> limiterRho_neg;
    Field<vector> limiterRhoU_pos;
    Field<vector> limiterRhoU_neg;
    Field<scalar> limiterRhoE_pos;
    Field<scalar> limiterRhoE_neg;
    getMatTSp
    (
	matT_pos,
	matT_neg,
	matT2_pos,
	matT2_neg,
	matT3_pos,
	matT3_neg,
	limiterRho_pos,
	limiterRho_neg,
	limiterRhoU_pos,
	limiterRhoU_neg,
	limiterRhoE_pos,
	limiterRhoE_neg,
	mesh
    );

    SpMat matC(nqBdFaces, nqBdFaces);
    getMatCSp
    (
	matC,
	thermo.rho()(),
	U,
	T,
	p,
	thermo.gamma()(),
	nqBdFaces/5,
	mesh
    );
    SpMat matD(nqIntCells, nqIntCells);
    getMatDSp
    (
	matD,
	thermo.rho()(),
	U,
	T,
	p,
	thermo.gamma()(),
	nqIntCells/5
    );

    SpMat matF(nqBdFaces, nqIntCells);
    getMatFSp
    (
	matF,
	p,
	U,
	T,
	nqBdFaces/5,
	mesh
    );

//    Eigen::MatrixXd matCden(matC);
//    SpMat matCinv(matCden.inverse().sparseView());
    SpMat matCinv(nqBdFaces, nqBdFaces);
    getMatCinvSp
    (
	matCinv,
	thermo.rho()(),
	U,
	T,
	p,
	thermo.gamma()(),
	nqBdFaces/5,
	mesh
    );

    SpMat matSP2_pos = matT2_pos + (matT3_pos*(matCinv*matF*matD));
    SpMat matSP2_neg = matT2_neg + (matT3_neg*(matCinv*matF*matD));
    
    matJcb = ((matDiv*0.5*(matA_pos*(matT_pos+matSP2_pos)+matA_neg*(matT_neg+matSP2_neg)-matA_abs*((matT_neg+matSP2_neg)-(matT_pos+matSP2_pos)))) + (matDivBd*matB*(matCinv*matF*matD)) );
//    myMatWriteDen(matJcb, "matJcb0");
//    Info << "matJcb0" << endl;
//    matJcb = ((matDiv*0.5*(matA_pos*(matT_pos+matSP2_pos)+matA_neg*(matT_neg+matSP2_neg)-matA_abs*((matT_neg+matSP2_neg)-(matT_pos+matSP2_pos)))) + (matDivBd*matB*(matCinv*matF*matD)) );
/*    myMatWriteDen(matJcb, "matJcb");
    Info << "matJcb" << endl;
    myMatWriteDen(matOdr, "matOdr");
    Info << "matOdr" << endl;
    Info << "done0" << endl;
    cin.get();*/
/*
    myMatWriteDen(matDiv, "matDiv");
    Info << "matDiv" << endl;
    myMatWriteDen(matDivBd, "matDivBd");
    Info << "matDivBd" << endl;
    myMatWriteDen(matA_abs, "matA_abs");
    Info << "matB" << endl;
    myMatWriteDen(matB, "matB");
    Info << "matA_abs" << endl;
    myMatWriteDen(matA_pos, "matA_pos");
    Info << "matA_pos" << endl;
    myMatWriteDen(matA_neg, "matA_neg");
    Info << "matA_neg" << endl;
    myMatWriteDen(matT_pos, "matT_pos");
    Info << "matT_pos" << endl;
    myMatWriteDen(matT_neg, "matT_neg");
    Info << "matT_neg" << endl;
    myMatWriteDen(matT2_pos, "matT2_pos");
    Info << "matT2_pos" << endl;
    myMatWriteDen(matT2_neg, "matT2_neg");
    Info << "matT2_neg" << endl;
    myMatWriteDen(matT3_pos, "matT3_pos");
    Info << "matT3_pos" << endl;
    myMatWriteDen(matT3_neg, "matT3_neg");
    Info << "matT3_neg" << endl;
    myMatWriteDen(matSP2_pos, "matSP2_pos");
    Info << "matSP2_pos" << endl;
    myMatWriteDen(matSP2_neg, "matSP2_neg");
    Info << "matSP2_neg" << endl;
    myMatWriteDen(matC, "matC");
    Info << "matC" << endl;
    myMatWriteDen(matCinv, "matCinv");
    Info << "matCinv" << endl;
    myMatWriteDen(matF, "matF");
    Info << "matF" << endl;
    myMatWriteDen(matD, "matD");
    Info << "matD" << endl;
    myMatWriteDen(matJcb, "matJcb");
    Info << "matJcb" << endl;
    Info << "done0" << endl;
    cin.get();*/
    return 1;
}

int myDdtEulerImp
(
    SpMat& A,
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const psiThermo& thermo,
    const surfaceVectorField& U_pos,
    const surfaceScalarField& H_pos,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& H_neg,
    const surfaceVectorField& U_tilt,
    const surfaceScalarField& H_tilt,
    const surfaceScalarField& gamma,
    double cfl,
    const fvMesh& mesh
)
{
    volScalarField rPsi(1.0/thermo.psi());
    volScalarField H(rhoE/rho+rPsi);

    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    TripletList tList;
    int n = mesh.nCells();
    int nq = 5*n;
    tList.reserve(nq);

    for(int i = 0; i < nq; ++i) {
	tList.push_back(Triplet(i, i, 1));
    }
/*    tList.reserve(nq);

    double rdt = 1.0/mesh.time().deltaTValue();
    for(int i = 0; i < nq; ++i) {
	tList.push_back(Triplet(i, i, rdt));
    }*/
    /////////////////////////////////////////////////
/*    tList.reserve(nq);
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    scalarField contraVSumX(mesh.nCells(), 0);
    scalarField contraVSumY(mesh.nCells(), 0);
    scalarField contraVSumZ(mesh.nCells(), 0);
    forAll(owner, fi) {
	double contraVX = .5*std::fabs(mesh.Sf()[fi].x());
	double contraVY = .5*std::fabs(mesh.Sf()[fi].y());
	double contraVZ = .5*std::fabs(mesh.Sf()[fi].z());
	contraVSumX[owner[fi]] += contraVX;
	contraVSumY[owner[fi]] += contraVY;
	contraVSumZ[owner[fi]] += contraVZ;
	contraVSumX[neighbour[fi]] += contraVX;
	contraVSumY[neighbour[fi]] += contraVY;
	contraVSumZ[neighbour[fi]] += contraVZ;
    }
    volScalarField c(sqrt(thermo.Cp()/thermo.Cv()*rPsi));
    scalarField dTLoc
    (
	cfl*mesh.V()/
	(
	    (mag(U.internalField().component(0))+c.internalField())*contraVSumX 
	    + (mag(U.internalField().component(1))+c.internalField())*contraVSumY 
	    + (mag(U.internalField().component(2))+c.internalField())*contraVSumZ 
	)
    );
    for(int i = 0; i < n; ++i) {
	for(int j = 0; j < 5; ++j) {
	    tList.push_back(Triplet(i*5+j, i*5+j, 1/dTLoc[i]));
//	    tList.push_back(Triplet(i*5+j, i*5+j, 1/1e-5));
	}
    }
    scalar minDtLoc = min(dTLoc);*/
//    Info << dTLoc << endl;
//    Info << minDtLoc << endl;
//    cin.get();

    /////////////////////////////////////////////////
    
    SpMat matJcb(nq, nq);
    getJcb
    (
	U_pos,
	H_pos,
	U_neg,
	H_neg,
	U_tilt,
	H_tilt,
	gamma,
	U,
	H,
	p,
	T,
	thermo,
	mesh,
	matJcb
    );

    SpMat diag(nq, nq);
    diag.setFromTriplets(tList.begin(), tList.end());
//    myMatWriteDen(A, "test1");
//    Info << "done1" << endl;
//    A = (diag);
    A = (matJcb);
//    A = (diag + matJcb);
/*    myMatWriteDen(A, "A");
    Info << "done" << endl;
    cin.get();*/

    return 1;
}

int myDdt
(
    SpMat& A,
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const psiThermo& thermo,
    const surfaceVectorField& U_pos,
    const surfaceScalarField& H_pos,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& H_neg,
    const surfaceVectorField& U_tilt,
    const surfaceScalarField& H_tilt,
    const surfaceScalarField& gamma,
    double cfl,
    const word scheme
)
{
    if(scheme == "EulerExplicit") {
	const fvMesh& mesh = rho.mesh();
	myDdtEulerExp(A, mesh);
    } else if(scheme == "EulerImplicit") {
	const fvMesh& mesh = rho.mesh();
	myDdtEulerImp
	(
	    A,
	    rho,
	    rhoU,
	    rhoE,
	    p,
	    U,
	    T,
	    thermo,
	    U_pos,
	    H_pos,
	    U_neg,
	    H_neg,
	    U_tilt,
	    H_tilt,
	    gamma,
	    cfl,
	    mesh
	);
    } else {
	Info << "SZhao::No available scheme found" << endl;
	cin.get();
    }

    return 1;
}

/*
    typename vector::labelType validComponents
    (
        pow
        (
            rhoU.mesh().solutionD(),
            pTraits<typename powProduct<Vector<label>, vector::rank>::type>::zero
        )
    );

    std::cout << x << std::endl;
    Info << rho.internalField() << endl;
    Info << rhoU.internalField() << endl;
    Info << rhoE.internalField() << endl;
    cin.get();

        value = k.x();
	tList.push_back(Triplet(0, pn+1, value));
	tList.push_back(Triplet(0, nn+1, -value));
        value = k.y();
	tList.push_back(Triplet(0, pn+2, value));
	tList.push_back(Triplet(0, nn+2, -value));
        value = k.z();
	tList.push_back(Triplet(0, pn+3, value));
	tList.push_back(Triplet(0, nn+3, -value));
        value = 0;
	tList.push_back(Triplet(0, pn+4, value));
	tList.push_back(Triplet(0, nn+4, -value));
	// r2
	value = k.x()*phi_sqr - u*theta;
	tList.push_back(Triplet(1, pn+0, value));
	tList.push_back(Triplet(1, nn+0, -value));
        value = theta - k.x()*(gamma-2)*u;
	tList.push_back(Triplet(1, pn+1, value));
	tList.push_back(Triplet(1, nn+1, -value));
        value = k.y()*u - k.x()*(gamma-1)*v;
	tList.push_back(Triplet(1, pn+2, value));
	tList.push_back(Triplet(1, nn+2, -value));
        value = k.z()*u - k.x()*(gamma-1)*w;
	tList.push_back(Triplet(1, pn+3, value));
	tList.push_back(Triplet(1, nn+3, -value));
        value = k.x()*(gamma-1);
	tList.push_back(Triplet(1, pn+4, value));
	tList.push_back(Triplet(1, nn+4, -value));
	// r3
	value = k.y()*phi_sqr - v*theta;
	tList.push_back(Triplet(2, pn+0, value));
	tList.push_back(Triplet(2, nn+0, -value));
        value = k.x()*v - k.y()*(gamma-1)*u;
	tList.push_back(Triplet(2, pn+1, value));
	tList.push_back(Triplet(2, nn+1, -value));
        value = theta - k.y()*(gamma-2)*v;
	tList.push_back(Triplet(2, pn+2, value));
	tList.push_back(Triplet(2, nn+2, -value));
        value = k.z()*v - k.y()*(gamma-1)*w;
	tList.push_back(Triplet(2, pn+3, value));
	tList.push_back(Triplet(2, nn+3, -value));
        value = k.y()*(gamma-1);
	tList.push_back(Triplet(2, pn+4, value));
	tList.push_back(Triplet(2, nn+4, -value));
	// r4                          
	value = k.z()*phi_sqr - w*theta;
	tList.push_back(Triplet(3, pn+0, value));
	tList.push_back(Triplet(3, nn+0, -value));
        value = k.x()*w - k.z()*(gamma-1)*u;
	tList.push_back(Triplet(3, pn+1, value));
	tList.push_back(Triplet(3, nn+1, -value));
        value = k.y()*w - k.z()*(gamma-1)*v;
	tList.push_back(Triplet(3, pn+2, value));
	tList.push_back(Triplet(3, nn+2, -value));
        value = theta - k.z()*(gamma-2)*w;
	tList.push_back(Triplet(3, pn+3, value));
	tList.push_back(Triplet(3, nn+3, -value));
        value = k.z()*(gamma-1);
	tList.push_back(Triplet(3, pn+4, value));
	tList.push_back(Triplet(3, nn+4, -value));
	// r5                          
	value = theta*(phi_sqr-H);
	tList.push_back(Triplet(4, pn+0, value));
	tList.push_back(Triplet(4, nn+0, -value));
        value = k.x()*H - (gamma-1)*u*theta;
	tList.push_back(Triplet(4, pn+1, value));
	tList.push_back(Triplet(4, nn+1, -value));
        value = k.y()*H - (gamma-1)*v*theta;
	tList.push_back(Triplet(4, pn+2, value));
	tList.push_back(Triplet(4, nn+2, -value));
        value = k.z()*H - (gamma-1)*w*theta;
	tList.push_back(Triplet(4, pn+3, value));
	tList.push_back(Triplet(4, nn+3, -value));
        value = gamma*theta;
	tList.push_back(Triplet(4, pn+4, value));
	tList.push_back(Triplet(4, nn+4, -value));
*/
