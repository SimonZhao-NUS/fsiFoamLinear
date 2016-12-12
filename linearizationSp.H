///////////linearization ///////////

unsigned int nIntCells = mesh.nCells();
unsigned int nIntFaces = mesh.owner().size();
unsigned int nBdFaces = 0;

forAll (mesh.boundary(), patchi)
{
    nBdFaces += mesh.boundary()[patchi].size();
}

unsigned int nqIntCells = nIntCells*5;
unsigned int nqIntFaces = nIntFaces*5;
unsigned int nqBdFaces =  nBdFaces*5;

SpMat matT_pos(nqIntFaces, nqIntCells);
SpMat matT_neg(nqIntFaces, nqIntCells);
SpMat matT2_pos(nqIntFaces, nqIntCells);
SpMat matT2_neg(nqIntFaces, nqIntCells);
SpMat matT3_pos(nqIntFaces, nqBdFaces);
SpMat matT3_neg(nqIntFaces, nqBdFaces);
getMatTSp
(
    matT_pos,
    matT_neg,
    matT2_pos,
    matT2_neg,
    matT3_pos,
    matT3_neg,
    rhoLimiter,
    rhoULimiter,
    rhoELimiter,
    mesh
);

SpMat matL_pos(nqIntFaces, nqIntCells);
SpMat matL_neg(nqIntFaces, nqIntCells);
SpMat matLB_pos(nqIntFaces, nqBdFaces);
SpMat matLB_neg(nqIntFaces, nqBdFaces);
getMatLSp
(
    matL_pos,
    matL_neg,
    matLB_pos,
    matLB_neg,
    rhoLimiter,
    rhoULimiter,
    rhoELimiter,
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
    gammaSurface,
    mesh
);

SpMat matC(nqBdFaces, nqBdFaces);
getMatCSp
(
    matC,
    rho,
    U,
    T,
    p,
    gamma,
    nBdFaces,
    mesh
);
SpMat matD(nqIntCells, nqIntCells);
getMatDSp
(
    matD,
    rho,
    U,
    T,
    p,
    gamma,
    nIntCells
);

SpMat matB(nqBdFaces, nqBdFaces);
getMatBSp
(
    matB,
    U,
//    H,
    h+.5*magSqr(U),
    gamma,
    nBdFaces,
    mesh
);

SpMat matBf(nqBdFaces, 4);
getMatBfSp
(
    matBf,
    rho,
    U,
    p,
    nBdFaces,
    mesh
);

SpMat matBdyn1(nqIntFaces, 4);
getMatBdyn1
(
    matBdyn1,
    rho_pos,
    U_pos,
    H_pos,
    p_pos,
    rho_neg,
    U_neg,
    H_neg,
    p_neg
);

SpMat matF(nqBdFaces, nqIntCells);
getMatFSp
(
    matF,
    rho,
    U,
    p,
    sqrt(gamma*p/rho),
    nBdFaces,
    mesh
);

#include "readFarFieldValue.H"
double U_freestream = 0;
double c_freestream = 0;
double rho_freestream = 0;
double p_freestream = 0;
double Mach_freestream = 0;

U_freestream = mag(U_inf);
rho_freestream = rho_inf;
p_freestream = p_inf;
c_freestream = std::sqrt(gamma[0]*p_inf/rho_inf);
Mach_freestream = U_freestream/c_freestream;

Info <<"Free stream veloctiy magnitude:" << U_freestream << endl;
Info <<"Free stream density:" << rho_freestream << endl;
Info <<"Free stream pressure:" << p_freestream << endl;
Info <<"Free stream speed of sound:" << c_freestream << endl;
Info <<"Free stream Mach:" << Mach_freestream << endl;
myMat matFreestream, matFreestream_dim;
matFreestream.push_back(std::vector<double>(1, U_freestream));
matFreestream.push_back(std::vector<double>(1, c_freestream));
matFreestream.push_back(std::vector<double>(1, rho_freestream));
matFreestream.push_back(std::vector<double>(1, p_freestream));
/*
matFreestream_dim.push_back(std::vector<double>(1, cSqr_freestream_dim));
matFreestream_dim.push_back(std::vector<double>(1, rho_freestream_dim));
matFreestream_dim.push_back(std::vector<double>(1, p_freestream_dim));
matFreestream_dim.push_back(std::vector<double>(1, U_freestream_dim));
matFreestream_dim.push_back(std::vector<double>(1, T_freestream_dim));*/

SpMat matH(nqBdFaces, 1);
getMatHSp
(
    matH,
    rho,
    U,
    U_freestream,
    nBdFaces,
    mesh
);

SpMat matK(2, nqIntCells);
getMatKSp
(
    matK,
    mesh
);

int nWallSf = 0;
forAll(mesh.boundary(), patchi) {
    if (mesh.boundary()[patchi].name() == "wing") {
	nWallSf = mesh.boundary()[patchi].size();
    }
}

SpMat matKp(nWallSf, nqIntCells);
getMatKpSp
(
    matKp,
    mesh
);

SpMat matCf(1, nWallSf);
getMatCfSp
(
    matCf,
    mesh
);

SpMat matSf(1, nWallSf);
getMatSfSp
(
    matSf,
    mesh
);

SpMat matDivBd(nqIntCells, nqBdFaces);
getMatDivBdSp
(
    matDivBd,
    nBdFaces,
    mesh
);

SpMat matDiv(nqIntCells, nqIntFaces);
getMatDivSp
(
    matDiv,
    mesh
);

SpMat matQcSp(nqIntCells,1);
TripletList QcSpList;
QcSpList.reserve(nqIntCells);
myMat matG, matQc, matFsum_ben, matQf_pos, matQf_neg;
for(unsigned int ni = 0; ni < nIntCells; ni++) {
    matG.push_back(std::vector<double>(1,p.internalField()[ni]));
    matG.push_back(std::vector<double>(1,U.internalField()[ni].x()));
    matG.push_back(std::vector<double>(1,U.internalField()[ni].y()));
    matG.push_back(std::vector<double>(1,U.internalField()[ni].z()));
    matG.push_back(std::vector<double>(1,T.internalField()[ni]));
    
/*    matQ.push_back(std::vector<double>(1,rho.internalField()[ni]));
    matQ.push_back(std::vector<double>(1,rhoU.internalField()[ni].x()));
    matQ.push_back(std::vector<double>(1,rhoU.internalField()[ni].y()));
    matQ.push_back(std::vector<double>(1,rhoU.internalField()[ni].z()));
    matQ.push_back(std::vector<double>(1,rhoE.internalField()[ni]));*/
    
    matQc.push_back(std::vector<double>(1,0));
    matQc.push_back(std::vector<double>(1,-U_freestream*std::sin(0/57.296)*rho.internalField()[ni]));
    matQc.push_back(std::vector<double>(1,U_freestream*std::cos(0/57.296)*rho.internalField()[ni]));
    matQc.push_back(std::vector<double>(1,0));
    matQc.push_back(std::vector<double>(1,0));                                                              
    QcSpList.push_back(Triplet(ni*5+1, 0,-U_freestream*std::sin(0/57.296)*rho.internalField()[ni]));
    QcSpList.push_back(Triplet(ni*5+2, 0,U_freestream*std::cos(0/57.296)*rho.internalField()[ni]));

    matFsum_ben.push_back(std::vector<double>(1,(fvc::div(phi))().internalField()[ni]));
    matFsum_ben.push_back(std::vector<double>(1,(fvc::div(phiUp))().internalField()[ni].x()));
    matFsum_ben.push_back(std::vector<double>(1,(fvc::div(phiUp))().internalField()[ni].y()));
    matFsum_ben.push_back(std::vector<double>(1,(fvc::div(phiUp))().internalField()[ni].z()));
    matFsum_ben.push_back(std::vector<double>(1,(fvc::div(phiEp))().internalField()[ni]));
}
matQcSp.setFromTriplets(QcSpList.begin(), QcSpList.end());
myMatWrite(matFsum_ben, "matFsum_ben");

myMat matGb;
forAll(mesh.boundary(), patchi) {
    const fvPatchField<scalar>& pb = p.boundaryField()[patchi];
    const fvPatchField<vector>& Ub = U.boundaryField()[patchi];
    const fvPatchField<scalar>& Tb = T.boundaryField()[patchi];
    forAll(mesh.boundary()[patchi], fi) {
	matGb.push_back(std::vector<double>(1,pb[fi]));
	matGb.push_back(std::vector<double>(1,Ub[fi].x()));
	matGb.push_back(std::vector<double>(1,Ub[fi].y()));
	matGb.push_back(std::vector<double>(1,Ub[fi].z()));
	matGb.push_back(std::vector<double>(1,Tb[fi]));    
    }
}
myMatWrite(matGb, "matGb");

myMat wallCenter
(
    1,
    std::vector<double>(nIntCells, 0)
);

int nBdIndex = 0;
forAll(mesh.boundary(), patchi) {
    const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
    forAll(mesh.boundary()[patchi], facei) {
	if (mesh.boundary()[patchi].name() == "wing") {
	    wallCenter[0][faceCells[facei]] = mesh.Cf().boundaryField()[patchi][facei].x();
	}

	nBdIndex++;
    }
}

SpMat matE(nqBdFaces, 1);
getMatESp
(
    matE,
    rho,
    p,
    U,
    T,
    nBdFaces,
    mesh,
    U_freestream
);

SpMat matQic(nqIntCells, 1), matQb(nqBdFaces, 1), matQif_pos(nqIntFaces, 1), matQif_neg(nqIntFaces, 1);
getMatQsSp
(
    matQic,
    matQif_pos,
    matQif_neg,
    matQb,
    rho,
    rhoU,
    rhoE,
    rho_pos,
    rho_neg,
    rhoU_pos,
    rhoU_neg,
    rhoE_pos,
    rhoE_neg
);

Info << "Writing matrices, wait..." << endl;
/*myMatWriteDenSp(matA_pos, "matA_pos");
myMatWriteDenSp(matA_neg, "matA_neg");
//myMatWriteDenSp(matA_abs_pos, "matA_abs_pos");
//myMatWriteDenSp(matA_abs_neg, "matA_abs_neg");
myMatWriteDenSp(matA_abs, "matA_abs");
myMatWriteDenSp(matT_pos, "matT_pos");
myMatWriteDenSp(matT_neg, "matT_neg");
myMatWriteDenSp(matT2_pos, "matT2_pos");
myMatWriteDenSp(matT2_neg, "matT2_neg");
myMatWriteDenSp(matT3_pos, "matT3_pos");
myMatWriteDenSp(matT3_neg, "matT3_neg");
myMatWriteDenSp(matB, "matB");
myMatWriteDenSp(matC, "matC");
myMatWriteDenSp(matF, "matF");
myMatWriteDenSp(matD, "matD");
myMatWriteDenSp(matDiv, "matDiv");
myMatWriteDenSp(matDivBd, "matDivBd");
myMatWriteDenSp(matH, "matH");
myMatWriteDenSp(matK, "matK");
myMatWrite(wallCenter, "matCenter");
myMatWrite(matG, "matG");
myMatWrite(matQ, "matQ");
myMatWrite(matQc, "matQc");
myMatWriteDenSp(matE, "matE");
myMatWrite(matFreestream, "matFreestream");
myMatWrite(matFreestream_dim, "matFreestream_dim");
myMatWriteDenSp(matQic, "matQic");
myMatWriteDenSp(matQb, "matQb");
myMatWriteDenSp(matQif_pos, "matQif_pos");
myMatWriteDenSp(matQif_neg, "matQif_neg");
myMatWriteDenSp(matLimiter_pos, "matLimiter_pos");
myMatWriteDenSp(matLimiter_neg, "matLimiter_neg");*/

//Eigen::MatrixXd matCden(matC);
//SpMat matCinv(matCden.inverse().sparseView());
SpMat matCinv(nqBdFaces, nqBdFaces);
getMatCinvSp
(
    matCinv,
    rho,
    U,
    T,
    p,
    gamma,
    nqBdFaces/5,
    mesh
);

SpMat matSP2_pos = matT2_pos + (matT3_pos*(matCinv*matF*matD));
SpMat matSP2_neg = matT2_neg + (matT3_neg*(matCinv*matF*matD));

////
TripletList wallList;
int nWall = 0;
forAll(mesh.boundary(), patchi) {
    const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
    forAll(mesh.boundary()[patchi], facei) {
	if (mesh.boundary()[patchi].name() == "wing") {
	    for(int i = 0; i < 5; ++i) {
		wallList.push_back(Triplet(nWall, 0, faceCells[facei]*5+i+1));
		++nWall;
	    }
	}
    }
}
SpMat matWall(nWall,1);
matWall.setFromTriplets(wallList.begin(), wallList.end());
myMatWriteDenSp(matWall, "matWall");
Info << nWall << endl;
////
TripletList xWallList, pWallList;
int nW = 0;
forAll(mesh.boundary(), patchi) {
    const labelUList &faceCells = mesh.boundary()[patchi].faceCells();
    const vectorField &cp = mesh.Cf().boundaryField()[patchi];
    forAll(mesh.boundary()[patchi], facei) {
	if (mesh.boundary()[patchi].name() == "wing") {
	    xWallList.push_back(Triplet(nW, 0, cp[facei].x()));
	    pWallList.push_back(Triplet(nW, 0, p.internalField()[faceCells[facei]]));
	    ++nW;
	}
    }
}
SpMat matxWall(nW,1);
SpMat matpWall(nW,1);
matxWall.setFromTriplets(xWallList.begin(), xWallList.end());
matpWall.setFromTriplets(pWallList.begin(), pWallList.end());
myMatWriteDenSp(matxWall, "matxWall");
myMatWriteDenSp(matpWall, "matpWall");
////2nd Jacobian
//SpMat matAeq = -((matDiv*0.5*(matA_pos*(matT_pos+matSP2_pos)+matA_neg*(matT_neg+matSP2_neg)-matA_abs*((matT_neg+matSP2_neg)-(matT_pos+matSP2_pos)))) + (matDivBd*matB*(matCinv*matF*matD)));

SpMat matAeq = -((matDiv*0.5*(matA_pos*(matT_pos+matSP2_pos)+matA_neg*(matT_neg+matSP2_neg)-matA_abs*((matT_neg+matSP2_neg)-(matT_pos+matSP2_pos)))) + (matDivBd*matB*(matCinv*matF*matD)));

SpMat matAeq_inner = -(matDiv*0.5*(matA_pos*(matT_pos+matSP2_pos)+matA_neg*(matT_neg+matSP2_neg)-matA_abs*((matT_neg+matSP2_neg)-(matT_pos+matSP2_pos))));
SpMat matAeq_bd = -(matDivBd*matB*(matCinv*matF*matD));

SpMat matBeq =  -matDivBd*matBf;
SpMat matBeq_dyn1 =  -matDiv*matBdyn1;

SpMat matCeq = matK*matD;

SpMat matCp = matKp*matD;

myMatWriteDenSp(matAeq, "matAeq");
myMatWriteDenSp(matA_pos, "matA_pos");
myMatWriteDenSp(matA_neg, "matA_neg");
myMatWriteDenSp(matA_abs, "matA_abs");
myMatWriteDenSp(matAeq_inner, "matAeq_inner");
myMatWriteDenSp(matAeq_bd, "matAeq_bd");
myMatWriteDenSp(matBeq, "matBeq");
myMatWriteDenSp(matCeq, "matCeq");
myMatWriteDenSp(matCp, "matCp");
myMatWriteDenSp(matBeq_dyn1, "matBeq_dyn1");
//myMatWriteDenSp(matBeq_0, "matBeq_0");
//myMatWriteDenSp(matBeq_1, "matBeq_1");
//myMatWriteDenSp(matBeq_2, "matBeq_2");
myMatWriteDenSp(matQcSp, "matQc");
myMatWriteDenSp(matB, "matB");
myMatWriteDenSp(matE, "matE");
myMatWriteDenSp(matCf, "matCf");
myMatWriteDenSp(matK, "matK");
myMatWriteDenSp(matKp, "matKp");
myMatWriteDenSp(matSf, "matSf");
myMatWriteDenSp(matT_pos, "matT_pos");
myMatWriteDenSp(matT_neg, "matT_neg");
myMatWriteDenSp(matT2_pos, "matT2_pos");
myMatWriteDenSp(matT2_neg, "matT2_neg");
myMatWriteDenSp(matT3_pos, "matT3_pos");
myMatWriteDenSp(matT3_neg, "matT3_neg");
myMatWriteDenSp(matL_pos, "matL_pos");
myMatWriteDenSp(matL_neg, "matL_neg");
myMatWriteDenSp(matLB_pos, "matLB_pos");
myMatWriteDenSp(matLB_neg, "matLB_neg");
myMatWriteDenSp(matSP2_pos, "matSP2_pos");
myMatWriteDenSp(matSP2_neg, "matSP2_neg");
myMatWriteDenSp(matCinv, "matCinv");
myMatWriteDenSp(matF, "matF");
myMatWriteDenSp(matD, "matD");
myMatWriteDenSp(matDiv, "matDiv");
myMatWriteDenSp(matDivBd, "matDivBd");
SpMat matS(nqBdFaces,1);
getMatSSp
(
    matS,
    mesh
);

myMatWriteDenSp(matS, "matS");
//myMatWriteDenSp(matQcSp, "matQcSp");
//myMatWriteDen(matBeq, "matBeq_d");

Info << "Done, martrices written" << endl;

//break;
////////linearization end////////////////////	

