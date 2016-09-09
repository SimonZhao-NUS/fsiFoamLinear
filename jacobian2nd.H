typedef std::vector<std::vector<std::vector<double> > >  myTensor;
#include "jacobian2ndFuncs.H"

myMat vec2DiagMat
(
    myMat diag
)
{
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );

    //r1
    mat[0][0] = diag[0][0];
    mat[0][1] = 0;
    mat[0][2] = 0;
    mat[0][3] = 0;
    mat[0][4] = 0;
    //r2
    mat[1][0] = 0;
    mat[1][1] = diag[1][0];
    mat[1][2] = 0;
    mat[1][3] = 0;
    mat[1][4] = 0;
    //r3
    mat[2][0] = 0;
    mat[2][1] = 0;
    mat[2][2] = diag[2][0];
    mat[2][3] = 0;
    mat[2][4] = 0;
    //r4
    mat[3][0] = 0;
    mat[3][1] = 0;
    mat[3][2] = 0;
    mat[3][3] = diag[3][0];
    mat[3][4] = 0;
    //r5
    mat[4][0] = 0;
    mat[4][1] = 0;
    mat[4][2] = 0;
    mat[4][3] = 0;
    mat[4][4] = diag[4][0];

    return mat;
}

myMat transpose
(
    myMat matT
)
{
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    //r1
    mat[0][0] = matT[0][0];
    mat[0][1] = matT[1][0];
    mat[0][2] = matT[2][0];
    mat[0][3] = matT[3][0];
    mat[0][4] = matT[4][0];
    //r2
    mat[1][0] = matT[0][1];
    mat[1][1] = matT[1][1];
    mat[1][2] = matT[2][1];
    mat[1][3] = matT[3][1];
    mat[1][4] = matT[4][1];
    //r3
    mat[2][0] = matT[0][2];
    mat[2][1] = matT[1][2];
    mat[2][2] = matT[2][2];
    mat[2][3] = matT[3][2];
    mat[2][4] = matT[4][2];
    //r4
    mat[3][0] = matT[0][3];
    mat[3][1] = matT[1][3];
    mat[3][2] = matT[2][3];
    mat[3][3] = matT[3][3];
    mat[3][4] = matT[4][3];
    //r5
    mat[4][0] = matT[0][4];
    mat[4][1] = matT[1][4];
    mat[4][2] = matT[2][4];
    mat[4][3] = matT[3][4];
    mat[4][4] = matT[4][4];

    return mat;
}
myMat colCmb
(
    myMat col0,
    myMat col1,
    myMat col2,
    myMat col3,
    myMat col4
)
{
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = col0[0][0];
    mat[0][1] = col1[0][0];
    mat[0][2] = col2[0][0];
    mat[0][3] = col3[0][0];
    mat[0][4] = col4[0][0];
    //r2
    mat[1][0] = col0[1][0];
    mat[1][1] = col1[1][0];
    mat[1][2] = col2[1][0];
    mat[1][3] = col3[1][0];
    mat[1][4] = col4[1][0];
    //r3
    mat[2][0] = col0[2][0];
    mat[2][1] = col1[2][0];
    mat[2][2] = col2[2][0];
    mat[2][3] = col3[2][0];
    mat[2][4] = col4[2][0];
    //r4
    mat[3][0] = col0[3][0];
    mat[3][1] = col1[3][0];
    mat[3][2] = col2[3][0];
    mat[3][3] = col3[3][0];
    mat[3][4] = col4[3][0];
    //r5
    mat[4][0] = col0[4][0];
    mat[4][1] = col1[4][0];
    mat[4][2] = col2[4][0];
    mat[4][3] = col3[4][0];
    mat[4][4] = col4[4][0];
    

    return mat;
}
std::vector<double> outterProd(scalar a, std::vector<double> vec1)
{
    unsigned int nr1 = vec1.size();

    std::vector<double> vec(nr1, 0);
    
    for (unsigned int r = 0; r < nr1; ++r) {
	    vec[r] = a*vec1[r];
    }	

    return vec;
}

myMat outterProd(scalar a, myMat mat1)
{
    unsigned int nr1 = mat1.size();
    unsigned int nc1 = mat1[0].size();

    myMat mat
    (
	nr1,
	std::vector<double>(nc1, 0)
    );
    
    for (unsigned int r = 0; r < nr1; ++r) {
	for (unsigned int c = 0; c < nc1; ++c) {
	    mat[r][c] = a*mat1[r][c];
	}
    }	

    return mat;
}

myMat add(myMat mat1, myMat mat2)
{
    unsigned int nr1 = mat1.size();
    unsigned int nc1 = mat1[0].size();
    unsigned int nr2 = mat2.size();
    unsigned int nc2 = mat2[0].size();

    myMat mat
    (
	nr1,
	std::vector<double>(nc1, 0)
    );
    
    if (nr1 == nr2 && nc1 == nc2 ) {
	for (unsigned int r = 0; r < nr1; ++r) {
	    for (unsigned int c = 0; c < nc1; ++c) {
		mat[r][c] = mat1[r][c] + mat2[r][c];
	    }
	}	
    }

    return mat;
}

myMat calcRoeDeriv()
{
    myMat mat;


    return mat;
}


myMat calcRoeDamping
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar H,
    const scalar rho_pos,
    const scalar rho_neg,
    const vector rhoU_pos,
    const vector rhoU_neg,
    const scalar rhoE_pos,
    const scalar rhoE_neg,
    const scalar gamma
)
{
    myMat jcb_roe ( 
	innProd(
	    calcEigVecR(sf, U, H, gamma),
	    innProd(
		calcEigVal_abs(sf, U, H, gamma),
		calcEigVecL(sf, U, H, gamma)
	    )
	)
    );
    myMat test ( 
	innProd(
	    calcEigVecR(sf, U, H, gamma),
	    calcEigVecL(sf, U, H, gamma)
	)
    );

    myPrintMat(jcb_roe);
    cin.get();
    myPrintMat(test);
    cin.get();

    myMat deltaW
    (
	5,
	std::vector<double>(1,0)
    );
    
    deltaW[0][0]= rho_neg - rho_pos;
    deltaW[1][0]= rhoU_neg.x() - rhoU_pos.x();
    deltaW[2][0]= rhoU_neg.y() - rhoU_pos.y();
    deltaW[3][0]= rhoU_neg.z() - rhoU_pos.z();
    deltaW[4][0]= rhoE_neg - rhoE_pos;

    myPrintMat(deltaW);
    myMat roeDamping 
    ( 
	innProd(
	    jcb_roe, 
	    deltaW
	)
    );

    myPrintMat(roeDamping);
    cin.get();
    Info<< "print end"<<endl;
    cin.get();

    return roeDamping;
}

myMat dM_bardW_pos
(
    const scalar rho_pos,
    const scalar rho_neg,
    const vector U_pos,
    const vector U_neg,
    const scalar H_pos,
    const scalar H_neg,
    const scalar gamma
)
{
    const scalar u_pos = U_pos.x();
    const scalar v_pos = U_pos.y();
    const scalar w_pos = U_pos.z();
    const scalar u_neg = U_neg.x();
    const scalar v_neg = U_neg.y();
    const scalar w_neg = U_neg.z();
     
    const scalar rho_sqrt_pos = std::sqrt(rho_pos);
    const scalar rho_sqrt_neg = std::sqrt(rho_neg);
    const scalar common1 = 1/(rho_sqrt_pos+rho_sqrt_neg);
    const scalar common2 = common1/(2*rho_sqrt_pos);
    const scalar q2_pos = magSqr(U_pos);
    const scalar H_roe = (rho_sqrt_pos*H_pos +rho_sqrt_neg*H_neg)*common1;

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = rho_sqrt_neg/(2*rho_sqrt_pos);
    mat[0][1] = 0; 
    mat[0][2] = 0;
    mat[0][3] = 0;
    mat[0][4] = 0; 
    //r2
    mat[1][0] = common2*(u_pos - (rho_sqrt_pos*u_pos+rho_sqrt_neg*u_neg)*common1);
    mat[1][1] = common2*2*rho_pos;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = 0; 
    //r3
    mat[2][0] = common2*(v_pos - (rho_sqrt_pos*v_pos+rho_sqrt_neg*v_neg)*common1);
    mat[2][1] = 0;
    mat[2][2] = common2*2*rho_pos;           
    mat[2][3] = 0; 
    mat[2][4] = 0; 
    //r4
    mat[3][0] = common2*(w_pos - (rho_sqrt_pos*w_pos+rho_sqrt_neg*w_neg)*common1);
    mat[3][1] = 0;
    mat[3][2] = 0; 
    mat[3][3] = common2*2*rho_pos; 
    mat[3][4] = 0; 
    //r5
    mat[4][0] = common2*(q2_pos-H_pos-H_roe);
    mat[4][1] = common2*(2*rho_pos*u_pos);                 
    mat[4][2] = common2*(2*rho_pos*v_pos);                 
    mat[4][3] = common2*(2*rho_pos*w_pos);     
    mat[4][4] = common2*(2*gamma/(gamma-1)); 
    return mat;
}

myMat dM_bardW_neg
(
    const scalar rho_pos,
    const scalar rho_neg,
    const vector U_pos,
    const vector U_neg,
    const scalar H_pos,
    const scalar H_neg,
    const scalar gamma
)
{
    return dM_bardW_pos
	   (
	       rho_neg,
	       rho_pos,
	       U_neg,
	       U_pos,
	       H_neg,
	       H_pos,
	       gamma
	   );
}

myMat dMdM_bar
(
    const scalar rho,
    const vector U,
    const scalar H,
    const scalar gamma
)
{
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();

    const scalar gamma_p = (gamma -1)/(gamma);
    const scalar q2 = magSqr(U);
     
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = 1;
    mat[0][1] = 0; 
    mat[0][2] = 0;
    mat[0][3] = 0;
    mat[0][4] = 0; 
    //r2
    mat[1][0] = 0;
    mat[1][1] = 1;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = 0; 
    //r3
    mat[2][0] = 0;
    mat[2][1] = 0;
    mat[2][2] = 1;
    mat[2][3] = 0; 
    mat[2][4] = 0; 
    //r4
    mat[3][0] = 0;
    mat[3][1] = 0;
    mat[3][2] = 0; 
    mat[3][3] = 1;
    mat[3][4] = 0; 
    //r5
    mat[4][0] = (H-.5*q2)*gamma_p;
    mat[4][1] = -rho*u*gamma_p;
    mat[4][2] = -rho*v*gamma_p;
    mat[4][3] = -rho*w*gamma_p;
    mat[4][4] = rho*gamma_p; 

    return mat;
}

myMat dWdQ
(
    const scalar rho,
    const vector U,
    const scalar gamma
)
{
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();

    const scalar gamma_s = gamma -1;
    const scalar q2 = magSqr(U);
     
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = 1;
    mat[0][1] = 0; 
    mat[0][2] = 0;
    mat[0][3] = 0;
    mat[0][4] = 0; 
    //r2
    mat[1][0] = -u/rho;
    mat[1][1] = 1/rho;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = 0; 
    //r3
    mat[2][0] = -v/rho;
    mat[2][1] = 0;
    mat[2][2] = 1/rho;
    mat[2][3] = 0; 
    mat[2][4] = 0; 
    //r4
    mat[3][0] = -w/rho;
    mat[3][1] = 0;
    mat[3][2] = 0; 
    mat[3][3] = 1/rho;
    mat[3][4] = 0; 
    //r5
    mat[4][0] = .5*gamma_s*q2;
    mat[4][1] = -u*gamma_s;
    mat[4][2] = -v*gamma_s;
    mat[4][3] = -w*gamma_s;
    mat[4][4] = gamma_s; 

    return mat;
}

myMat dMdQ_pos
(
    const scalar rho_tilt,
    const vector U_tilt,
    const scalar H_tilt,
    const scalar rho_pos,
    const scalar rho_neg,
    const vector U_pos,
    const vector U_neg,
    const scalar H_pos,
    const scalar H_neg,
    const scalar gamma
)
{
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    myMat tmp =
	innProd
	(
	    dMdM_bar
	    (
		rho_tilt,
		U_tilt,
		H_tilt,
		gamma
	    ),
	    dM_bardW_pos
	    (
		rho_pos,
		rho_neg,
		U_pos,
		U_neg,
		H_pos,
		H_neg,
		gamma
	    )
	);
    mat = 
	innProd
	(
	    tmp,
	    dWdQ
	    (
		rho_pos,
		U_pos,
		gamma
	    )
	);
    return mat;
}

myTensor dP_invdQ_pos
(
    const vector sf,
    const scalar rho_tilt,
    const vector U_tilt,
    const scalar H_tilt,
    const scalar rho_pos,
    const scalar rho_neg,
    const vector U_pos,
    const vector U_neg,
    const scalar H_pos,
    const scalar H_neg,
    const scalar gamma
)
{
    const scalar p_tilt = (gamma-1)/gamma*rho_tilt*(H_tilt-.5*magSqr(U_tilt));
    myTensor tensor
    (
	5,
	myMat
	(
	    5,
	    std::vector<double>(5,0)
	)
    );

    myMat matdMdQ_pos = 
	    dMdQ_pos(
		rho_tilt,
		U_tilt,
		H_tilt,
		rho_pos,
		rho_neg,
		U_pos,
		U_neg,
		H_pos,
		H_neg,
		gamma
	    );

    tensor[0] = 
	innProd
	(
	    dR1dM(sf, U_tilt),
	    matdMdQ_pos
	);
    tensor[1] = 
	innProd
	(
	    dR2dM(sf, U_tilt),
	    matdMdQ_pos
	);
    tensor[2] = 
	innProd
	(
	    dR3dM(sf, U_tilt),
	    matdMdQ_pos
	);
    tensor[3] = 
	innProd
	(
	    dR4dM(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);
    tensor[4] = 
	innProd
	(
	    dR5dM(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);

    return tensor;
}

myMat dLambda_absdM
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const scalar sf1 = sf.x();
    const scalar sf2 = sf.y();
    const scalar sf3 = sf.z();
    const scalar cSf = std::sqrt(gamma*p/rho)*mag(sf);

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );

    scalar sign1 = sign(U&sf);
    scalar sign3 = sign((U&sf)+cSf);
    scalar sign4 = sign((U&sf)-cSf);

    //r1
    mat[0][0] = 0;
    mat[0][1] = sign1*sf1; 
    mat[0][2] = sign1*sf2;
    mat[0][3] = sign1*sf3;
    mat[0][4] = 0; 
    //r2
    mat[1][0] = 0; 
    mat[1][1] = sign1*sf1;
    mat[1][2] = sign1*sf2;
    mat[1][3] = sign1*sf3;
    mat[1][4] = 0; 
    //r3
    mat[2][0] = 0; 
    mat[2][1] = sign1*sf1;
    mat[2][2] = sign1*sf2;
    mat[2][3] = sign1*sf3;
    mat[2][4] = 0; 
    //r4
    mat[3][0] = sign3*(-cSf/(2*rho));
    mat[3][1] = sign3*sf1;        
    mat[3][2] = sign3*sf2;        
    mat[3][3] = sign3*sf3;        
    mat[3][4] = sign3*cSf/(2*p); 
    //r5
    mat[4][0] = sign4*cSf/(2*rho);
    mat[4][1] = sign4*sf1;
    mat[4][2] = sign4*sf2; 
    mat[4][3] = sign4*sf3; 
    mat[4][4] = sign4*(-cSf/(2*p));
    
    return mat;
}

myMat dLambda_absdQ_pos
(
    const vector sf,
    const scalar rho_tilt,
    const vector U_tilt,
    const scalar H_tilt,
    const scalar rho_pos,
    const scalar rho_neg,
    const vector U_pos,
    const vector U_neg,
    const scalar H_pos,
    const scalar H_neg,
    const scalar gamma
)
{
    const scalar p_tilt = (gamma-1)/gamma*rho_tilt*(H_tilt-.5*magSqr(U_tilt));
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );

    myMat matdLambda_absdM
    (
	dLambda_absdM
	(
	    sf,
	    rho_tilt,
	    U_tilt,
	    p_tilt,
	    gamma
	)
    );

    myMat matdMdQ_pos = 
	    dMdQ_pos(
		rho_tilt,
		U_tilt,
		H_tilt,
		rho_pos,
		rho_neg,
		U_pos,
		U_neg,
		H_pos,
		H_neg,
		gamma
	    );
    mat = 
	innProd
	(
	    matdLambda_absdM,
	    matdMdQ_pos
	);

    return mat;
}

myTensor dPdQ_pos
(
    const vector sf,
    const scalar rho_tilt,
    const vector U_tilt,
    const scalar H_tilt,
    const scalar rho_pos,
    const scalar rho_neg,
    const vector U_pos,
    const vector U_neg,
    const scalar H_pos,
    const scalar H_neg,
    const scalar gamma
)
{
    const scalar p_tilt = (gamma-1)/gamma*rho_tilt*(H_tilt-.5*magSqr(U_tilt));
    myTensor tensor
    (
	5,
	myMat
	(
	    5,
	    std::vector<double>(5,0)
	)
    );

    myMat matdMdQ_pos = 
	    dMdQ_pos(
		rho_tilt,
		U_tilt,
		H_tilt,
		rho_pos,
		rho_neg,
		U_pos,
		U_neg,
		H_pos,
		H_neg,
		gamma
	    );

    tensor[0] = 
	innProd
	(
	    dL1dM(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);
    tensor[1] = 
	innProd
	(
	    dL2dM(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);
    tensor[2] = 
	innProd
	(
	    dL3dM(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);
    tensor[3] = 
	innProd
	(
	    dL4dM(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);
    tensor[4] = 
	innProd
	(
	    dL5dM(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);

    return tensor;
}
myTensor dPdQ_posT
(
    const vector sf,
    const scalar rho_tilt,
    const vector U_tilt,
    const scalar H_tilt,
    const scalar rho_pos,
    const scalar rho_neg,
    const vector U_pos,
    const vector U_neg,
    const scalar H_pos,
    const scalar H_neg,
    const scalar gamma
)
{
    const scalar p_tilt = (gamma-1)/gamma*rho_tilt*(H_tilt-.5*magSqr(U_tilt));
    myTensor tensor
    (
	5,
	myMat
	(
	    5,
	    std::vector<double>(5,0)
	)
    );

    myMat matdMdQ_pos = 
	    dMdQ_pos(
		rho_tilt,
		U_tilt,
		H_tilt,
		rho_pos,
		rho_neg,
		U_pos,
		U_neg,
		H_pos,
		H_neg,
		gamma
	    );

    tensor[0] = 
	innProd
	(
	    dL1dMT(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);
    tensor[1] = 
	innProd
	(
	    dL2dMT(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);
    tensor[2] = 
	innProd
	(
	    dL3dMT(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);
    tensor[3] = 
	innProd
	(
	    dL4dMT(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);
    tensor[4] = 
	innProd
	(
	    dL5dMT(sf, rho_tilt, U_tilt, p_tilt, gamma),
	    matdMdQ_pos
	);

    return tensor;
}

myMat tensorInnerProd3T
(
    myMat P_invLambda_abs,
    myTensor tensordPdQ_pos,
    myMat deltaQ
)
{
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );

    myMat tmp
    (
	5,
	std::vector<double>(5,0)
    );
    
    scalar a0 = deltaQ[0][0];
    scalar a1 = deltaQ[1][0];
    scalar a2 = deltaQ[2][0];
    scalar a3 = deltaQ[3][0];
    scalar a4 = deltaQ[4][0];
    myMat& mat0 = tensordPdQ_pos[0];
    myMat& mat1 = tensordPdQ_pos[1];
    myMat& mat2 = tensordPdQ_pos[2];
    myMat& mat3 = tensordPdQ_pos[3];
    myMat& mat4 = tensordPdQ_pos[4];

    tmp = 
	add(
	    outterProd(a0,mat0),
	    add(
		outterProd(a1,mat1),
		add(
		    outterProd(a2,mat2),
		    add(
			outterProd(a3,mat3),
			outterProd(a4,mat4)
		    )
		)
	    )
	);

    mat = innProd(P_invLambda_abs, tmp);

    return mat;
}

myMat tensorInnerProd1
(
    myTensor tensorIn,
    myMat matIn
)
{
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
   
//     mat += matIn[0][0]*tensorIn[0]
    scalar a0 = matIn[0][0];
    scalar a1 = matIn[1][0];
    scalar a2 = matIn[2][0];
    scalar a3 = matIn[3][0];
    scalar a4 = matIn[4][0];
    myMat& mat0 = tensorIn[0];
    myMat& mat1 = tensorIn[1];
    myMat& mat2 = tensorIn[2];
    myMat& mat3 = tensorIn[3];
    myMat& mat4 = tensorIn[4];

    mat = 
	add(
	    outterProd(a0,mat0),
	    add(
		outterProd(a1,mat1),
		add(
		    outterProd(a2,mat2),
		    add(
			outterProd(a3,mat3),
			outterProd(a4,mat4)
		    )
		)
	    )
	);
       	 
    return mat;
}

myMat tensorInnerProd2
(
    myMat P_inv,
    myMat matdLambda_absdQ_pos,
    myMat PdeltaQ
)
{
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );

    myMat tmp
    (
	5,
	std::vector<double>(5,0)
    );
    
    tmp[0] = outterProd(PdeltaQ[0][0], matdLambda_absdQ_pos[0]);
    tmp[1] = outterProd(PdeltaQ[1][0], matdLambda_absdQ_pos[1]);
    tmp[2] = outterProd(PdeltaQ[2][0], matdLambda_absdQ_pos[2]);
    tmp[3] = outterProd(PdeltaQ[3][0], matdLambda_absdQ_pos[3]);
    tmp[4] = outterProd(PdeltaQ[4][0], matdLambda_absdQ_pos[4]);

    mat = innProd(P_inv, tmp);

    return mat;
}

myMat tensorInnerProd3
(
    myMat P_invLambda_abs,
    myTensor tensordPdQ_pos,
    myMat deltaQ
)
{
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    myMat deltaQT
    (
	1,
	std::vector<double>(5,0)
    );

    myMat tmp
    (
	5,
	std::vector<double>(5,0)
    );

    deltaQT[0][0] = deltaQ[0][0];
    deltaQT[0][1] = deltaQ[1][0];
    deltaQT[0][2] = deltaQ[2][0];
    deltaQT[0][3] = deltaQ[3][0];
    deltaQT[0][4] = deltaQ[4][0];

    myMat mat0 = innProd(deltaQT, tensordPdQ_pos[0]);
    myMat mat1 = innProd(deltaQT, tensordPdQ_pos[1]);
    myMat mat2 = innProd(deltaQT, tensordPdQ_pos[2]);
    myMat mat3 = innProd(deltaQT, tensordPdQ_pos[3]);
    myMat mat4 = innProd(deltaQT, tensordPdQ_pos[4]);

    tmp[0] = mat0[0];
    tmp[1] = mat1[0];
    tmp[2] = mat2[0];
    tmp[3] = mat3[0];
    tmp[4] = mat4[0];

    mat = innProd(P_invLambda_abs, tmp);

    return mat;
}

myMat dRoeDampingdQ_pos
(
    const vector sf,
    const scalar rho_tilt,
    const vector U_tilt,
    const scalar H_tilt,
    const scalar rho_pos,
    const scalar rho_neg,
    const vector rhoU_pos,
    const vector rhoU_neg,
    const scalar rhoE_pos,
    const scalar rhoE_neg,
    const scalar gamma
)
{
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    myMat deltaQ
    (
	5,
	std::vector<double>(1,0)
    );
    
    deltaQ[0][0]= rho_neg - rho_pos;
    deltaQ[1][0]= rhoU_neg.x() - rhoU_pos.x();
    deltaQ[2][0]= rhoU_neg.y() - rhoU_pos.y();
    deltaQ[3][0]= rhoU_neg.z() - rhoU_pos.z();
    deltaQ[4][0]= rhoE_neg - rhoE_pos;
   
    myMat P_inv = calcEigVecR(sf, U_tilt, H_tilt, gamma);
    myMat P = calcEigVecL(sf, U_tilt, H_tilt, gamma);
    myMat lambda_abs = calcEigVal_abs(sf, U_tilt, H_tilt, gamma);
    
    const vector U_pos = rhoU_pos/rho_pos;
    const scalar H_pos = gamma*rhoE_pos/rho_pos-.5*(gamma-1)*magSqr(U_pos);
    const vector U_neg = rhoU_neg/rho_neg;
    const scalar H_neg = gamma*rhoE_neg/rho_neg-.5*(gamma-1)*magSqr(U_neg);
    
    //mat1
    myTensor tensordP_invdQ_pos = 
	dP_invdQ_pos 
	(
	    sf,
	    rho_tilt,
	    U_tilt,
	    H_tilt,
	    rho_pos,
	    rho_neg,
	    U_pos,
	    U_neg,
	    H_pos,
	    H_neg,
	    gamma
	);

    myMat lambda_absPdeltaQ = 
	innProd(
	    lambda_abs,
	    innProd(
		P,
		deltaQ
	    )
	);

    myMat mat1(tensorInnerProd1(tensordP_invdQ_pos, lambda_absPdeltaQ));
    //mat2
    myMat PdeltaQ = innProd(P, deltaQ);
    myMat matdLambda_absdQ_pos =
	dLambda_absdQ_pos
	(
	    sf,
	    rho_tilt,
	    U_tilt,
	    H_tilt,
	    rho_pos,
	    rho_neg,
	    U_pos,
	    U_neg,
	    H_pos,
	    H_neg,
	    gamma
	);
    myMat mat2(tensorInnerProd2(P_inv, matdLambda_absdQ_pos, PdeltaQ));
    //mat3 
    myMat P_invLambda_abs = innProd(P_inv, lambda_abs);
    myTensor tensordPdQ_pos = 
	dPdQ_pos 
	(
	    sf,
	    rho_tilt,
	    U_tilt,
	    H_tilt,
	    rho_pos,
	    rho_neg,
	    U_pos,
	    U_neg,
	    H_pos,
	    H_neg,
	    gamma
	);
      
    myMat mat3(tensorInnerProd3(P_invLambda_abs, tensordPdQ_pos, deltaQ));
    //mat4
    myMat mat4 = 
	innProd(
	   P_invLambda_abs,
	   P
	);

    mat =
	add
	(
	    mat1,
	    add
	    (
		mat2,
		add
		(
		    mat3,
		    outterProd(-1,mat4)
		)
	    )
	);
    
    return mat;
//    return outterProd(-1,mat4);
}

myMat dRoeDampingdQ_neg
(
    const vector sf,
    const scalar rho_tilt,
    const vector U_tilt,
    const scalar H_tilt,
    const scalar rho_pos,
    const scalar rho_neg,
    const vector rhoU_pos,
    const vector rhoU_neg,
    const scalar rhoE_pos,
    const scalar rhoE_neg,
    const scalar gamma
)
{

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    myMat deltaQ
    (
	5,
	std::vector<double>(1,0)
    );
    
    deltaQ[0][0]= rho_neg - rho_pos;
    deltaQ[1][0]= rhoU_neg.x() - rhoU_pos.x();
    deltaQ[2][0]= rhoU_neg.y() - rhoU_pos.y();
    deltaQ[3][0]= rhoU_neg.z() - rhoU_pos.z();
    deltaQ[4][0]= rhoE_neg - rhoE_pos;
   
    myMat P_inv = calcEigVecR(sf, U_tilt, H_tilt, gamma);
    myMat P = calcEigVecL(sf, U_tilt, H_tilt, gamma);
    myMat lambda_abs = calcEigVal_abs(sf, U_tilt, H_tilt, gamma);
    
    const vector U_pos = rhoU_pos/rho_pos;
    const scalar H_pos = gamma*rhoE_pos/rho_pos-.5*(gamma-1)*magSqr(U_pos);
    const vector U_neg = rhoU_neg/rho_neg;
    const scalar H_neg = gamma*rhoE_neg/rho_neg-.5*(gamma-1)*magSqr(U_neg);
    
    //mat1
    myTensor tensordP_invdQ_pos = 
	dP_invdQ_pos 
	(
	    sf,
	    rho_tilt,
	    U_tilt,
	    H_tilt,
	    rho_neg,
	    rho_pos,
	    U_neg,
	    U_pos,
	    H_neg,
	    H_pos,
	    gamma
	);

    myMat lambda_absPdeltaQ = 
	innProd(
	    lambda_abs,
	    innProd(
		P,
		deltaQ
	    )
	);

    myMat mat1(tensorInnerProd1(tensordP_invdQ_pos, lambda_absPdeltaQ));
    //mat2
    myMat PdeltaQ = innProd(P, deltaQ);
    myMat matdLambda_absdQ_pos =
	dLambda_absdQ_pos
	(
	    sf,
	    rho_tilt,
	    U_tilt,
	    H_tilt,
	    rho_neg,
	    rho_pos,
	    U_neg,
	    U_pos,
	    H_neg,
	    H_pos,
	    gamma
	);
    myMat mat2(tensorInnerProd2(P_inv, matdLambda_absdQ_pos, PdeltaQ));
    //mat3 
    myMat P_invLambda_abs = innProd(P_inv, lambda_abs);
    myTensor tensordPdQ_pos = 
	dPdQ_pos 
	(
	    sf,
	    rho_tilt,
	    U_tilt,
	    H_tilt,
	    rho_neg,
	    rho_pos,
	    U_neg,
	    U_pos,
	    H_neg,
	    H_pos,
	    gamma
	);
      
    myMat mat3(tensorInnerProd3(P_invLambda_abs, tensordPdQ_pos, deltaQ));
    //mat4
    myMat mat4 = 
	innProd(
	   P_invLambda_abs,
	   P
	);

    mat =
	add
	(
	    mat1,
	    add
	    (
		mat2,
		add
		(
		    mat3,
		    mat4
		)
	    )
	);
    
    
    return mat;
}

