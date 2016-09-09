myMat dR1dM
(
    const vector sf,
    const vector U
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    
    myMat mat
    (
        5,
        std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = 0;
    mat[0][1] = 0;
    mat[0][2] = 0;
    mat[0][3] = 0;
    mat[0][4] = 0;
    //r2
    mat[1][0] = 0; 
    mat[1][1] = n1; 
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = 0; 
    //r3
    mat[2][0] = 0; 
    mat[2][1] = 0; 
    mat[2][2] = n1; 
    mat[2][3] = 0; 
    mat[2][4] = 0; 
    //r4
    mat[3][0] = 0; 
    mat[3][1] = 0; 
    mat[3][2] = 0; 
    mat[3][3] = n1; 
    mat[3][4] = 0; 
    //r5
    mat[4][0] = 0; 
    mat[4][1] = u*n1; 
    mat[4][2] = v*n1+n3; 
    mat[4][3] = w*n1-n2; 
    mat[4][4] = 0; 

    return mat;
}

myMat dR2dM
(
    const vector sf,
    const vector U
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = 0;
    mat[0][1] = 0;
    mat[0][2] = 0;
    mat[0][3] = 0;
    mat[0][4] = 0;
    //r2
    mat[1][0] = 0; 
    mat[1][1] = n2;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = 0; 
    //r3
    mat[2][0] = 0; 
    mat[2][1] = 0; 
    mat[2][2] = n2; 
    mat[2][3] = 0; 
    mat[2][4] = 0; 
    //r4
    mat[3][0] = 0; 
    mat[3][1] = 0; 
    mat[3][2] = 0; 
    mat[3][3] = n2; 
    mat[3][4] = 0; 
    //r5
    mat[4][0] = 0; 
    mat[4][1] = u*n2-n3; 
    mat[4][2] = v*n2; 
    mat[4][3] = w*n2+n1; 
    mat[4][4] = 0; 

    return mat;
}

myMat dR3dM
(
    const vector sf,
    const vector U
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = 0;
    mat[0][1] = 0;
    mat[0][2] = 0;
    mat[0][3] = 0;
    mat[0][4] = 0;
    //r2
    mat[1][0] = 0; 
    mat[1][1] = n3;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = 0; 
    //r3
    mat[2][0] = 0; 
    mat[2][1] = 0; 
    mat[2][2] = n3; 
    mat[2][3] = 0; 
    mat[2][4] = 0; 
    //r4
    mat[3][0] = 0; 
    mat[3][1] = 0; 
    mat[3][2] = 0; 
    mat[3][3] = n3; 
    mat[3][4] = 0; 
    //r5
    mat[4][0] = 0; 
    mat[4][1] = u*n3+n2; 
    mat[4][2] = v*n3-n1; 
    mat[4][3] = w*n3; 
    mat[4][4] = 0; 

    return mat;
}

myMat dR4dM
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar q2 = magSqr(U);
    const scalar c2 = c*c;
    const scalar rhoC = rho*c;
    const scalar rhoC2 = rho*c2;
    const scalar pC = p*c;
    const scalar pC2 = p*c2;
    
    
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = 1/(2*gamma*p);
    mat[0][1] = 0;
    mat[0][2] = 0;
    mat[0][3] = 0;
    mat[0][4] = -rho/(2*gamma*p*p);
    //r2
    mat[1][0] = u/(2*rhoC2) + n1/(4*rhoC); 
    mat[1][1] = 1/(2*c2);
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = -u/(2*pC2)-n1/(4*pC); 
    //r3
    mat[2][0] = v/(2*rhoC2) + n2/(4*rhoC);
    mat[2][1] = 0; 
    mat[2][2] = 1/(2*c2); 
    mat[2][3] = 0; 
    mat[2][4] = -v/(2*pC2)-n2/(4*pC); 
    //r4
    mat[3][0] = w/(2*rhoC2) + n3/(4*rhoC);
    mat[3][1] = 0;                        
    mat[3][2] = 0;                 
    mat[3][3] = 1/(2*c2);
    mat[3][4] = -w/(2*pC2)-n3/(4*pC);     
    //r5
    mat[4][0] = q2/(4*rhoC2) + (n&U)/(4*rhoC);
    mat[4][1] = u/(2*c2) + n1/(2*c);
    mat[4][2] = v/(2*c2) + n2/(2*c);
    mat[4][3] = w/(2*c2) + n3/(2*c);
    mat[4][4] = -q2/(4*pC2) - (n&U)/(4*pC); 

    return mat;
}

myMat dR5dM
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar q2 = magSqr(U);
    const scalar c2 = c*c;
    const scalar rhoC = rho*c;
    const scalar rhoC2 = rho*c2;
    const scalar pC = p*c;
    const scalar pC2 = p*c2;


    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = 1/(2*gamma*p);
    mat[0][1] = 0;
    mat[0][2] = 0;
    mat[0][3] = 0;
    mat[0][4] = -rho/(2*gamma*p*p);
    //r2
    mat[1][0] = u/(2*rhoC2) - n1/(4*rhoC); 
    mat[1][1] = 1/(2*c2);
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = -u/(2*pC2) + n1/(4*pC); 
    //r3
    mat[2][0] = v/(2*rhoC2) - n2/(4*rhoC);
    mat[2][1] = 0; 
    mat[2][2] = 1/(2*c2); 
    mat[2][3] = 0; 
    mat[2][4] = -v/(2*pC2) + n2/(4*pC); 
    //r4
    mat[3][0] = w/(2*rhoC2) - n3/(4*rhoC);
    mat[3][1] = 0;                        
    mat[3][2] = 0;                 
    mat[3][3] = 1/(2*c2);
    mat[3][4] = -w/(2*pC2) + n3/(4*pC);     
    //r5
    mat[4][0] = q2/(4*gamma*p) - (n&U)/(4*rhoC);
    mat[4][1] = u/(2*c2) - n1/(2*c);
    mat[4][2] = v/(2*c2) - n2/(2*c);
    mat[4][3] = w/(2*c2) - n3/(2*c);
    mat[4][4] = -q2/(4*pC2) + (n&U)/(4*pC); 

    return mat;
}

myMat dL1dM
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar gamma_s = gamma -1;
    const scalar phi2 = .5*gamma_s*magSqr(U);
    const scalar c2 = c*c;
    const scalar rhoC2 = rho*c2;
    const scalar pC2 = p*c2;
     
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = -n1*phi2/(rhoC2);
    mat[0][1] = -n1*gamma_s*u/c2;
    mat[0][2] = -n1*gamma_s*v/c2 - n3;
    mat[0][3] = -n1*gamma_s*w/c2 + n2;
    mat[0][4] =  n1*phi2/(pC2);       
    //r2
    mat[1][0] = n1*gamma_s*u/rhoC2; 
    mat[1][1] = n1*gamma_s/c2;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = -n1*u*gamma_s/pC2; 
    //r3
    mat[2][0] = n1*gamma_s*v/rhoC2;
    mat[2][1] = 0;     
    mat[2][2] = n1*gamma_s/c2;
    mat[2][3] = 0;                 
    mat[2][4] = -n1*v*gamma_s/pC2; 
    //r4
    mat[3][0] = n1*gamma_s*w/rhoC2;
    mat[3][1] = 0;                 
    mat[3][2] = 0;     
    mat[3][3] = n1*gamma_s/c2; 
    mat[3][4] = -n1*w*gamma_s/pC2; 
    //r5
    mat[4][0] = -n1*gamma_s/rhoC2;
    mat[4][1] = 0;                 
    mat[4][2] = 0;                 
    mat[4][3] = 0;     
    mat[4][4] = n1*gamma_s/pC2; 

    return mat;
}

myMat dL2dM
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar gamma_s = gamma -1;
    const scalar phi2 = .5*gamma_s*magSqr(U);
    const scalar c2 = c*c;
    const scalar rhoC2 = rho*c2;
    const scalar pC2 = p*c2;

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = -n2*phi2/(rhoC2);
    mat[0][1] = -n2*gamma_s*u/c2 + n3;
    mat[0][2] = -n2*gamma_s*v/c2;
    mat[0][3] = -n2*gamma_s*w/c2 - n1;
    mat[0][4] =  n2*phi2/(pC2);       
    //r2
    mat[1][0] = n2*gamma_s*u/rhoC2; 
    mat[1][1] = n2*gamma_s/c2;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = -n2*u*gamma_s/pC2; 
    //r3
    mat[2][0] = n2*gamma_s*v/rhoC2;
    mat[2][1] = 0;     
    mat[2][2] = n2*gamma_s/c2;
    mat[2][3] = 0;                 
    mat[2][4] = -n2*v*gamma_s/pC2; 
    //r4
    mat[3][0] = n2*gamma_s*w/rhoC2;
    mat[3][1] = 0;                 
    mat[3][2] = 0;     
    mat[3][3] = n2*gamma_s/c2; 
    mat[3][4] = -n2*w*gamma_s/pC2; 
    //r5
    mat[4][0] = -n2*gamma_s/rhoC2;
    mat[4][1] = 0;                 
    mat[4][2] = 0;                 
    mat[4][3] = 0;     
    mat[4][4] = n2*gamma_s/pC2; 

    return mat;
}

myMat dL3dM
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar gamma_s = gamma -1;
    const scalar phi2 = .5*gamma_s*magSqr(U);
    const scalar c2 = c*c;
    const scalar rhoC2 = rho*c2;
    const scalar pC2 = p*c2;

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = -n3*phi2/(rhoC2);
    mat[0][1] = -n3*gamma_s*u/c2 - n2;
    mat[0][2] = -n3*gamma_s*v/c2 + n1;
    mat[0][3] = -n3*gamma_s*w/c2;
    mat[0][4] =  n3*phi2/(pC2);       
    //r2
    mat[1][0] = n3*gamma_s*u/rhoC2; 
    mat[1][1] = n3*gamma_s/c2;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = -n3*u*gamma_s/pC2; 
    //r3
    mat[2][0] = n3*gamma_s*v/rhoC2;
    mat[2][1] = 0;     
    mat[2][2] = n3*gamma_s/c2;
    mat[2][3] = 0;                 
    mat[2][4] = -n3*v*gamma_s/pC2; 
    //r4
    mat[3][0] = n3*gamma_s*w/rhoC2;
    mat[3][1] = 0;                 
    mat[3][2] = 0;     
    mat[3][3] = n3*gamma_s/c2; 
    mat[3][4] = -n3*w*gamma_s/pC2; 
    //r5
    mat[4][0] = -n3*gamma_s/rhoC2;
    mat[4][1] = 0;                 
    mat[4][2] = 0;                 
    mat[4][3] = 0;     
    mat[4][4] = n3*gamma_s/pC2; 

    return mat;
}

myMat dL4dM
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar gamma_s = gamma -1;
    
    myMat mat
    (
        5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = c*(U&n)/(2*rho);
    mat[0][1] = gamma_s*u - c*n1;
    mat[0][2] = gamma_s*v - c*n2;
    mat[0][3] = gamma_s*w - c*n3;
    mat[0][4] = -c*(U&n)/(2*p);  
    //r2
    mat[1][0] = -c*n1/(2*rho);
    mat[1][1] = -gamma_s;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = c*n1/(2*p); 
    //r3
    mat[2][0] = -c*n2/(2*rho);
    mat[2][1] = 0;    
    mat[2][2] = -gamma_s;           
    mat[2][3] = 0;           
    mat[2][4] = c*n2/(2*p);  
    //r4
    mat[3][0] = -c*n3/(2*rho);
    mat[3][1] = 0;            
    mat[3][2] = 0;
    mat[3][3] = -gamma_s;     
    mat[3][4] = c*n3/(2*p);   
    //r5
    mat[4][0] = 0;
    mat[4][1] = 0;                 
    mat[4][2] = 0;                 
    mat[4][3] = 0;     
    mat[4][4] = 0; 

    return mat;
}

myMat dL5dM
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar gamma_s = gamma -1;

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = -c*(U&n)/(2*rho);
    mat[0][1] = gamma_s*u + c*n1;
    mat[0][2] = gamma_s*v + c*n2;
    mat[0][3] = gamma_s*w + c*n3;
    mat[0][4] = c*(U&n)/(2*p);   
    //r2
    mat[1][0] = c*n1/(2*rho);
    mat[1][1] = -gamma_s;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = -c*n1/(2*p); 
    //r3
    mat[2][0] = c*n2/(2*rho);
    mat[2][1] = 0;    
    mat[2][2] = -gamma_s;           
    mat[2][3] = 0;           
    mat[2][4] = -c*n2/(2*p);  
    //r4
    mat[3][0] = c*n3/(2*rho);
    mat[3][1] = 0;            
    mat[3][2] = 0;
    mat[3][3] = -gamma_s;     
    mat[3][4] = -c*n3/(2*p);   
    //r5
    mat[4][0] = 0;
    mat[4][1] = 0;                 
    mat[4][2] = 0;                 
    mat[4][3] = 0;     
    mat[4][4] = 0; 

    return mat;
}


myMat dL1dMT
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar gamma_s = gamma -1;
    const scalar phi2 = .5*gamma_s*magSqr(U);
    const scalar c2 = c*c;
    const scalar rhoC2 = rho*c2;
    const scalar pC2 = p*c2;
     
    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = -n1*phi2/(rhoC2); 
    mat[0][1] = -n1*gamma_s*u/c2;
    mat[0][2] = -n1*gamma_s*v/c2 - n3;
    mat[0][3] = -n1*gamma_s*w/c2 + n2;
    mat[0][4] =  n1*phi2/(pC2);       
    //r2
    mat[1][0] = -n2*phi2/(rhoC2); 
    mat[1][1] = -n2*gamma_s*u/c2 + n3;
    mat[1][2] = -n2*gamma_s*v/c2;
    mat[1][3] = -n2*gamma_s*w/c2 - n1;
    mat[1][4] =  n2*phi2/(pC2);       
    //r3
    mat[2][0] = -n3*phi2/(rhoC2);
    mat[2][1] = -n3*gamma_s*u/c2 - n2;
    mat[2][2] = -n3*gamma_s*v/c2 + n1;
    mat[2][3] = -n3*gamma_s*w/c2;
    mat[2][4] =  n3*phi2/(pC2); 
    //r4
    mat[3][0] = c*(U&n)/(2*rho); 
    mat[3][1] = gamma_s*u - c*n1;
    mat[3][2] = gamma_s*v - c*n2;
    mat[3][3] = gamma_s*w - c*n3;
    mat[3][4] = -c*(U&n)/(2*p);  
    //r5
    mat[4][0] = -c*(U&n)/(2*rho); 
    mat[4][1] = gamma_s*u + c*n1;
    mat[4][2] = gamma_s*v + c*n2;
    mat[4][3] = gamma_s*w + c*n3;
    mat[4][4] = c*(U&n)/(2*p);   

    return mat;
}

myMat dL2dMT
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar gamma_s = gamma -1;
    const scalar phi2 = .5*gamma_s*magSqr(U);
    const scalar c2 = c*c;
    const scalar rhoC2 = rho*c2;
    const scalar pC2 = p*c2;

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = n1*gamma_s*u/rhoC2;  
    mat[0][1] = n1*gamma_s/c2;
    mat[0][2] = 0; 
    mat[0][3] = 0; 
    mat[0][4] = -n1*u*gamma_s/pC2; 
    //r2
    mat[1][0] = n2*gamma_s*u/rhoC2;  
    mat[1][1] = n2*gamma_s/c2;
    mat[1][2] = 0; 
    mat[1][3] = 0; 
    mat[1][4] = -n2*u*gamma_s/pC2; 
    //r3
    mat[2][0] = n3*gamma_s*u/rhoC2;  
    mat[2][1] = n3*gamma_s/c2;
    mat[2][2] = 0; 
    mat[2][3] = 0; 
    mat[2][4] = -n3*u*gamma_s/pC2; 
    //r4
    mat[3][0] = -c*n1/(2*rho); 
    mat[3][1] = -gamma_s;
    mat[3][2] = 0; 
    mat[3][3] = 0; 
    mat[3][4] = c*n1/(2*p); 
    //r5
    mat[4][0] = c*n1/(2*rho); 
    mat[4][1] = -gamma_s;
    mat[4][2] = 0; 
    mat[4][3] = 0; 
    mat[4][4] = -c*n1/(2*p); 

    return mat;
}

myMat dL3dMT
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar gamma_s = gamma -1;
    const scalar phi2 = .5*gamma_s*magSqr(U);
    const scalar c2 = c*c;
    const scalar rhoC2 = rho*c2;
    const scalar pC2 = p*c2;

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = n1*gamma_s*v/rhoC2; 
    mat[0][1] = 0;     
    mat[0][2] = n1*gamma_s/c2;
    mat[0][3] = 0;                 
    mat[0][4] = -n1*v*gamma_s/pC2; 
    //r2
    mat[1][0] = n2*gamma_s*v/rhoC2; 
    mat[1][1] = 0;     
    mat[1][2] = n2*gamma_s/c2;
    mat[1][3] = 0;                 
    mat[1][4] = -n2*v*gamma_s/pC2; 
    //r3
    mat[2][0] = n3*gamma_s*v/rhoC2; 
    mat[2][1] = 0;     
    mat[2][2] = n3*gamma_s/c2;
    mat[2][3] = 0;                 
    mat[2][4] = -n3*v*gamma_s/pC2; 
    //r4
    mat[3][0] = -c*n2/(2*rho); 
    mat[3][1] = 0;    
    mat[3][2] = -gamma_s;     
    mat[3][3] = 0;           
    mat[3][4] = c*n2/(2*p);  
    //r5
    mat[4][0] = c*n2/(2*rho); 
    mat[4][1] = 0;    
    mat[4][2] = -gamma_s;    
    mat[4][3] = 0;           
    mat[4][4] = -c*n2/(2*p); 

    return mat;
}

myMat dL4dMT
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar gamma_s = gamma -1;
    const scalar c2 = c*c;
    const scalar rhoC2 = rho*c2;
    const scalar pC2 = p*c2;

    myMat mat
    (
        5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = n1*gamma_s*w/rhoC2; 
    mat[0][1] = 0;                 
    mat[0][2] = 0;     
    mat[0][3] = n1*gamma_s/c2; 
    mat[0][4] = -n1*w*gamma_s/pC2; 
    //r2
    mat[1][0] = n2*gamma_s*w/rhoC2; 
    mat[1][1] = 0;                 
    mat[1][2] = 0;     
    mat[1][3] = n2*gamma_s/c2; 
    mat[1][4] = -n2*w*gamma_s/pC2; 
    //r3
    mat[2][0] = n3*gamma_s*w/rhoC2; 
    mat[2][1] = 0;                 
    mat[2][2] = 0;     
    mat[2][3] = n3*gamma_s/c2; 
    mat[2][4] = -n3*w*gamma_s/pC2; 
    //r4
    mat[3][0] = -c*n3/(2*rho); 
    mat[3][1] = 0;            
    mat[3][2] = 0;
    mat[3][3] = -gamma_s;     
    mat[3][4] = c*n3/(2*p);   
    //r5
    mat[4][0] = c*n3/(2*rho); 
    mat[4][1] = 0;            
    mat[4][2] = 0;
    mat[4][3] = -gamma_s;     
    mat[4][4] = -c*n3/(2*p);   

    return mat;
}

myMat dL5dMT
(
    const vector sf,
    const scalar rho,
    const vector U,
    const scalar p,
    const scalar gamma
)
{
    const vector n = sf/mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const scalar c = std::sqrt(gamma*p/rho);
     
    const scalar gamma_s = gamma -1;
    const scalar c2 = c*c;
    const scalar rhoC2 = rho*c2;
    const scalar pC2 = p*c2;

    myMat mat
    (
	5,
	std::vector<double>(5,0)
    );
    
    //r1
    mat[0][0] = -n1*gamma_s/rhoC2; 
    mat[0][1] = 0;                 
    mat[0][2] = 0;                 
    mat[0][3] = 0;     
    mat[0][4] = n1*gamma_s/pC2; 
    //r2
    mat[1][0] = -n2*gamma_s/rhoC2; 
    mat[1][1] = 0;                 
    mat[1][2] = 0;                 
    mat[1][3] = 0;     
    mat[1][4] = n2*gamma_s/pC2; 
    //r3
    mat[2][0] = -n3*gamma_s/rhoC2; 
    mat[2][1] = 0;                 
    mat[2][2] = 0;                 
    mat[2][3] = 0;     
    mat[2][4] = n3*gamma_s/pC2; 
    //r4
    mat[3][0] = 0; 
    mat[3][1] = 0;
    mat[3][2] = 0;
    mat[3][3] = 0;
    mat[3][4] = 0;
    //r5
    mat[4][0] = 0; 
    mat[4][1] = 0;
    mat[4][2] = 0;
    mat[4][3] = 0;
    mat[4][4] = 0;
    
    return mat;
}
