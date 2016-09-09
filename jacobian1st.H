myMat calcEigVecR
(
    const vector sf,
    const vector U,
    const scalar H,
    const scalar gamma
)
{
    const vector n = sf / mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar U_sqr = magSqr(U);
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const vector e1(1, 0, 0);
    const vector e2(0, 1, 0);
    const vector e3(0, 0, 1);
    const scalar c = std::sqrt((gamma-1)*(H-.5*magSqr(U)));
    const scalar c_sqr = c*c;
//  scalar p = rho*c_sqr/gamma;
/*    sqrt(
          sf.x()*sf.x() +
	    sf.y()*sf.y() +
	    sf.z()*sf.z()
	);
*/
    myMat mat
    (
	5,
	std::vector<double>(5, 0)
    );

    //r1
    mat[0][0] = n1;
    mat[0][1] = n2;
    mat[0][2] = n3;
    mat[0][3] = 1/(2*c_sqr);
    mat[0][4] = 1/(2*c_sqr);
    //r2
    mat[1][0] = u*n1;
    mat[1][1] = u*n2 - n3;
    mat[1][2] = u*n3 + n2;
    mat[1][3] = (u+c*n1)/(2*c_sqr);
    mat[1][4] = (u-c*n1)/(2*c_sqr);
    //r3
    mat[2][0] = v*n1 + n3;
    mat[2][1] = v*n2;
    mat[2][2] = v*n3 - n1;
    mat[2][3] = (v+c*n2)/(2*c_sqr);
    mat[2][4] = (v-c*n2)/(2*c_sqr);
    //r4
    mat[3][0] = w*n1 - n2;
    mat[3][1] = w*n2 + n1;
    mat[3][2] = w*n3;
    mat[3][3] = (w+c*n3)/(2*c_sqr);
    mat[3][4] = (w-c*n3)/(2*c_sqr);
    //r5
    mat[4][0] = (U_sqr*n/2 - (n^U))&e1;
    mat[4][1] = (U_sqr*n/2 - (n^U))&e2;
    mat[4][2] = (U_sqr*n/2 - (n^U))&e3;
    mat[4][3] = (H+c*(U&n))/(2*c_sqr);
    mat[4][4] = (H-c*(U&n))/(2*c_sqr);

    return mat;
}

myMat calcEigVecL
(
    const vector sf,
    const vector U,
    const scalar H,
    const scalar gamma
)
{
    const vector n = sf / mag(sf);
    const scalar n1 = n.x();
    const scalar n2 = n.y();
    const scalar n3 = n.z();
    const scalar U_sqr = magSqr(U);
    const scalar u = U.x();
    const scalar v = U.y();
    const scalar w = U.z();
    const vector e1(1, 0, 0);
    const vector e2(0, 1, 0);
    const vector e3(0, 0, 1);
    const scalar gamma_s = gamma - 1;
    const scalar c = std::sqrt((gamma-1)*(H-.5*magSqr(U)));
    const scalar c_sqr = c*c;
//    scalar p = rho*c_sqr/gamma;
    
    myMat mat
    (
	5,
	std::vector<double>(5, 0)
    );

    //r1
    mat[0][0] = (n-(gamma_s*U_sqr)/(2*c_sqr)*n+(n^U))&e1;
    mat[0][1] = n1*gamma_s*u/(c_sqr);
    mat[0][2] = n3 + n1*gamma_s*v/(c_sqr);
    mat[0][3] = -n2 + n1*gamma_s*w/(c_sqr);
    mat[0][4] = -n1*gamma_s/(c_sqr);
    //r2
    mat[1][0] = (n-(gamma_s*U_sqr)/(2*c_sqr)*n+(n^U))&e2;
    mat[1][1] = -n3 + n2*gamma_s*u/(c_sqr);                  
    mat[1][2] = n2*gamma_s*v/(c_sqr);             
    mat[1][3] = n1 + n2*gamma_s*w/(c_sqr);            
    mat[1][4] = -n2*gamma_s/(c_sqr);                   
    //r3
    mat[2][0] = (n-(gamma_s*U_sqr)/(2*c_sqr)*n+(n^U))&e3;
    mat[2][1] = n2 + n3*gamma_s*u/(c_sqr);            
    mat[2][2] = -n1 + n3*gamma_s*v/(c_sqr);                  
    mat[2][3] = n3*gamma_s*w/(c_sqr);             
    mat[2][4] = -n3*gamma_s/(c_sqr);                   
    //r4
    mat[3][0] = (gamma_s*U_sqr)/2 - c*(U&n);
    mat[3][1] = c*n1 - gamma_s*u;              
    mat[3][2] = c*n2 - gamma_s*v;             
    mat[3][3] = c*n3 - gamma_s*w;                   
    mat[3][4] = gamma_s;                    
    //r5
    mat[4][0] = (gamma_s*U_sqr)/2 + c*(U&n);
    mat[4][1] = -c*n1 - gamma_s*u;         
    mat[4][2] = -c*n2 - gamma_s*v;         
    mat[4][3] = -c*n3 - gamma_s*w;         
    mat[4][4] = gamma_s;                  

    return mat;
}

myMat calcEigVal
(
    const vector sf,
    const vector U,
    const scalar H,
    const scalar gamma
)
{
    std::vector<std::vector<double> > lambda
    (
	5,
	std::vector<double>(5,0)
    );
    const scalar magSf = mag(sf);
    const scalar c = std::sqrt((gamma-1)*(H-magSqr(U)/2));

    lambda[0][0] = sf&U;
    lambda[1][1] = sf&U;
    lambda[2][2] = sf&U;
    lambda[3][3] = (sf&U)+c*magSf;
    lambda[4][4] = (sf&U)-c*magSf;

    return lambda;
}

myMat calcEigVal_abs
(
    const vector sf,
    const vector U,
    const scalar H,
    const scalar gamma
)
{
    std::vector<std::vector<double> > lambda
    (
	5,
	std::vector<double>(5,0)
    );
    const scalar magSf = mag(sf);
    const scalar c = std::sqrt((gamma-1)*(H-magSqr(U)/2));

    lambda[0][0] = std::fabs(sf&U);
    lambda[1][1] = std::fabs(sf&U);
    lambda[2][2] = std::fabs(sf&U);
    lambda[3][3] = std::fabs((sf&U)+c*magSf);
    lambda[4][4] = std::fabs((sf&U)-c*magSf);

    return lambda;
}

