#include <limits>

typedef std::numeric_limits< double > dbl;

typedef std::vector<std::vector<double> > myMat;
template<class myType>
void myPrintMat(std::vector<std::vector<myType> > mat)
{
    int nr = mat.size();
    int nc = mat[0].size();	
    cout << "row #:" << nr << "\n";
    cout << "column #:" << nc << "\n";
//    cout.precision(17);
    for (unsigned int r = 0; r < mat.size(); r++) {
	for (unsigned int c = 0; c < mat[r].size(); c++) {
	    cout << mat[r][c] << " ";
	}
	cout << "\n";
    }
    return;
}

myMat innProd(myMat mat1, myMat mat2)
{
    unsigned int nr1 = mat1.size();
    unsigned int nc1 = mat1[0].size();
    unsigned int nr2 = mat2.size();
    unsigned int nc2 = mat2[0].size();

    myMat mat3
    (
	nr1,
	std::vector<double>(nc2, 0)
    );
    
    if (nc1 == nr2) {
	for (unsigned int r = 0; r < nr1; ++r) {
	    for (unsigned int c = 0; c < nc2; ++c) {
		for (unsigned int i = 0; i < nr2; ++i) {
		    mat3[r][c] += mat1[r][i]*mat2[i][c];
		}
	    }
	}	
    }

    return mat3;
}

void myMatWrite(const myMat& mat, const char* fileName)
{
    std::ofstream out(fileName);
//    dbl::digits10+2
    out.precision(17);
    for (unsigned int r = 0; r < mat.size(); r++) {
	for (unsigned int c = 0; c < mat[r].size(); c++) {
	    out << mat[r][c] << "    ";
	}
	out << "\n";
    }
    return;
}

void myMatWriteDen(const SpMat& mat, const char* fileName)
{
    std::ofstream out(fileName);
//    dbl::digits10+2
    out.precision(17);
    Eigen::MatrixXd dMat(mat);
    out << dMat << std::endl;    
}

void myMatWriteSp(const SpMat& mat, const char* fileName)
{
    std::ofstream out(fileName);
//    dbl::digits10+2
    out.precision(17);
    out << mat << std::endl;
}

void myEigen2Matlab(const SpMat& mat, const char* fileName)
{
    std::ofstream out(fileName);
    out.precision(17);
//    int lastRow = 0;
//    int lastCol = 0;
    for (int k=0; k<mat.outerSize(); ++k) {
	for (SpMat::InnerIterator it(mat,k); it; ++it) {
//	    lastRow = (it.row()+1);
//	    lastCol = (it.col()+1);
	    out << (it.row()+1) << " "
		<< (it.col()+1) << " "
		<< it.value() << " "
		<< "\n";
	    /*it.value();
	    it.row();   // row index
	    it.col();   // col index (here it is equal to k)
	    it.index(); // inner index, here it is equal to it.row()*/
	}
    }

    out << mat.rows() << " "
	<< mat.cols() << " "
	<< 0 << " "
	<< "\n";
}

void myMatWriteDenSp(const SpMat& mat, const char* fileName)
{
//    myMatWriteDen(mat, fileName);
//    myMatWriteSp(mat, fileName);
    myEigen2Matlab(mat, fileName);
}

