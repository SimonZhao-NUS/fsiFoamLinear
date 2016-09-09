#include <vector>
#include <cmath>
#include <fstream>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Triplet;
typedef std::vector<Triplet> TripletList;

#include "myFuncSp.H"
#include "myInterpolate.H"
#include "linearizedRoe.H"
#include "roeAverage.H"
#include "getMatFuncs.H"
//#include "linearSolvers.H"
#include "farFieldBoundary.C"
