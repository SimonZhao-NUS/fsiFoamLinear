/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    szhaoFoam::limiter

Group
    NA

Description
    Venkatakrishnan's slope limiter

SourceFiles
    vkSlopeLimiter.C

\*---------------------------------------------------------------------------*/

#ifndef vkSlopeLimiter_H
#define vkSlopeLimiter_H

#include "slopeLimiter.H"

namespace szhaoFoam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class Type>
class vkSlopeLimiter
:
    public slopeLimiter<Type>
{
//private properties
private:
    Field<Type> Umax_;
    Field<Type> Umin_;
    Field<Type> delta1max_;
    Field<Type> delta1min_;
    Field<Type> delta2P_;
    Field<Type> delta2N_;
    double K_;

private:
    void findMax(const GeometricField<Type, fvPatchField, volMesh>& U);
    
    void findMin(const GeometricField<Type, fvPatchField, volMesh>& U);
    
    void setDelta
    (
	const GeometricField<Type, fvPatchField, volMesh>& U,
	const GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>& Ugrad
    );

    void setValue
    (
	const GeometricField<Type, fvPatchField, volMesh>& U,
	const GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>& Ugrad
    );

    double func
    (
	const double delta1max,
	const double delta1min,
	const double delta2,
	double epsilonSqr
    ) const;

public:
    //constructor
    vkSlopeLimiter
    (
	const GeometricField<Type, fvPatchField, volMesh>& U,
	const GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>& Ugrad,
	const double K
    )
    :
	slopeLimiter<Type>(U.size(), Foam::one()),
	Umax_(U.size()),
	Umin_(U.size()),
	delta1max_(U.size()),
	delta1min_(U.size()),
	delta2P_(U.mesh().owner().size()),
	delta2N_(U.mesh().neighbour().size()),
	K_(K)
    {
	setValue(U, Ugrad);
    };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace szhaoFoam

//realization
namespace szhaoFoam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void szhaoFoam::vkSlopeLimiter<Type>::findMax(const GeometricField<Type, fvPatchField, volMesh>& U)
{
    const fvMesh& mesh = U.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    
    Umax_ = U;
    forAll(P, fi) {
	for(direction cmpt=0; cmpt<Type::nComponents; cmpt++) {
	    scalar tmp = U[N[fi]].component(cmpt);
	    if (Umax_[P[fi]].component(cmpt) < tmp) {
		Umax_[P[fi]].component(cmpt) = tmp;
	    }

	    tmp = U[P[fi]].component(cmpt);
	    if (Umax_[N[fi]].component(cmpt) < tmp) {
		Umax_[N[fi]].component(cmpt) = tmp;
	    }
	}
    }

    forAll(mesh.boundary(), patchi)
    {
	const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
	const fvPatchField<Type>& Ubp = U.boundaryField()[patchi];

	forAll(mesh.boundary()[patchi], facei)
	{
	    for(direction cmpt=0; cmpt<Type::nComponents; cmpt++) {
		double tmp = Ubp[facei].component(cmpt);
		if (Umax_[pFaceCells[facei]].component(cmpt) < tmp) {
		    Umax_[pFaceCells[facei]].component(cmpt) = tmp;
		}
	    }
	}
    }
    
};

template<>
void szhaoFoam::vkSlopeLimiter<double>::findMax(const GeometricField<double, fvPatchField, volMesh>& U)
{
    const fvMesh& mesh = U.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    
    Umax_ = U;
    forAll(P, fi) {
	scalar tmp = U[N[fi]];
	if (Umax_[P[fi]] < tmp) {
	    Umax_[P[fi]] = tmp;
	}

	tmp = U[P[fi]];
	if (Umax_[N[fi]] < tmp) {
	    Umax_[N[fi]] = tmp;
	}
    }
    
    forAll(mesh.boundary(), patchi)
    {
	const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
	const fvPatchField<double>& Ubp = U.boundaryField()[patchi];

	forAll(mesh.boundary()[patchi], facei)
	{
	    double tmp = Ubp[facei];
	    if (Umax_[pFaceCells[facei]] < tmp) {
		Umax_[pFaceCells[facei]] = tmp;
	    }
	}
    }
};

template<class Type>
void szhaoFoam::vkSlopeLimiter<Type>::findMin(const GeometricField<Type, fvPatchField, volMesh>& U)
{
    const fvMesh& mesh = U.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    
    Umin_ = U;
    forAll(P, fi) {
	for(direction cmpt=0; cmpt<Type::nComponents; cmpt++) {
	    scalar tmp = U[N[fi]].component(cmpt);
	    if (Umin_[P[fi]].component(cmpt) > tmp) {
		Umin_[P[fi]].component(cmpt) = tmp;
	    }

	    tmp = U[P[fi]].component(cmpt);
	    if (Umin_[N[fi]].component(cmpt) > tmp) {
		Umin_[N[fi]].component(cmpt) = tmp;
	    }
	}
    }
    
    forAll(mesh.boundary(), patchi)
    {
	const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
	const fvPatchField<Type>& Ubp = U.boundaryField()[patchi];

	forAll(mesh.boundary()[patchi], facei)
	{
	    for(direction cmpt=0; cmpt<Type::nComponents; cmpt++) {
		double tmp = Ubp[facei].component(cmpt);
		if (Umin_[pFaceCells[facei]].component(cmpt) > tmp) {
		    Umin_[pFaceCells[facei]].component(cmpt) = tmp;
		}
	    }
	}
    }
};

template<>
void szhaoFoam::vkSlopeLimiter<double>::findMin(const GeometricField<double, fvPatchField, volMesh>& U)
{
    const fvMesh& mesh = U.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    
    Umin_ = U;
    forAll(P, fi) {
	scalar tmp = U[N[fi]];
	if (Umin_[P[fi]] > tmp) {
	    Umin_[P[fi]] = tmp;
	}

	tmp = U[P[fi]];
	if (Umin_[N[fi]] > tmp) {
	    Umin_[N[fi]] = tmp;
	}
    }
    
    forAll(mesh.boundary(), patchi)
    {
	const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
	const fvPatchField<double>& Ubp = U.boundaryField()[patchi];

	forAll(mesh.boundary()[patchi], facei)
	{
	    double tmp = Ubp[facei];
	    if (Umin_[pFaceCells[facei]] > tmp) {
		Umin_[pFaceCells[facei]] = tmp;
	    }
	}
    }
};

template<class Type>
void szhaoFoam::vkSlopeLimiter<Type>::setDelta
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>& Ugrad
)
{
///////////find Umax and Umin///////////
    findMax(U);
    findMin(U);
///////////set delta1//////////////
    delta1max_ = Umax_ - U.internalField();
    delta1min_ = Umin_ - U.internalField();
///////////set delta2//////////////
    const surfaceVectorField &faceCenter = U.mesh().Cf();
    const volVectorField &cellCenter = U.mesh().C();
    const labelUList& P = U.mesh().owner();
    const labelUList& N = U.mesh().neighbour();

    forAll(P, fi) {
	vector drt = faceCenter[fi] - cellCenter[P[fi]];
	delta2P_[fi] = drt&Ugrad.internalField()[P[fi]];
/*	tensor grad = Ugrad.internalField()[P[fi]];
    
	delta2P_[fi].x() = drt.x()*grad.xx() + drt.y()*grad.yx();
	delta2P_[fi].y() = drt.x()*grad.xy() + drt.y()*grad.yy();
	delta2P_[fi].z() = 0.0;*/
    }
    
    forAll(N, fi) {
	vector drt = faceCenter[fi] - cellCenter[N[fi]];
	delta2N_[fi] = drt&Ugrad.internalField()[N[fi]];
/*	tensor grad = Ugrad.internalField()[N[fi]];
    
	delta2N_[fi].x() = drt.x()*grad.xx() + drt.y()*grad.yx();
	delta2N_[fi].y() = drt.x()*grad.xy() + drt.y()*grad.yy();
	delta2N_[fi].z() = 0.0;*/
    }
};

template<>
void szhaoFoam::vkSlopeLimiter<double>::setDelta
(
    const GeometricField<double, fvPatchField, volMesh>& U,
    const GeometricField<typename outerProduct<vector, double>::type, fvPatchField, volMesh>& Ugrad
)
{
///////////find Umax and Umin///////////
    findMax(U);
    findMin(U);
///////////set delta1//////////////
    delta1max_ = Umax_ - U.internalField();
    delta1min_ = Umin_ - U.internalField();
///////////set delta2//////////////
    const surfaceVectorField &faceCenter = U.mesh().Cf();
    const volVectorField &cellCenter = U.mesh().C();
    const labelUList& P = U.mesh().owner();
    const labelUList& N = U.mesh().neighbour();

    forAll(P, fi) {
	vector drt = faceCenter[fi] - cellCenter[P[fi]];
	delta2P_[fi] = drt&Ugrad.internalField()[P[fi]];
//	vector grad = Ugrad.internalField()[P[fi]];
//	delta2P_[fi] = drt.x()*grad.x() + drt.y()*grad.y();
    }
    
    forAll(N, fi) {
	vector drt = faceCenter[fi] - cellCenter[N[fi]];
	delta2N_[fi] = drt&Ugrad.internalField()[N[fi]];
//	vector grad = Ugrad.internalField()[N[fi]];
//	delta2N_[fi] = drt.x()*grad.x() + drt.y()*grad.y();
    }
};

template<class Type>
void szhaoFoam::vkSlopeLimiter<Type>::setValue
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>& Ugrad
)
{
    ///////////set Delta values (ref. Computational Fluid Dynamics - Principles and Applications - J. Blazek P166-167)/////////////
    setDelta(U, Ugrad);
    
    const fvMesh& mesh = U.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    const scalarField& vol = mesh.V();
    
    ///////////find min j/////////////
    ///////////func(delta1max, delta1min, delta2) is the segment function///////////
    Field<Type>&f = *this;
    forAll(P, fi) {
	for(direction cmpt=0; cmpt<Type::nComponents; cmpt++) {
	    double psi;
	    psi = 
		func
		(
		    delta1max_[P[fi]].component(cmpt),
		    delta1min_[P[fi]].component(cmpt),
		    delta2P_[fi].component(cmpt),
		    pow(K_*vol[P[fi]],3)
		);
	    if (psi<f[P[fi]].component(cmpt)) {
		f[P[fi]].component(cmpt) = psi;
	    }
	}
    }
    
    forAll(N, fi) {
	for(direction cmpt=0; cmpt<Type::nComponents; cmpt++) {
	    double psi;
	    psi = 
		func
		(
		    delta1max_[N[fi]].component(cmpt),
		    delta1min_[N[fi]].component(cmpt),
		    delta2N_[fi].component(cmpt),
		    pow(K_*vol[N[fi]],3)
		);
	    if (psi<f[N[fi]].component(cmpt)) {
		f[N[fi]].component(cmpt) = psi;
	    }
	}
    }

    const vectorField& pCellCenter = mesh.C().internalField();
    forAll(mesh.boundary(), patchi)
    {
	const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
	const vectorField& pFaceCenter = mesh.Cf().boundaryField()[patchi];
	
	forAll(mesh.boundary()[patchi], facei)
	{
//	    Type delta2B = (pFaceCenter[facei]-pCellCenter[pFaceCells[facei]]) & Ugrad.internalField()[pFaceCells[facei]];

	    vector drt = (pFaceCenter[facei]-pCellCenter[pFaceCells[facei]]);
	    tensor grad = (Ugrad.internalField()[pFaceCells[facei]]);
/*	
	    vector delta2B
	    (
		drt.x()*grad.xx() + drt.y()*grad.yx(),
		drt.x()*grad.xy() + drt.y()*grad.yy(),
		0.0
	    );
*/
	    vector delta2B = drt & grad;

	    for(direction cmpt=0; cmpt<Type::nComponents; cmpt++) {
		double psi = 
		    func
		    (
			delta1max_[pFaceCells[facei]].component(cmpt),
			delta1min_[pFaceCells[facei]].component(cmpt),
			delta2B.component(cmpt),
			pow(K_*vol[pFaceCells[facei]],3)
		    );
		if (psi<f[pFaceCells[facei]].component(cmpt)) {
		    f[pFaceCells[facei]].component(cmpt) = psi;
		}
	    }
	}
    }
};

template<>
void szhaoFoam::vkSlopeLimiter<double>::setValue
(
    const GeometricField<double, fvPatchField, volMesh>& U,
    const GeometricField<vector, fvPatchField, volMesh>& Ugrad
)
{
    ///////////set Delta values (ref. Computational Fluid Dynamics - Principles and Applications - J. Blazek P166-167)/////////////
    setDelta(U, Ugrad);

    const fvMesh& mesh = U.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();
    const scalarField& vol = mesh.V();
    
    ///////////find min j/////////////
    ///////////func(delta1max, delta1min, delta2) is the segment function///////////
    Field<double>&f = *this;
    forAll(P, fi) {
	double psi;
	psi = func(delta1max_[P[fi]], delta1min_[P[fi]], delta2P_[fi], pow(K_*vol[P[fi]],3));
	if (psi<f[P[fi]]) {
	    f[P[fi]] = psi;
	}
    }
    
    forAll(N, fi) {
	double psi;
	psi = func(delta1max_[N[fi]], delta1min_[N[fi]], delta2N_[fi], pow(K_*vol[N[fi]],3));
	if (psi<f[N[fi]]) {
	    f[N[fi]] = psi;
	}
    }
    
    const vectorField& pCellCenter = mesh.C().internalField();
    forAll(mesh.boundary(), patchi)
    {
	const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
	const vectorField& pFaceCenter = mesh.Cf().boundaryField()[patchi];
	
	forAll(mesh.boundary()[patchi], facei)
	{
//	    double delta2B = (pFaceCenter[facei]-pCellCenter[pFaceCells[facei]]) & (Ugrad.internalField()[pFaceCells[facei]]);

	    //when doing inner product operation "&", the third dimenstion is participating in it. However, for a 2D case, this make induce error, in that the mesh is not perfectly 3D. Sf = (x,y, epsilon)
	    vector drt = (pFaceCenter[facei])-(pCellCenter[pFaceCells[facei]]);
	    vector grad = (Ugrad.internalField()[pFaceCells[facei]]);
//	    double delta2B = drt.x()*grad.x() + drt.y()*grad.y();
	    double delta2B = drt & grad;

	    double psi = 
		func
		(
		    delta1max_[pFaceCells[facei]],
		    delta1min_[pFaceCells[facei]],
		    delta2B,
		    pow(K_*vol[pFaceCells[facei]],3)
		);

	    if (psi<f[pFaceCells[facei]]) {
		f[pFaceCells[facei]] = psi;
	    }
	}
    }
/*	    typename vector::labelType validComponents
	    (
		mesh.template validComponents<vector>()
	    );

	    Info << validComponents << endl;
	    cin.get();*/
};

template<class Type>
double szhaoFoam::vkSlopeLimiter<Type>::func
(
    const double delta1max,
    const double delta1min,
    const double delta2,
    double epsilonSqr
) const
{
    ///////////func(delta1max, delta1min, delta2) is the segment function///////////
    double psi = 1;
//    double epsilonSqr = Foam::doubleScalarSMALL;
    double delta2Sqr = delta2*delta2;
    double delta1maxSqr = delta1max*delta1max;
    double delta1minSqr = delta1min*delta1min;

/*    if (delta2>Foam::doubleScalarSMALL) {
	psi = 
	    (
		((delta1maxSqr+epsilonSqr)+2*delta2*delta1max)
		/(delta1maxSqr+2*delta2Sqr+delta1max*delta2+epsilonSqr)
	    );
    } else if (delta2<Foam::doubleScalarSMALL) {
	psi = 
	    (
		((delta1minSqr+epsilonSqr)+2*delta2*delta1min)
		/(delta1minSqr+2*delta2Sqr+delta1min*delta2+epsilonSqr)
	    );
    } else {
	psi = 1;
    }*/

    if (delta2>0) {
	psi = min(1, (delta1max+Foam::doubleScalarSMALL)/(delta2+Foam::doubleScalarSMALL));
    } else if (delta2<0) {
	psi = min(1, (delta1min-Foam::doubleScalarSMALL)/(delta2-Foam::doubleScalarSMALL));
    } else {
	psi = 1;
    }
    
    if (psi < 0) {
	Info << psi << endl;
	Info << delta1max << endl;
	Info << delta1min << endl;
	Info << delta2 << endl;
	Info << epsilonSqr << endl;
	Info << "-----------" << endl;
	cin.get();
    }

    return psi;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace szhaoFoam
#endif
