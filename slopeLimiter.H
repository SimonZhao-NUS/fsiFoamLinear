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
    A slope limiter

SourceFiles
    slopeLimiter.C

\*---------------------------------------------------------------------------*/

#ifndef slopeLimiter_H
#define slopeLimiter_H

#include "Field.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace szhaoFoam
{

//template<class Type>
//class slopeLimiter;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class slopeLimiter
:
    public Field<Type>
{

public:
    slopeLimiter(const Field<Type>& field);

    slopeLimiter(const label size, const one)
    :
	Field<Type>(size, Zero)
    {
	Field<Type>& f = *this;
	for(label n=0; n<size; n++) {
	    for(direction cmpt=0; cmpt<Type::nComponents; cmpt++) {
		f[n].component(cmpt) = 1;
	    }
	}
    };

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace szhaoFoam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace szhaoFoam
{
//inspired by ~/OpenFOAM/OpenFOAM-v1606+/src/finiteVolume/lnInclude/fvScalarMatrix.H
template<>
slopeLimiter<double>::slopeLimiter(const label size, const one)
:
    Field<double>(size, Zero)
{
	Field<double>& f = *this;
	for(label n=0; n<size; n++) {
	    f[n] = 1;
	}
};

}
#endif

