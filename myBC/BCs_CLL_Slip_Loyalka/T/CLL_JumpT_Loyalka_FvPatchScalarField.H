/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::CLL_JumpT_Loyalka_FvPatchScalarField

Description
    CLL_ temperature jump boundary condition

SourceFiles
    CLL_JumpT_Loyalka_FvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef CLL_JumpT_Loyalka_FvPatchScalarField_H
#define CLL_JumpT_Loyalka_FvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class CLL_JumpT_Loyalka_FvPatch Declaration
\*---------------------------------------------------------------------------*/

class CLL_JumpT_Loyalka_FvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private Data

        //- Velocity field name, default = "U"
        word UName_;

        //- Density field name, default = "rho"
        word rhoName_;

        //- Compressibility field name, default = "thermo:psi"
        word psiName_;

        //- Dynamic viscosity field name, default = "thermo:mu"
        word muName_;

        //- Tangential (Momentum) Accommodation coefficient
        scalar accommodationCoeff_t_;

        //- Normal (Energy) Accommodation coefficient
        scalar accommodationCoeff_n_;

        //- Wall surface temperature
        scalarField Twall_;

        //- Heat capacity ratio (default 1.4)
        scalar gamma_;

public:

    //- Runtime type information
    TypeName("CLL_JumpT_Loyalka");


    // Constructors

        //- Construct from patch and internal field
        CLL_JumpT_Loyalka_FvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        CLL_JumpT_Loyalka_FvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given CLL_JumpT_Loyalka_FvPatchScalarField
        //  onto a new patch
        CLL_JumpT_Loyalka_FvPatchScalarField
        (
            const CLL_JumpT_Loyalka_FvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        CLL_JumpT_Loyalka_FvPatchScalarField
        (
            const CLL_JumpT_Loyalka_FvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new CLL_JumpT_Loyalka_FvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
