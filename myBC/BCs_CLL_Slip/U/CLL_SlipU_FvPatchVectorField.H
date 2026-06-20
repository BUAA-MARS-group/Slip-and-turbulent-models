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
    Foam::CLL_SlipU_FvPatchVectorField

Description
    CLL_ slip boundary condition including thermal creep and surface
    curvature terms that can be optionally switched off.

SourceFiles
    fixedRhoFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef CLL_SlipU_FvPatchVectorField_H
#define CLL_SlipU_FvPatchVectorField_H

#include "mixedFixedValueSlipFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CLL_SlipUFvPatch Declaration
\*---------------------------------------------------------------------------*/

class CLL_SlipU_FvPatchVectorField
:
    public mixedFixedValueSlipFvPatchVectorField
{
    // Private Data

        //- Temperature field name, default = "T"
        word TName_;

        //- Density field name, default = "rho"
        word rhoName_;

        //- Compressibility field name, default = "thermo:psi"
        word psiName_;

        //- Dynamic viscosity field name, default = "thermo:mu"
        word muName_;

        //- tauMC field name, default = "tauMC"
        word tauMCName_;

        // Tangential (Momentum) Accommodation coefficient
        scalar accommodationCoeff_t_;

        // Normal (Energy) Accommodation coefficient
        scalar accommodationCoeff_n_;

        // Wall velocity
        vectorField Uwall_;

        // Include thermal creep term (default on)
        Switch thermalCreep_;

        // Include boundary curvature term (default on)
        Switch curvature_;


public:

    //- Runtime type information
    TypeName("CLL_SlipU");


    // Constructors

        //- Construct from patch and internal field
        CLL_SlipU_FvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        CLL_SlipU_FvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  CLL_SlipU_FvPatchVectorField onto a new patch
        CLL_SlipU_FvPatchVectorField
        (
            const CLL_SlipU_FvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        CLL_SlipU_FvPatchVectorField
        (
            const CLL_SlipU_FvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new CLL_SlipU_FvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchVectorField&, const labelList&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchVectorField&);


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
