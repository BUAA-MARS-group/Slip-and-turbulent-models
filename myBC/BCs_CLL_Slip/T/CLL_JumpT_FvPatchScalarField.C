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

\*---------------------------------------------------------------------------*/

#include "CLL_JumpT_FvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CLL_JumpT_FvPatchScalarField::CLL_JumpT_FvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UName_("U"),
    rhoName_("rho"),
    psiName_("thermo:psi"),
    muName_("thermo:mu"),
    accommodationCoeff_t_(1.0),
    accommodationCoeff_n_(1.0),
    Twall_(p.size(), 0.0),
    gamma_(1.4)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::CLL_JumpT_FvPatchScalarField::CLL_JumpT_FvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    psiName_(dict.lookupOrDefault<word>("psi", "thermo:psi")),
    muName_(dict.lookupOrDefault<word>("mu", "thermo:mu")),
    accommodationCoeff_t_(dict.lookup<scalar>("accommodationCoeff_t")),
    accommodationCoeff_n_(dict.lookup<scalar>("accommodationCoeff_n")),
    Twall_("Twall", dict, p.size()),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 1.4))
{
    if
    (
        mag(accommodationCoeff_t_) < small
     || mag(accommodationCoeff_t_) > 2.0
    )
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "unphysical accommodationCoeff_t_ specified"
            << "(0 < accommodationCoeff_t_ <= 1)" << endl
            << exit(FatalIOError);
    }

    if
    (
        mag(accommodationCoeff_n_) < small
     || mag(accommodationCoeff_n_) > 2.0
    )
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "unphysical accommodationCoeff_n_ specified"
            << "(0 < accommodationCoeff_n_ <= 1)" << endl
            << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::CLL_JumpT_FvPatchScalarField::CLL_JumpT_FvPatchScalarField
(
    const CLL_JumpT_FvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    muName_(ptf.muName_),
    accommodationCoeff_t_(ptf.accommodationCoeff_t_),
    accommodationCoeff_n_(ptf.accommodationCoeff_n_),
    Twall_(mapper(ptf.Twall_)),
    gamma_(ptf.gamma_)
{}


Foam::CLL_JumpT_FvPatchScalarField::CLL_JumpT_FvPatchScalarField
(
    const CLL_JumpT_FvPatchScalarField& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    accommodationCoeff_t_(ptpsf.accommodationCoeff_t_),
    accommodationCoeff_n_(ptpsf.accommodationCoeff_n_),
    Twall_(ptpsf.Twall_),
    gamma_(ptpsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CLL_JumpT_FvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    m(Twall_, Twall_);
}


void Foam::CLL_JumpT_FvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const CLL_JumpT_FvPatchScalarField& ptpsf =
        refCast<const CLL_JumpT_FvPatchScalarField>(ptf);

    Twall_.rmap(ptpsf.Twall_, addr);
}


void Foam::CLL_JumpT_FvPatchScalarField::reset
(
    const fvPatchField<scalar>& ptf
)
{
    mixedFvPatchField<scalar>::reset(ptf);

    const CLL_JumpT_FvPatchScalarField& ptpsf =
        refCast<const CLL_JumpT_FvPatchScalarField>(ptf);

    Twall_.reset(ptpsf.Twall_);
}


void Foam::CLL_JumpT_FvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& pmu =
        patch().lookupPatchField<volScalarField, scalar>(muName_);
    const fvPatchScalarField& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);
    const fvPatchField<scalar>& ppsi =
        patch().lookupPatchField<volScalarField, scalar>(psiName_);
    const fvPatchVectorField& pU =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    // Prandtl number reading consistent with rhoCentralFoam
    const dictionary& thermophysicalProperties =
        db().lookupObject<IOdictionary>(physicalProperties::typeName);

    dimensionedScalar Pr
    (
        "Pr",
        dimless,
        thermophysicalProperties.subDict("mixture").subDict("transport")
        .lookup("Pr")
    );

    Field<scalar> C2
    (
        pmu/prho
        *sqrt(ppsi*constant::mathematical::piByTwo)
        *2.0*gamma_/Pr.value()/(gamma_ + 1.0)
        *2.0*(10.0-3.0*accommodationCoeff_n_-2.0*accommodationCoeff_t_*(2.0-accommodationCoeff_t_))
        /(5.0*(accommodationCoeff_n_+accommodationCoeff_t_*(2.0-accommodationCoeff_t_)))
    );

    Field<scalar> aCoeff(prho.snGrad() - prho/C2);
    Field<scalar> KEbyRho(0.5*magSqr(pU));

    valueFraction() = (1.0/(1.0 + patch().deltaCoeffs()*C2));
    refValue() = Twall_;
    refGrad() = 0.0;

    // ------------------------------------------------------------------------
    // // 输出温度跳跃系数
    // Field<scalar> Ct
    // (
    //     sqrt(1/(2.0*ppsi)) * C2 / pmu * prho
    // );

    // Info << "Twall_:" << Twall_ << endl;       
    // Info << "Ct:" << *Ct.begin() << endl;
    // -------------------------------------------------------------------------

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::CLL_JumpT_FvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os, "psi", "thermo:psi", psiName_);
    writeEntryIfDifferent<word>(os, "mu", "thermo:mu", muName_);

    writeEntry(os, "accommodationCoeff_t", accommodationCoeff_t_);
    writeEntry(os, "accommodationCoeff_n", accommodationCoeff_n_);
    writeEntry(os, "Twall", Twall_);
    writeEntry(os, "gamma", gamma_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        CLL_JumpT_FvPatchScalarField
    );
}


// ************************************************************************* //
