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
    Foam::RASModels::nonlinearCLSwithoutK_Zhang

Description
    Craft cubic non-linear low-Reynolds k-epsilon turbulence models for
    incompressible and compressible flows.

    This turbulence model is described in:
    \verbatim
        T.J. CRAFT, H. IACOVIDES and J.H. YOON (1999).
        Progress in the Use of Non-Linear Two-Equation Models in the 
        Computation of Convective Heat-Transfer in Impinging and Separated 
        Flows. Flow, Turbulence and Combustion 63: 59�?0, 1999.
    \endverbatim

SourceFiles
    nonlinearCLSwithoutK_Zhang.C

\*---------------------------------------------------------------------------*/

#ifndef nonlinearCLSwithoutK_Zhang_H
#define nonlinearCLSwithoutK_Zhang_H

#include "RASModel.H"
// #include "eddyViscosity.H"
#include "nonlinearEddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class nonlinearCLSwithoutK_Zhang Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class nonlinearCLSwithoutK_Zhang
:
    public nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>
{
protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar C30_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar C4_;
            dimensionedScalar C5_;
            dimensionedScalar C6_;
            dimensionedScalar C7_;



        // Fields

            volScalarField k_;
            volScalarField epsilon_;
            volSymmTensorField nonlinearStress_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;

    // Private Member Functions

        tmp<volScalarField> fMu() const;
        tmp<volScalarField> f2() const;

        virtual void correctNut();
        virtual void correctNonlinearStress(const volTensorField& gradU);
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("nonlinearCLSwithoutK_Zhang");


    // Constructors

        //- Construct from components
        nonlinearCLSwithoutK_Zhang
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        nonlinearCLSwithoutK_Zhang(const nonlinearCLSwithoutK_Zhang&) = delete;


    //- Destructor
    virtual ~nonlinearCLSwithoutK_Zhang()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                (this->nut_/sigmak_ + this->nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                (this->nut_/sigmaEps_ + this->nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            // the epsilon returned here is changed 
            // from epsison tilda to "real" epsilon
            // By Haoyuan Zhang at 09/11/2020
            return epsilon_ + 2.0*this->nu()*magSqr(fvc::grad(sqrt(k_)));
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return volScalarField::New
            (
                "omega",
                epsilon_/(Cmu_*k_),
                epsilon_.boundaryField().types()
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nonlinearCLSwithoutK_Zhang&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonlinearCLSwithoutK_Zhang.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
