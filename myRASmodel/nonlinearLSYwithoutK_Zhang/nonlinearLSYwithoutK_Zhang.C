/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "nonlinearLSYwithoutK_Zhang.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> nonlinearLSYwithoutK_Zhang<BasicMomentumTransportModel>::fMu() const
{
    return exp(-3.4/sqr(scalar(1) + sqr(k_)/(this->nu()*epsilon_)/50.0));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> nonlinearLSYwithoutK_Zhang<BasicMomentumTransportModel>::f2() const
{
    return
        scalar(1)
      - 0.3*exp(-min(sqr(sqr(k_)/(this->nu()*epsilon_)), scalar(50.0)));
}


template<class BasicMomentumTransportModel>
void nonlinearLSYwithoutK_Zhang<BasicMomentumTransportModel>::correctNut()
{
    correctNonlinearStress(fvc::grad(this->U_));
}

template<class BasicMomentumTransportModel>
void nonlinearLSYwithoutK_Zhang<BasicMomentumTransportModel>::correctNonlinearStress(const volTensorField& gradU)
{
    dimensionedScalar rhoMin
    (
        "rhoMin",
        dimensionSet(1,-3,0,0,0,0,0),
        SMALL
    );
    
    volSymmTensorField S(2.0*dev(symm(gradU)));
    volTensorField W(2.0*skew(gradU));

    volVectorField gradK(fvc::grad(k_));
    volSymmTensorField gradKgradK(symm(gradK*gradK));
    volSymmTensorField gradGradK(symm(fvc::grad(gradK)));

    volVectorField gradRho(fvc::grad(this->rho_));
    volSymmTensorField gradKgradRho(symm(gradK*gradRho));
    
    this->nut_ = Cmu_*fMu()*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
    volScalarField nut= this->nut_;

    this->nonlinearStress_ =
        sqr(nut)/k_
       *(

            // Quadratic terms
            (
                C1_*dev(innerSqr(S))
              + C2_*twoSymm(S&W)
              + C3_*dev(symm(W&W))
            )

            // Extra TKE terms
        //   + (
        //         C4_/k_*(dev(gradKgradK))
        //       + C5_*(dev(gradGradK))
        //       + C6_/max(this->rho_, rhoMin)*dev(twoSymm(gradKgradRho))               
        //     )
        );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
nonlinearLSYwithoutK_Zhang<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
nonlinearLSYwithoutK_Zhang<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
nonlinearLSYwithoutK_Zhang<BasicMomentumTransportModel>::nonlinearLSYwithoutK_Zhang
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),


    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    Ceps1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1",
            this->coeffDict_,
            1.44
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            1.92
        )
    ),
    C30_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C30",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),

    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.5
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            0.75
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0.0
        )
    ),
    C4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4",
            this->coeffDict_,
            9.0
        )
    ),
    C5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C5",
            this->coeffDict_,
            2.0
        )
    ),
    C6_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C6",
            this->coeffDict_,
            1.0
        )
    ),
    C7_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C7",
            this->coeffDict_,
            0.0
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    nonlinearStress_
    (
        IOobject
        (
            IOobject::groupName("nonlinearStress", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedSymmTensor
        (
            "nonlinearStress",
            sqr(dimVelocity),
            Zero
        )
    ),

    y_
    (
        wallDist::New(this->mesh_).y()
    )

{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool nonlinearLSYwithoutK_Zhang<BasicMomentumTransportModel>::read()
{
    if (nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        Ceps1_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        C30_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        C4_.readIfPresent(this->coeffDict());
        C5_.readIfPresent(this->coeffDict());
        C6_.readIfPresent(this->coeffDict());
        C7_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void nonlinearLSYwithoutK_Zhang<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    // eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();
    nonlinearEddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    // Calculate Yap correction term

    volScalarField lByle((pow(k_,1.5)/epsilon_)/(2.5*y_)); 
    volScalarField Yap_coff // 17/05/2019
    (
        max((lByle - scalar(1.0))*sqr(lByle),scalar(1.0e-16))
    );
    volScalarField Yap 
    (
    	"Yap",
       0.83*Yap_coff*sqr(epsilon_)/k_
    );   

    // Calculate parameters and coefficients for Launder-Sharma low-Reynolds
    // number model

    volScalarField E(2.0*this->nu()*nut*fvc::magSqrGradGrad(U));
    volScalarField D(2.0*this->nu()*magSqr(fvc::grad(sqrt(k_)))); 

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G
    (
        this->GName(), 
        (nut*dev(twoSymm(tgradU())) - this->nonlinearStress_)&&tgradU()  // Added nonlinearStress
        // (nut*dev(twoSymm(gradU)) - this->nonlinearStress_) && gradU
    );

    // Calculate Zhang correction term
    //--------------------------------------
    volSymmTensorField S(2.0*dev(symm(tgradU())));
    volTensorField W(2.0*skew(tgradU()));
    volScalarField sBar((k_/epsilon_)*sqrt(0.5)*mag(S));
    volScalarField wBar((k_/epsilon_)*sqrt(0.5)*mag(W));    
    volScalarField eta(max(sBar,wBar));

    volScalarField lamda(max(k_/epsilon_*divU,-0.5));
    volScalarField Snew(3.0*lamda*sqr(epsilon_)/k_*(tanh(2.0*(eta-3.0))-1.0));
    //--------------------------------------
    tgradU.clear();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        Ceps1_*alpha*rho*G*epsilon_/k_
      - fvm::SuSp(((2.0/3.0)*Ceps1_ - C30_)*alpha*rho*divU, epsilon_)
      - fvm::Sp(Ceps2_*f2()*alpha*rho*epsilon_/k_, epsilon_)
      + alpha*rho*E
      + alpha*rho*Yap  // Added Yap correction term
      + alpha*rho*Snew
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);

    fvConstraints.constrain(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G - fvm::SuSp(2.0/3.0*alpha*rho*divU, k_)
      - fvm::Sp(alpha*rho*(epsilon_ + D)/k_, k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    Info<<"nonlinearLSYwithoutK_Zhang20240125"<<nl<<endl;  

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
