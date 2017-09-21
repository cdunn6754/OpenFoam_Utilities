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

Utility
    parcelInfo

Description
   Utility to write out the coal particle mass fraction and char burnout where 
   char burnout is determined by (Y_ash - Y_ash_0)/(1.0 - Y_ash_0).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "basicThermoCloud.H"
#include "DaemCoalCloud.H"
#include "psiCombustionModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    // #include "createControl.H"
    // #include "createTimeControls.H"
    // #include "createRDeltaT.H"
    // #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info <<  mesh.boundary()["sec_inlet"].Cn() << endl;
  
    
  // // Create scalar field to hold Parcel mass fraction

  //   tmp<volScalarField> tparcelMassFraction
  //     (
  //      new volScalarField
  //      (
  // 	IOobject
  // 	(
  // 	 "parcelMassFraction",
  // 	 runTime.timeName(),
  // 	 mesh,
  // 	 IOobject::NO_READ,
  // 	 IOobject::AUTO_WRITE,
  // 	 false
  // 	 ),
  // 	mesh,
  // 	dimensionedScalar("zero", dimMass, 0.0)
  // 	)
  //      );

  //   scalarField& parcelMassFraction = tparcelMassFraction.ref().primitiveFieldRef();

  // // Create scalar field to hold average char burnout
  // volScalarField charBurnout
  //   (
  //    IOobject
  //      (
  //         "charBurnout",
  // 	  runTime.timeName(),
  //         mesh,
  //         IOobject::NO_READ,
  //         IOobject::AUTO_WRITE
  // 	),
  //    mesh,
  //    scalar(0)
  //    );

  // // Need to keep track of number of parcels accounted for 
  // // in each cell to enable an average without reiterating
  // scalarField parcelCount(charBurnout.size(), 0.0);

  // // get the id labels from the thermo composition stuff
  // //const label idGas = coalParcels.composition().idGas();
  // const label idSolid = coalParcels.composition().idSolid();

  // // Mass fractions of Gas and solid that every Particle starts with
  // // const scalar YGas0 = coalParcels.composition().YMixture0()[idGas];
  // // const scalar YSolid0 = coalParcels.composition().YMixture0()[idSolid];

  // // Mass fractions again but itemized by specie
  // // const scalarField& GasFractions0 = coalParcels.composition().Y0(idGas);
  // const scalarField& SolidFractions0 = coalParcels.composition().Y0(idSolid);

  // // Get the ash id from within the SolidFraction
  // const label AshId = coalParcels.composition().localId(idSolid, "ash");

  // // The initial mass fraction of ash within the solid mass fraction
  // const scalar AshFraction0 = SolidFractions0[AshId];


  // // iterate through all parcels in the cloud

  // Info << "Looping through parcels" << endl;
  // forAllIter(typename Cloud<DaemCoalParcel>, coalParcels, pIter)
  //   {
  //     // current parcel
  //     DaemCoalParcel& p = pIter();
  //     // cell of p
  //     const label cellP = p.cell();
  //     //increment the averaging counter
  //     parcelCount[cellP] ++;

  //     // Char Burnout
  //     // ash fraction of current parcel
  //     const scalar AshFraction = p.YSolid()[AshId];
  //     // char burnout fraction
  //     const scalar pCharBurnout = max((AshFraction - AshFraction0) 
  // 				      / (1.0 - AshFraction0), 0.0) ;
  //     // accumulate the average char burnout
  //     // mean(x)_n = mean(x1,...,x_{n-1}) * (n-1)/(n) + (x_n)/n
  //     charBurnout[cellP] = charBurnout[cellP] * 
  // 	(parcelCount[cellP] - 1)/(parcelCount[cellP])
  // 	+ pCharBurnout/parcelCount[cellP];

  //     // Parcel Mass fraction
  //     parcelMassFraction[cellP] += p.nParticle()*p.mass();
  //   }

  // // Adjust to find mass fraction in each cell
  // parcelMassFraction /= (mesh.V() * rho) + tparcelMassFraction.ref();
  

  // Info << "Writing Parcel Mass Fraction to last time step" << endl;
  // tparcelMassFraction.ref().write();

  // Info << "Writing Average Char Burnout to last time step" << endl;
  // charBurnout.write();

  return 0;
}

// ************************************************************************* //
