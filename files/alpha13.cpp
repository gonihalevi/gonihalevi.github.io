//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file alpha13.cpp
//  \brief implementation of functions in class ChemNetwork, using a 13-isotope
//  alpha-chain network
//======================================================================================

// this class header
#include "alpha13.hpp"

//athena++ header
#include "network.hpp"
#include "../../scalars/scalars.hpp"
#include "../../parameter_input.hpp"       //ParameterInput
#include "../../mesh/mesh.hpp"
#include "../../hydro/hydro.hpp"
#include "../utils/chemistry_utils.hpp"
#include "../../defs.hpp"
#include "../../eos/eos.hpp"
#include "../utils/thermo.hpp"

//c++ header
#include <sstream>    // stringstream
#include <iostream>   // endl
#include <limits>    //inf
#include <fstream>   //file()
#include <stdio.h>    // c style file
#include <algorithm>    // std::find()
#include <iterator>     // std::distance()
#include <cmath>       //M_PI

//species names
const std::string ChemNetwork::species_names[NSCALARS] = 
{"4He", "12C", "16O", "20Ne", "24Mg", "28Si", "32S", 
"36Ar", "40Ca", "44Ti", "48Cr", "52Fe", "56Ni"};

const int ChemNetwork::iHe_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "4He");
const int ChemNetwork::iC_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "12C");
const int ChemNetwork::iO_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "16O");
const int ChemNetwork::iNe_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "20Ne");
const int ChemNetwork::iMg_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "24Mg");
const int ChemNetwork::iSi_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "28Si");
const int ChemNetwork::iS_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "32S");
const int ChemNetwork::iAr_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "36Ar");
const int ChemNetwork::iCa_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "40Ca");
const int ChemNetwork::iTi_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "44Ti");
const int ChemNetwork::iCr_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "48Cr");
const int ChemNetwork::iFe_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "52Fe");
const int ChemNetwork::iNi_ =
  ChemistryUtility::FindStrIndex(species_names, NSCALARS, "56Ni");

static const int NISO = 13;
static const int NEQN = 14;
static const int NREAC = 18;
static const int NALP = 38;

/* Nuclear data */
/* Pre-exponential factor in partition function */
static Real g0[NISO];
/* Term in exponential part of partition function divided by T */
static Real apf[NISO];
/* Constant term in exponential part of partition function */
static Real bpf[NISO];
/* Term in exponential part of partition function multiplied by T */
static Real cpf[NISO];
/* Term in screening coefficient */
static Real gscr[NISO];
/* Binding energy */
static Real q[NISO];

/* Reaction rate data */
/* Temperature polynomial coefficients in reaction rate factor */
static Real calp[NALP][7];
/* Constant term in exponential part of reaction rate factor */
static Real ca[NREAC];
/* Energy of reaction; term in exponential part of reaction rate factor
divided by T */
static Real cb[NREAC];


ChemNetwork::ChemNetwork(MeshBlock *pmb, ParameterInput *pin) {
	//number of species and a list of name of species
  pmy_spec_ = pmb->pscalars;
	pmy_mb_ = pmb;
  // gm1_ = pmy_mb_->peos->GetGamma() - 1;

	//set the parameters from input file
  /* Number of the alpha-chain isotope, which is considered "fuel". This isotope
     is part of the cold fuel and it is predominantly consumed. Thus the rate of
     its consumption can be used to calculate quantities, such as the global
     turbulent flame speed. */
  NISOfuel = pin->GetOrAddInteger("chemistry", "NISOfuel", 1);
  /* Increment of energy epsder*e is used in rhands to calculate df(i)/de
     numerically. */
  alphanet_epsder = pin->GetOrAddReal("chemistry","alphanet_epsder",1.e-5);
  //units
	unit_density = pin->GetOrAddReal("chemistry", "unit_density",1.);
	unit_length_in_cm_ = pin->GetOrAddReal("chemistry", "unit_length_in_cm", 1.);
	unit_vel_in_cms_ = pin->GetOrAddReal("chemistry", "unit_vel_in_cms",1.);
  unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  unit_E_in_cgs_ = 1.67e-24 * (gm1_ + 1) * unit_density
                           * unit_vel_in_cms_ * unit_vel_in_cms_;
}

ChemNetwork::~ChemNetwork() {}

void ChemNetwork::InitializeNextStep(const int k, const int j, const int i) {
  const Real five_thirds = 5.0 / 3.0;
  const Real conv_factor = 9.867425e9;
  Real rho, rho_floor;
  //density
  rho = pmy_mb_->phydro->w(IDN, k, j, i);
  //apply density floor
  rho_floor = pmy_mb_->peos->GetDensityFloor();
  rho = (rho > rho_floor) ?  rho : rho_floor;
  //density in proper units
  rho_ =  rho * unit_density;

  int MAXLEN = 256;
  char buf[MAXLEN];
  int ai, ax;
  int l, m;
  Real z[NISO];
  std::string namei;
  std::string namex;

  /* Entering nuclear data table */
  std::ifstream nuc_data("/Users/ghalevi/athena/src/chemistry/network/alpnet.dat");
  if (!nuc_data) {
    std::cout << "### FATAL ERROR in ChemNetwork::InitializeNextStep" << "\n"
        << "Unable to open alpnet.dat" << "\n";
  }
  std::string line;
  int nline = 0;
  while (nline < NISO) {
    getline(nuc_data, line);
    // std::cout << nline << ": " << line << "\n";
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::copy(std::istream_iterator<std::string>(iss),
     std::istream_iterator<std::string>(),
     std::back_inserter(tokens));
    l = 1;
    if (nline == 0) {l++;}
    z[nline] = std::stod(tokens[l]);
    q[nline] = std::stod(tokens[l+1]);
    g0[nline] = std::stod(tokens[l+2]);
    apf[nline] = std::stod(tokens[l+3]);
    bpf[nline] = std::stod(tokens[l+4]);
    cpf[nline] = std::stod(tokens[l+5]);
    nline++;
  }
  // for (i = 0; i < NISO; ++i) {
  //   std::cout << z[i] << ", " << q[i] << ", " << g0[i] << ", " << apf[i] << ", " <<
  //     bpf[i] << ", " << cpf[i] << "\n \n";
  // }

  /* Entering reaction rate data table */
  nline = 0;
  m = 0;
  while (getline(nuc_data, line)) {
    // std::cout << nline << ": " << line << "\n";
    if (nline == 0) {
      /* 3He ==> C */
      for (l = 0; l < 6; ++l) {
        calp[0][l] = std::stod(line.substr(l*13, 13));
        // std::cout << line.substr(i*13, 13) << "\n";
      }
    }
    else if (nline == 1) {
      calp[0][6] = std::stod(line.substr(0,13));
      calp[0][0] -= log(6.0); // to make sure it only happens once
      // for (i = 0; i < 7; ++i) {
      //   std::cout << calp[0][i] << ' ';
      // }
      // std::cout << "\n";
    }
    else {
      // /* a(x,y)b reactions */
      if (nline % 2 == 0) {
      //   /* Read the first line of this entry */
        namex = line.substr(1,2).c_str();
        ax = atoi(line.substr(3,2).c_str());
        namei = line.substr(6,2).c_str();
        ai = atoi(line.substr(8,2).c_str());
        for (l = 0; l < 4; ++l) {
          calp[m][l] = std::stod(line.substr(20+l*14,14));
        }
      }
      else if (nline %2 != 0) {
        for (l = 0; l < 3; ++l) {
          calp[m][l+4] = std::stod(line.substr(20+l*14,14));
        }
        if ((ax == ai) && (namex == namei)) {
          calp[m][0] -= log(2.0);
        }

        // for (i = 0 ; i < 7; i++) {
        //     std::cout << calp[j][i] << "\n";
        // }
        m++;
      }
    }
    nline++;
  }

  nuc_data.close();

  // /* Calculation of the reverse reaction rate coefficients */

  /* C ==> 3 He */
  ca[0] = log(1.199252e21);
  cb[0] = 3.0 * q[0] - q[1];

  /* Ne + He ==> C + C */
  ca[1] = 1.5 * log(12.0 * 12.0 / (4.0 * 20.0));
  cb[1] = q[1] + q[1] - q[0] - q[3];

  /* Mg ==> C + C */
  ca[2] = log(conv_factor * pow(12.0 * 12.0 / 24, 1.5));
  cb[2] = q[1] + q[1] - q[4];

  /* Mg + He ==> C + O */
  ca[3] = 1.5 * log(12.0 * 16.0 / (4.0 * 24.0));
  cb[3] = q[1] + q[2] - q[0] - q[4];

  /* Si ==> C + O */
  ca[4] = log(conv_factor * pow(12.0 * 16.0 / 28.0, 1.5));
  cb[4] = q[1] + q[2] - q[5];

  /* Si + He ==> O + O */
  ca[5] = 1.5 * log(16.0 * 16.0 / (4.0 * 28.0));
  cb[5] = q[2] + q[2] - q[0] - q[5];

  /* S ==> O + O */
  ca[6] = log(conv_factor * pow(16.0 * 16.0 / 32.0, 1.5));
  cb[6] = q[2] + q[2] - q[6];

  /* O ==> He + C */
  ca[7] = log(conv_factor * pow(12.0 * 4.0 / 16.0, 1.5));
  cb[7] = q[0] + q[1] - q[2];

  /* Ne ==> He + O */
  ca[8] = log(conv_factor * pow(16.0 * 4.0 / 20.0, 1.5));
  cb[8] = q[0] + q[2] - q[3];

  /* Mg ==> He + Ne */
  ca[9] = log(conv_factor * pow(20.0 * 4.0 / 24.0, 1.5));
  cb[9] = q[0] + q[3] - q[4];

  /* Si ==> He + Mg */
  ca[10] = log(conv_factor * pow(24.0 * 4.0 / 28.0, 1.5));
  cb[10] = q[0] + q[4] - q[5];

  /* S  ==> He + Si */
  ca[11] = log(conv_factor * pow(28.0 * 4.0 / 32.0, 1.5));
  cb[11] = q[0] + q[5] - q[6];

  /* Ar ==> He + S */
  ca[12] = log(conv_factor * pow(32.0 * 4.0 / 36.0, 1.5));
  cb[12] = q[0] + q[6] - q[7];

  /* Ca ==> He + Ar */
  ca[13] = log(conv_factor * pow(36.0 * 4.0 / 40.0, 1.5));
  cb[13] = q[0] + q[7] - q[8];

  /* Ti ==> He + Ca */
  ca[14] = log(conv_factor * pow(40.0 * 4.0 / 44.0, 1.5));
  cb[14] = q[0] + q[8] - q[9];

  /* Cr ==> He + Ti */
  ca[15] = log(conv_factor * pow(44.0 * 4.0 / 48.0, 1.5));
  cb[15] = q[0] + q[9] - q[10];

  /* Fe ==> He + Cr */
  ca[16] = log(conv_factor * pow(48.0 * 4.0 / 52.0, 1.5));
  cb[16] = q[0] + q[10] - q[11];

  /* Ni ==> He + Fe */
  ca[17] = log(conv_factor * pow(52.0 * 4.0 / 56.0, 1.5));
  cb[17] = q[0] + q[11] - q[12];

  /* Set screening coefficients */
  for (l = 0; l < NISO; ++l) {
    gscr[l] = 0.2275e-3 * pow(z[l], five_thirds);
  }
  return;
}

void ChemNetwork::CalculateRates(Real rho, Real tp, Real frv[NREAC], Real rev[NREAC]){
  /* Parameters for screening corrections */
  const Real a1 = -0.897744;
  const Real a2 =  4.0 * 0.95043;
  const Real a3 = -4.0 * 0.18956;
  const Real a4 = -0.81487;
  const Real a5 = -2.58020;
  const Real a6 = -0.57735;
  const Real a8 =  2.0160;
  const Real a7 =  0.29341 / a8;

  const Real one_third   = 1.0 /  3.0;
  const Real two_thirds  = 2.0 /  3.0;
  const Real one_twelfth = 1.0 / 12.0;

  int k;
  Real t9, t9i, t9l, t923, t9r, g1, gam, gam4;
  Real pf[NISO];   /* Partition functions */
  Real pf0_inv;    /* Inverse of partition function (formerly) at index 0 */
  Real falp[NALP];
  Real fscr[NISO]; /* Screening factors */

  /* Calculation of forward rates */
  t9   = fmax(0.01, 1.e-9 * tp);
  t9i  = 1.0 / t9;
  t9l  = log(t9);
  t923 = pow(t9, two_thirds);
  t9r  = 11.605 * t9i;

  for (k = 0; k < NALP; ++k) {
    falp[k] =   calp[k][0] + 
        t9i  * (calp[k][1] + 
        t923 * (calp[k][2] + 
        t923 * (calp[k][3] + 
        t923 * (calp[k][4] + 
        t923 * (calp[k][5] ))))) + 
        t9l  *  calp[k][6];
  }
  for (k = 0; k < NALP; ++k) {
    falp[k] = rho * exp(falp[k]);
  }

  frv[ 0] =  falp[0] * rho * one_twelfth;
  frv[ 1] =  falp[1] * 0.5;
  frv[ 2] = (falp[2] + falp[3]) * 0.5;
  frv[ 3] =  falp[4];
  frv[ 4] =  falp[5] + falp[6];
  frv[ 5] =  falp[7] * 0.5;
  frv[ 6] = (falp[8] + falp[9]) * 0.5;
  frv[ 7] =  falp[10] + falp[11];
  frv[ 8] =  falp[12] + falp[13];
  frv[ 9] =  falp[14] + falp[15] + falp[16] + falp[17] + falp[18];
  frv[10] =  falp[19] + falp[20] + falp[21] + falp[22] + falp[23];
  frv[11] =  falp[24] + falp[25];
  frv[12] =  falp[26] + falp[27];
  frv[13] =  falp[28] + falp[29];
  frv[14] =  falp[30] + falp[31];
  frv[15] =  falp[32] + falp[33];
  frv[16] =  falp[34] + falp[35];
  frv[17] =  falp[36] + falp[37];

  t9  = 1.e-9 * tp;
  t9i = 1.0 / t9;
  t9r = 11.605 * t9i;

  /* Screening corrections to the forward rates */
  g1 = t9i * pow(0.5 * rho, one_third);
  for (k = 0; k < NISO; ++k) {
    gam = fmin(150.0, g1 * gscr[k]);
    if (gam < 1.0) {
      fscr[k] = a6 * gam * sqrt(gam) + a7 * pow(gam, a8);
    } else {
      gam4 = sqrt(sqrt(gam));
      fscr[k] = a1 * gam + a2 * gam4 + a3 / gam4 + a4 * log(gam) + a5;
    }
  }

  frv[ 0] *= exp(3.0 * fscr[0]            - fscr[1]);
  frv[ 1] *= exp(2.0 * fscr[1]            - fscr[4]);
  frv[ 2] *= exp(2.0 * fscr[1]            - fscr[4]);
  frv[ 3] *= exp(      fscr[1] + fscr[ 2] - fscr[5]);
  frv[ 4] *= exp(      fscr[1] + fscr[ 2] - fscr[5]);
  frv[ 5] *= exp(2.0 * fscr[2]            - fscr[6]);
  frv[ 6] *= exp(2.0 * fscr[2]            - fscr[6]);
  frv[ 7] *= exp(      fscr[0] + fscr[ 1] - fscr[2]);
  frv[ 8] *= exp(      fscr[0] + fscr[ 2] - fscr[3]);
  frv[ 9] *= exp(      fscr[0] + fscr[ 3] - fscr[4]);
  frv[10] *= exp(      fscr[0] + fscr[ 4] - fscr[5]);
  frv[11] *= exp(      fscr[0] + fscr[ 5] - fscr[6]);
  frv[12] *= exp(      fscr[0] + fscr[ 6] - fscr[7]);
  frv[13] *= exp(      fscr[0] + fscr[ 7] - fscr[8]);
  frv[14] *= exp(      fscr[0] + fscr[ 8] - fscr[9]);
  frv[15] *= exp(      fscr[0] + fscr[ 9] - fscr[10]);
  frv[16] *= exp(      fscr[0] + fscr[10] - fscr[11]);
  frv[17] *= exp(      fscr[0] + fscr[11] - fscr[12]);

  /* Calculation of partition functions */
  for (k = 0; k < NISO; ++k) {
    pf[k] = g0[k] * (1.0 + exp(apf[k] * t9i + bpf[k] + t9 * cpf[k]));
  }
  pf0_inv = t9 * sqrt(t9) / rho;

  /* Calculation of the reverse rate coefficients */

  /* 3 He ==> C */
  rev[ 0] = exp(ca[ 0] + cb[ 0] * t9r + fscr[1] - 3.0 * fscr[0])
          * pf[0] * pf[0] * pf[0] * t9 * t9 * t9 / (pf[1] * rho * rho);

  /* Ne + He ==> C + C */
  rev[ 1] = exp(ca[ 1] + cb[ 1] * t9r + 2.0 * fscr[1] - fscr[0] - fscr[3])
          * pf[1] * pf[1] / (pf[3] * pf[0]);

  /* Mg ==> C + C */
  rev[ 2] = exp(ca[ 2] + cb[ 2] * t9r + 2.0 * fscr[1] - fscr[4])
          * pf[1] * pf[1] * pf0_inv / pf[4];

  /* Mg + He ==> C + O */
  rev[ 3] = exp(ca[ 3] + cb[ 3] * t9r + fscr[1] + fscr[2] - fscr[0] - fscr[4])
          * pf[1] * pf[2] / (pf[4] * pf[0]);

  /* Si ==> C + O */
  rev[ 4] = exp(ca[ 4] + cb[ 4] * t9r + fscr[1] + fscr[2] - fscr[5])
          * pf[1] * pf[2] * pf0_inv / pf[5];

  /* Si + He ==> O + O */
  rev[ 5] = exp(ca[ 5] + cb[ 5] * t9r + 2.0 * fscr[2] - fscr[0] - fscr[5])
          * pf[2] * pf[2] / (pf[5] * pf[0]);

  /* S ==> O + O */
  rev[ 6] = exp(ca[ 6] + cb[ 6] * t9r + 2.0 * fscr[2] - fscr[6])
          * pf[2] * pf[2] * pf0_inv / pf[6];

  /* O ==> C + He */
  rev[ 7] = exp(ca[ 7] + cb[ 7] * t9r + fscr[1] + fscr[0] - fscr[2])
          * pf[0] * pf0_inv * pf[1] / pf[2];

  /* Ne ==> O + He */
  rev[ 8] = exp(ca[ 8] + cb[ 8] * t9r + fscr[2] + fscr[0] - fscr[3])
          * pf[0] * pf0_inv * pf[2] / pf[3];

  /* Mg ==> Ne + He */
  rev[ 9] = exp(ca[ 9] + cb[ 9] * t9r + fscr[3] + fscr[0] - fscr[4])
          * pf[0] * pf0_inv * pf[3] / pf[4];

  /* Si ==> Mg + He */
  rev[10] = exp(ca[10] + cb[10] * t9r + fscr[4] + fscr[0] - fscr[5])
          * pf[0] * pf0_inv * pf[4] / pf[5];

  /* S  ==> Si + He */
  rev[11] = exp(ca[11] + cb[11] * t9r + fscr[5] + fscr[0] - fscr[6])
          * pf[0] * pf0_inv * pf[5] / pf[6];

  /* Ar ==> S  + He */
  rev[12] = exp(ca[12] + cb[12] * t9r + fscr[6] + fscr[0] - fscr[7])
          * pf[0] * pf0_inv * pf[6] / pf[7];

  /* Ca ==> Ar + He */
  rev[13] = exp(ca[13] + cb[13] * t9r + fscr[7] + fscr[0] - fscr[8])
          * pf[0] * pf0_inv * pf[7] / pf[8];

  /* Ti ==> Ca + He */
  rev[14] = exp(ca[14] + cb[14] * t9r + fscr[8] + fscr[0] - fscr[9])
          * pf[0] * pf0_inv * pf[8] / pf[9];

  /* Cr ==> Ti + He */
  rev[15] = exp(ca[15] + cb[15] * t9r + fscr[9] + fscr[0] - fscr[10])
          * pf[0] * pf0_inv * pf[9] / pf[10];

  /* Fe ==> Cr + He */
  rev[16] = exp(ca[16] + cb[16] * t9r + fscr[10] + fscr[0] - fscr[11])
          * pf[0] * pf0_inv * pf[10] / pf[11];

  /* Ni ==> Fe + He */
  rev[17] = exp(ca[17] + cb[17] * t9r + fscr[11] + fscr[0] - fscr[12])
          * pf[0] * pf0_inv * pf[11] / pf[12];

}

void ChemNetwork::RatesOfChange(const Real frv[NREAC], const Real rev[NREAC],
  const Real y[NSCALARS], Real f[NSCALARS])
{
  int i;
  Real r;    /* Reaction rate */

  /* 3He ==> C */
  r = frv[0] * (y[0] * y[0] * y[0] - rev[0] * y[1]);
  f[0] = -r * 3.0;
  f[1] =  r;

  /* C + C ==> Ne + He */
  r = frv[1] * (y[1] * y[1] - rev[1] * y[0] * y[3]);
  f[1] -= r * 2.0;
  f[0] += r;
  f[3]  = r;

  /* C + C ==> Mg */
  r = frv[2] * (y[1] * y[1] - rev[2] * y[4]);
  f[1] -= r * 2.0;
  f[4]  = r;

  /* C + O ==> He + Mg */
  r = frv[3] * (y[1] * y[2] - rev[3] * y[0] * y[4]);
  f[1] -=  r;
  f[2]  = -r;
  f[0] +=  r;
  f[4] +=  r;

  /* C + O ==> Si */
  r = frv[4] * (y[1] * y[2] - rev[4] * y[5]);
  f[1] -= r;
  f[2] -= r;
  f[5]  = r;

  /* O + O ==> Si + He */
  r = frv[5] * (y[2] * y[2] - rev[5] * y[0] * y[5]);
  f[2] -= r * 2.0;
  f[0] += r;
  f[5] += r;

  /* O + O ==> S */
  r = frv[6] * (y[2] * y[2] - rev[6] * y[6]);
  f[2] -= r * 2.0;
  f[6]  = r;

  /* He + C ==> O */
  r = frv[7] * (y[0] * y[1] - rev[7] * y[2]);
  f[0] -= r;
  f[1] -= r;
  f[2] += r;

  /* He + O ==> Ne */
  r = frv[8] * (y[0] * y[2] - rev[8] * y[3]);
  f[0] -= r;
  f[2] -= r;
  f[3] += r;

  /* He + Ne ==> Mg */
  r = frv[9] * (y[0] * y[3] - rev[9] * y[4]);
  f[0] -= r;
  f[3] -= r;
  f[4] += r;

  /* He + Mg ==> Si */
  r = frv[10] * (y[0] * y[4] - rev[10] * y[5]);
  f[0] -= r;
  f[4] -= r;
  f[5] += r;

  /* He + Si ==> S */
  r = frv[11] * (y[0] * y[5] - rev[11] * y[6]);
  f[0] -= r;
  f[5] -= r;
  f[6] += r;

  /* He + S ==> Ar */
  r = frv[12] * (y[0] * y[6] - rev[12] * y[7]);
  f[0] -=  r;
  f[6] -=  r;
  f[7]  =  r;

  /* He + Ar ==> Ca */
  r = frv[13] * (y[0] * y[7] - rev[13] * y[8]);
  f[0] -=  r;
  f[7] -=  r;
  f[8]  =  r;

  /* He + Ca ==> Ti */
  r = frv[14] * (y[0] * y[8] - rev[14] * y[9]);
  f[0] -=  r;
  f[8] -=  r;
  f[9]  =  r;

  /* He + Ti ==> Cr */
  r = frv[15] * (y[0] * y[9] - rev[15] * y[10]);
  f[0]  -=  r;
  f[9]  -=  r;
  f[10]  =  r;

  /* He + Cr ==> Fe */
  r = frv[16] * (y[0] * y[10] - rev[16] * y[11]);
  f[0]  -=  r;
  f[10] -=  r;
  f[11]  =  r;

  /* He + Fe ==> Ni */
  r = frv[17] * (y[0] * y[11] - rev[17] * y[12]);
  f[0]  -=  r;
  f[11] -=  r;
  f[12]  =  r;
}

void ChemNetwork::PartialDerivatives(const Real frv[NREAC], const Real rev[NREAC],
  const Real y[NSCALARS], Real f[NEQN], Real df[NEQN][NEQN])
{
 const Real conv_factor = 9.64867e17;

  /* Reaction rate */
  Real r;

  /* Partial derivatives of reaction rate wrt species */
  Real r_0, r_1, r_2, r_3, r_4, r_5, r_6, r_7, r_8, r_9, r_10, r_11, r_12;

  /* Reaction rate */
  Real edot;

  int j, k;

  for (j = 0; j < NEQN; ++j) {
    for (k = 0; k < NEQN; ++k) {
      df[j][k] = 0.0;
    }
  }

  /* 3He ==> C */
  r = frv[0] * (y[0] * y[0] * y[0] - rev[0] * y[1]);
  f[0] = -r * 3.0;
  f[1] =  r;
  r_0 =  frv[0] * y[0] * y[0];
  r_1 = -frv[0] * rev[0];
  df[0][0] -= r_0 * 9.0;
  df[1][0] -= r_1 * 3.0;
  df[0][1] += r_0 * 3.0;
  df[1][1] += r_1;

  /* C + C ==> Ne + He */
  r = frv[1] * (y[1] * y[1] - rev[1] * y[0] * y[3]);
  f[1] -= r * 2.0;
  f[0] += r;
  f[3]  = r;
  r_1 =  frv[1] * y[1];
  r_0 = -frv[1] * rev[1] * y[3];
  r_3 = -frv[1] * rev[1] * y[0];
  df[0][1] -= r_0 * 2.0;
  df[1][1] -= r_1 * 4.0;
  df[3][1] -= r_3 * 2.0;
  df[0][0] += r_0;
  df[1][0] += r_1 * 2.0;
  df[3][0] += r_3;
  df[0][3] += r_0;
  df[1][3] += r_1 * 2.0;
  df[3][3] += r_3;

  /* C + C ==> Mg */
  r = frv[2] * (y[1] * y[1] - rev[2] * y[4]);
  f[1] -= r * 2.0;
  f[4]  = r;
  r_1 =  frv[2] * y[1];
  r_4 = -frv[2] * rev[2];
  df[1][1] -= r_1 * 4.0;
  df[4][1] -= r_4 * 2.0;
  df[1][4] += r_1 * 2.0;
  df[4][4] += r_4;

  /* C + O ==> He + Mg */
  r = frv[3] * (y[1] * y[2] - rev[3] * y[0] * y[4]);
  f[1] -=  r;
  f[2]  = -r;
  f[0] +=  r;
  f[4] +=  r;
  r_1 =  frv[3] * y[2];
  r_2 =  frv[3] * y[1];
  r_0 = -frv[3] * rev[3] * y[4];
  r_4 = -frv[3] * rev[3] * y[0];
  df[1][1] -= r_1;
  df[2][1] -= r_2;
  df[0][1] -= r_0;
  df[4][1] -= r_4;
  df[1][2] -= r_1;
  df[2][2] -= r_2;
  df[0][2] -= r_0;
  df[4][2] -= r_4;
  df[1][0] += r_1;
  df[2][0] += r_2;
  df[0][0] += r_0;
  df[4][0] += r_4;
  df[1][4] += r_1;
  df[2][4] += r_2;
  df[0][4] += r_0;
  df[4][4] += r_4;

  /* C + O ==> Si */
  r = frv[4] * (y[1] * y[2] - rev[4] * y[5]);
  f[1] -= r;
  f[2] -= r;
  f[5]  = r;
  r_1 =  frv[4] * y[2];
  r_2 =  frv[4] * y[1];
  r_5 = -frv[4] * rev[4];
  df[1][1] -= r_1;
  df[2][1] -= r_2;
  df[5][1] -= r_5;
  df[1][2] -= r_1;
  df[2][2] -= r_2;
  df[5][2] -= r_5;
  df[1][5] += r_1;
  df[2][5] += r_2;
  df[5][5] += r_5;

  /* O + O ==> Si + He */
  r = frv[5] * (y[2] * y[2] - rev[5] * y[0] * y[5]);
  f[2] -= r * 2.0;
  f[0] += r;
  f[5] += r;
  r_2 =  frv[5] * y[2];
  r_0 = -frv[5] * rev[5] * y[5];
  r_5 = -frv[5] * rev[5] * y[0];
  df[2][2] -= r_2 * 4.0;
  df[0][2] -= r_0 * 2.0;
  df[5][2] -= r_5 * 2.0;
  df[2][0] += r_2 * 2.0;
  df[0][0] += r_0;
  df[5][0] += r_5;
  df[2][5] += r_2 * 2.0;
  df[0][5] += r_0;
  df[5][5] += r_5;

  /* O + O ==> S */
  r = frv[6] * (y[2] * y[2] - rev[6] * y[6]);
  f[2] -= r * 2.0;
  f[6]  = r;
  r_2 =  frv[6] * y[2];
  r_6 = -frv[6] * rev[6];
  df[2][2] -= r_2 * 4.0;
  df[6][2] -= r_6 * 2.0;
  df[2][6] += r_2 * 2.0;
  df[6][6] += r_6;

  /* He + C ==> O */
  r = frv[7] * (y[0] * y[1] - rev[7] * y[2]);
  f[0] -= r;
  f[1] -= r;
  f[2] += r;
  r_0 =  frv[7] * y[1];
  r_1 =  frv[7] * y[0];
  r_2 = -frv[7] * rev[7];
  df[0][0] -= r_0;
  df[1][0] -= r_1;
  df[2][0] -= r_2;
  df[0][1] -= r_0;
  df[1][1] -= r_1;
  df[2][1] -= r_2;
  df[0][2] += r_0;
  df[1][2] += r_1;
  df[2][2] += r_2;

  /* He + O ==> Ne */
  r = frv[8] * (y[0] * y[2] - rev[8] * y[3]);
  f[0] -= r;
  f[2] -= r;
  f[3] += r;
  r_0 =  frv[8] * y[2];
  r_2 =  frv[8] * y[0];
  r_3 = -frv[8] * rev[8];
  df[0][0] -= r_0;
  df[2][0] -= r_2;
  df[3][0] -= r_3;
  df[0][2] -= r_0;
  df[2][2] -= r_2;
  df[3][2] -= r_3;
  df[0][3] += r_0;
  df[2][3] += r_2;
  df[3][3] += r_3;

  /* He + Ne ==> Mg */
  r = frv[9] * (y[0] * y[3] - rev[9] * y[4]);
  f[0] -= r;
  f[3] -= r;
  f[4] += r;
  r_0 =  frv[9] * y[3];
  r_3 =  frv[9] * y[0];
  r_4 = -frv[9] * rev[9];
  df[0][0] -= r_0;
  df[3][0] -= r_3;
  df[4][0] -= r_4;
  df[0][3] -= r_0;
  df[3][3] -= r_3;
  df[4][3] -= r_4;
  df[0][4] += r_0;
  df[3][4] += r_3;
  df[4][4] += r_4;

  /* He + Mg ==> Si */
  r = frv[10] * (y[0] * y[4] - rev[10] * y[5]);
  f[0] -= r;
  f[4] -= r;
  f[5] += r;
  r_0 =  frv[10] * y[4];
  r_4 =  frv[10] * y[0];
  r_5 = -frv[10] * rev[10];
  df[0][0] -= r_0;
  df[4][0] -= r_4;
  df[5][0] -= r_5;
  df[0][4] -= r_0;
  df[4][4] -= r_4;
  df[5][4] -= r_5;
  df[0][5] += r_0;
  df[4][5] += r_4;
  df[5][5] += r_5;

  /* He + Si ==> S */
  r = frv[11] * (y[0] * y[5] - rev[11] * y[6]);
  f[0] -= r;
  f[5] -= r;
  f[6] += r;
  r_0 =  frv[11] * y[5];
  r_5 =  frv[11] * y[0];
  r_6 = -frv[11] * rev[11];
  df[0][0] -= r_0;
  df[5][0] -= r_5;
  df[6][0] -= r_6;
  df[0][5] -= r_0;
  df[5][5] -= r_5;
  df[6][5] -= r_6;
  df[0][6] += r_0;
  df[5][6] += r_5;
  df[6][6] += r_6;

  /* He + S ==> Ar */
  r = frv[12] * (y[0] * y[6] - rev[12] * y[7]);
  f[0] -=  r;
  f[6] -=  r;
  f[7]  =  r;
  r_0 =  frv[12] * y[6];
  r_6 =  frv[12] * y[0];
  r_7 = -frv[12] * rev[12];
  df[0][0] -=  r_0;
  df[6][0] -=  r_6;
  df[7][0] -=  r_7;
  df[0][6] -=  r_0;
  df[6][6] -=  r_6;
  df[7][6] -=  r_7;
  df[0][7] +=  r_0;
  df[6][7] +=  r_6;
  df[7][7] +=  r_7;

  /* He + Ar ==> Ca */
  r = frv[13] * (y[0] * y[7] - rev[13] * y[8]);
  f[0] -=  r;
  f[7] -=  r;
  f[8]  =  r;
  r_0 =  frv[13] * y[7];
  r_7 =  frv[13] * y[0];
  r_8 = - frv[13] * rev[13];
  df[0][0] -=  r_0;
  df[7][0] -=  r_7;
  df[8][0] -=  r_8;
  df[0][7] -=  r_0;
  df[7][7] -=  r_7;
  df[8][7] -=  r_8;
  df[0][8] +=  r_0;
  df[7][8] +=  r_7;
  df[8][8] +=  r_8;

  /* He + Ca ==> Ti */
  r = frv[14] * (y[0] * y[8] - rev[14] * y[9]);
  f[0] -=  r;
  f[8] -=  r;
  f[9]  =  r;
  r_0 =  frv[14] * y[8];
  r_8 =  frv[14] * y[0];
  r_9 = -frv[14] * rev[14];
  df[0][0] -=  r_0;
  df[8][0] -=  r_8;
  df[9][0] -=  r_9;
  df[0][8] -=  r_0;
  df[8][8] -=  r_8;
  df[9][8] -=  r_9;
  df[0][9] +=  r_0;
  df[8][9] +=  r_8;
  df[9][9] +=  r_9;

  /* He + Ti ==> Cr */
  r = frv[15] * (y[0] * y[9] - rev[15] * y[10]);
  f[ 0] -=  r;
  f[ 9] -=  r;
  f[10]  =  r;
  r_0  =  frv[15] * y[9];
  r_9  =  frv[15] * y[0];
  r_10 = -frv[15] * rev[15];
  df[ 0][ 0] -=  r_0;
  df[ 9][ 0] -=  r_9;
  df[10][ 0] -=  r_10;
  df[ 0][ 9] -=  r_0;
  df[ 9][ 9] -=  r_9;
  df[10][ 9] -=  r_10;
  df[ 0][10] +=  r_0;
  df[ 9][10] +=  r_9;
  df[10][10] +=  r_10;

  /* He + Cr ==> Fe */
  r = frv[16] * (y[0] * y[10] - rev[16] * y[11]);
  f[ 0] -=  r;
  f[10] -=  r;
  f[11]  =  r;
  r_0  =  frv[16] * y[10];
  r_10 =  frv[16] * y[0];
  r_11 = -frv[16] * rev[16];
  df[ 0][ 0] -=  r_0;
  df[10][ 0] -=  r_10;
  df[11][ 0] -=  r_11;
  df[ 0][10] -=  r_0;
  df[10][10] -=  r_10;
  df[11][10] -=  r_11;
  df[ 0][11] +=  r_0;
  df[10][11] +=  r_10;
  df[11][11] +=  r_11;

  /* He + Fe ==> Ni */
  r = frv[17] * (y[0] * y[11] - rev[17] * y[12]);
  f[ 0] -=  r;
  f[11] -=  r;
  f[12]  =  r;
  r_0  =  frv[17] * y[11];
  r_11 =  frv[17] * y[0];
  r_12 = -frv[17] * rev[17];
  df[0 ][ 0] -=  r_0;
  df[11][ 0] -=  r_11;
  df[12][ 0] -=  r_12;
  df[ 0][11] -=  r_0;
  df[11][11] -=  r_11;
  df[12][11] -=  r_12;
  df[ 0][12] +=  r_0;
  df[11][12] +=  r_11;
  df[12][12] +=  r_12;

  edot = 0.0;
  for (k = 0; k < NISO; ++k) {
    edot += q[k] * f[k];
  }
  f[NEQN-1] = conv_factor * edot;

  for (j = 0; j < NISO; ++j) {
    edot = 0.0;
    for (k = 0; k < NISO; ++k) {
      edot += q[k] * df[j][k];
    }
    df[j][NEQN-1] = conv_factor * edot;
  }

}


void ChemNetwork::RHS(const Real t, const Real y[NSCALARS], const Real ED,
                      Real ydot[NSCALARS])
{
  Real frv[NREAC]; /* Forward reaction rates */
  Real rev[NREAC]; /* Reverse reaction rates */
  Real f[NEQN]; //rates of change; last element is de/dt
  // assumes adiabatic
  const Real gm1 = pmy_mb_->peos->GetGamma() - 1;
  Real mn_ = 1.674920e-24;
  Real k_ = 1.380658e-16;
  Real temp = ED*gm1*mn_/(rho_*k_);
  if (t == 0) {
    printf("T_9 = %f \n", temp*1e-9);
  }
  Real y_corr[NSCALARS];
  Real y_floor = 0.0;
  for (int i=0; i<NSCALARS; i++) {
    if (y[i] < y_floor) {
      y_corr[i] = y_floor;
    } else {
      y_corr[i] = y[i];
    }
  }
  #ifdef DEBUG
    printf("BEFORE \n");
    for (int i = 0; i <NSCALARS; i++) {
      printf("y[%i] = %4.2f \n", i, y_corr[i]);
    }
  #endif

  CalculateRates(rho_, temp, frv, rev);
  RatesOfChange(frv, rev, y_corr, f);

	for (int i=0; i<NSCALARS; i++) {
    //return in code units
		ydot[i] = unit_time_in_s_*f[i];

    // ydot[i] = 0.0;
	}
  #ifdef DEBUG
    printf("AFTER \n");
    for (int i = 0; i <NREAC; i++) {
      printf("frv[%i] = %4.2f \n", i, frv[i]);
      printf("rev[%i] = %4.2f \n", i, rev[i]);
      printf("ydot[%i ] = %4.2f \n", i, ydot[i]);
    }
  #endif
  return;
}

Real ChemNetwork::Edot(const Real t, const Real y[NSCALARS], const Real ED){
  //isothermal
  if (!NON_BAROTROPIC_EOS) {
    return 0;
  }
  Real ydot[NSCALARS];
  RHS(t, y, ED, ydot);
  const Real conv_factor = 9.64867e17;
  Real edot = 0.0;
  for (int i = 0; i < NSCALARS; ++i) {
    edot += q[i] * ydot[i];
  }
  Real dEDdt;
  dEDdt = conv_factor * edot;
  Real e0_inv = 1.0 / ED;

  /* Scale energy derivative */
  dEDdt *= e0_inv;

  // but this is erg/g/s and we want erg/cm^3/s, right? so multiply by the density...?
  dEDdt *= rho_;
  printf("dEDdt = %4.2f erg\n",i,dEDdt);
  return dEDdt;
}

//calculate Jacobian with numerical differentiation 

// void ChemNetwork::Jacobian_numerical(const Real t, const Real y[NSCALARS],
//                                    const Real ydot[NSCALARS], const Real ED
//                                    AthenaArray<Real> &jac) {
//   int k;
//   Real dens[2];    /* Density [g/cm^3] */
//   Real rhoe[2];    /* Energy density [erg/cm^3] */
//   Real temp[2];    /* Temperature [K] */
//   Real e0_inv;     /* Inverse of energy scale (e0) [g/erg] */
//   Real fn[NEQN];   /* RHS corresponding to perturbed energy */
//   Real frv[NREAC];  /*Forward reaction rates */
//   Real rev[NREAC]; /* Reverse reaction rates */
//   Real e_diff_inv; /* Inverse of the scaled energy perturbation */

//   /* Find temperature of original and perturbed energy density in temp[0]
//      and temp[1], respectively */
//   dens[0] = rho_;
//   dens[1] = rho_;
//   /* Original energy density */
//   rhoe[0] = unit_E_in_cgs_ * ED;
//   // Perturbed energy density.
//   rhoe[1] = rhoe[0] * (1.0 + alphanet_epsder);
//   // assumes adiabatic
//   temp[0] = rhoe[0]*gm1_*conv_factor/(dens[0]);

//   /* Calculate f and jac, except for the partial derivatives wrt energy:
//      jac[NEQN-1][:] */
//   CalculateRates(dens[0], temp[0], frv, rev);
//   PartialDerivatives(frv, rev, y, f, df);

//   e0_inv = 1.0 / unit_E_in_cgs_;

//   /* Scale energy derivatives */
//   f[NEQN-1] *= e0_inv;
//   for (k = 0; k < NEQN - 1; ++k) {
//     jac[k][NEQN-1] *= e0_inv;
//   }

//   /* Calculate last row of jac[NEQN-1][:] numerically */
//   CalculateRates(dens[1], temp[1], frv, rev);
//   RatesOfChange(frv, rev, y, fn);

//   /* Scale energy derivative */
//   fn[NEQN-1] *= e0_inv;

//   /* Calculate numerical derivatives with respect to energy */
//   e_diff_inv = 1.0 / (alphanet_epsder * y[NEQN-1]);
//   for (k = 0; k < NEQN; ++k) {
//     jac[NEQN-1][k] = (fn[k] - f[k]) * e_diff_inv;
//   }
// }
