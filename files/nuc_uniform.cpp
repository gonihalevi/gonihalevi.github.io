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
//! \file nuc_uniform.cpp
//  \brief problem generator, uniform mesh with nuclear reactions
//======================================================================================

// c headers
#include <stdio.h>    // c style file
#include <string.h>   // strcmp()

// C++ headers
#include <algorithm>  // std::find()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // std::runtime_error()
#include <string>     // c_str()
#include <vector>     // vector container

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief initialize problem 
//======================================================================================

Real ReactionTimeStep(MeshBlock *pmb);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  EnrollUserTimeStepFunction(ReactionTimeStep);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  //dimensions of meshblock
  const int Nx = ie - is + 1;
  const int Ny = je - js + 1;
  const int Nz = ke - ks + 1;
  //read density and radiation field strength
  const Real rho = pin->GetReal("problem", "rho");
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real vx = pin->GetOrAddReal("problem", "vx", 0);
  const Real s_init = pin->GetOrAddReal("problem", "s_init", 0.);
  const Real pres = rho*SQR(iso_cs);
  const Real gm1  = peos->GetGamma() - 1.0;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        //density
        phydro->u(IDN, k, j, i) = rho;
        //velocity, x direction
        phydro->u(IM1, k, j, i) = rho*vx;
        //energy
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN, k, j, i) = pres/gm1 + 0.5*rho*SQR(vx);
        }
      }
    }
  }

  //intialize isotopic abundances
  if (NSCALARS > 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ispec=0; ispec < NSCALARS; ++ispec) {
            pscalars->s(ispec, k, j, i) = s_init*rho;
#ifdef INCLUDE_CHEMISTRY
            Real s_ispec = pin->GetOrAddReal("problem",
                "s_init_"+pscalars->chemnet.species_names[ispec], -1);
            // Real s_ispec = 1.;
            if (s_ispec >= 0.) {
              pscalars->s(ispec, k, j, i) = s_ispec*rho;
            }
#endif
          }
        }
      }
    }
  }

  return;
}

Real ReactionTimeStep(MeshBlock *pmb) {
    Real min_dt=FLT_MAX;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real dt;
        dt = 2*pmb->pscalars->h(k,j,i) + 1e-20;
        min_dt = std::min(min_dt, dt);
        // printf("min_dt = %.2e\n", min_dt);
      }
    }
  }
  return min_dt;
}
