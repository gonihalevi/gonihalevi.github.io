#ifndef ALPHA13_HPP
#define ALPHA13_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file H2.hpp
//  \brief definitions for a very simple chemical network with H2 formation on grains,
//  and H2 distruction by CR. This has an analytic solution.
//======================================================================================
//c++ headers
#include <string> //std::string

// Athena++ classes headers
#include "network.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

//! \class ChemNetwork
//  \brief Chemical Network that defines the reaction rates between species.
//  Note: This is a template for chemistry network.
//  When implementing a new chemistry network, all public functions should be
//  in the same form.
//  The internal calculations are in cgs units. The input and 
//  return of RHS and Edot must be in code units.
class ChemNetwork : public NetworkWrapper {
  //It would be convenient to know the species names in
  //initialization of chemical species in problem
  friend class MeshBlock; 
public:
  ChemNetwork(MeshBlock *pmb, ParameterInput *pin);
  ~ChemNetwork();

	//a list of species name, used in output
	static const std::string species_names[NSCALARS];

	//Set the rates of chemical reactions, eg. through density and radiation field.
  //k, j, i are the corresponding index of the grid
  void InitializeNextStep(const int k, const int j, const int i);

  //RHS: right-hand-side of ODE. dy/dt = ydot(t, y). Here y are the abundance
  //of species. details see CVODE package documentation.
  //all input/output variables are in code units
  void RHS(const Real t, const Real y[NSCALARS], const Real ED,
           Real ydot[NSCALARS]);
  
  //energy equation dED/dt, all input/output variables are in code units
  //(ED is the energy density)
  Real Edot(const Real t, const Real y[NSCALARS], const Real ED);

private:
  PassiveScalars *pmy_spec_;
	MeshBlock *pmy_mb_;

	std::string species_names_all_[NSCALARS];//all species
	//index of species
	static const int iHe_;
	static const int iC_;
  static const int iO_;
  static const int iNe_;
  static const int iMg_;
  static const int iSi_;
  static const int iS_;
  static const int iAr_;
  static const int iCa_;
  static const int iTi_;
  static const int iCr_;
  static const int iFe_;
  static const int iNi_;

  static const int NISO = 13;
  static const int NEQN = 14;
  static const int NREAC = 18;
  static const int NALP = 38;

  // isotope atomic weights
  static constexpr Real Aiso[NISO] =
  {4.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0, 44.0, 48.0, 52.0, 56.0};
  // isotope atomic numbers
  static constexpr Real Ziso[NISO] =
  {2.0,  6.0,  8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0};
  static constexpr Real alphanet13_Tcold = 2.e8; //Temp cutoff, below which plasma is assumed to be inert

	// //units 
  Real gm1_; //adiabatic gamm, read from input
	Real unit_density; //read from input
  int NISOfuel; //read from input
  Real alphanet_epsder;
	Real unit_length_in_cm_; //read from input
	Real unit_vel_in_cms_; //read from input
	Real unit_time_in_s_; //from length and velocity units
  //unit of energy density, in erg cm-3, from density and velocity units
  Real unit_E_in_cgs_; 
	Real rho_; //density, updated at InitializeNextStep from hydro variable

  /*-----------------------------------------------------------------------------
   * Calculate reaction rates
   *
   * Input:
   *     rho     - density (g/cm3)
   *     tp      - temperature (K)
   *
   * Output:
   *     frv[18] - forward reaction rates
   *     rev[18] - backward reaction rate coefficients
   *               (backward rate = frv * rev)
   *-----------------------------------------------------------------------------*/
  void CalculateRates(Real rho, Real tp, Real frv[NREAC], Real rev[NREAC]);

  /*-----------------------------------------------------------------------------
   * Calculate right hand sides of nuclear kinetic equations, and the energy
   * generation rate
   *
   * Input:
   *     frv[18] - forward reaction rates
   *     rev[18] - backward reaction rate coefficients
   *     y[14]   - alpha-nuclei mole fractions
   *                 y[0]  - 4He
   *                 y[1]  - 12C
   *                 ..........
   *                 ..........
   *                 y[12] - 56Ni
   *                 y[13] - scaled energy (not used here)
   *
   * Output:
   *     f[14]   - rates of change
   *                 f[0:12] - dy(i)/dt [1/sec]
   *                 f[13])  - de/dt    [ergs/gm/sec]
   *-----------------------------------------------------------------------------*/
  void RatesOfChange(const Real frv[NREAC], const Real rev[NREAC],
      const Real y[NEQN], Real f[NEQN]);

  /*-----------------------------------------------------------------------------
   * Calculate right hand sides, energy generation rate and their derivatives
   * with respect to mole fractions of reactants
   *
   * Input:
   *     frv[18]  - forward reaction rates
   *     rev[18]  - backward reaction rate coefficients
   *     y[14]    - alpha-nuclei mole fractions
   *                  y[0] - 4He
   *                  y[1] - 12C
   *                  ..........
   *                  ..........
   *                  y[12] - 32S
   *                  y[13] - scaled energy (not used here)
   *
   * Output:
   *     f[14]    - rates of change
   *                  f[0:12] - dy(i)/dt [1/sec]
   *                  f[13])  - de/dt    [ergs/gm/sec]
   *     df[14][14] - partial derivatives of f with respect to mole fractions
   *-----------------------------------------------------------------------------*/
  void PartialDerivatives(const Real frv[NREAC], const Real rev[NREAC],
           const Real y[NEQN], Real f[NEQN], Real df[NEQN][NEQN]);

  /*-----------------------------------------------------------------------------
   * Calculates production rates f(i) = dy(i)/dt along with the energy generation
   * rate dE/dt, and the Jacobian matrix of the right hand sides for YASS.
   *
   * Input:
   *     y[14]       - initial solution: isotope mole fractions y[0:12],
   *                   internal energy density scaled by e0, y[13] (= 1.0 on input)
   *     rdata[2]    - density, g/cm3 (=const)
   *                   e0, energy scale, erg/g (=const)
   *
   * Output:
   *     f[14]       - f[i] = dy[i]/dt, [1/sec]
   *     jac[14][14] - jac[j][i] = df[i]/dy[j]
   *
   * Return flag:
   *     0           - do not integrate this cell, it is too cold
   *     1           - integrate this cell
   *
   *  Global parameters:
   *     alphanet_epsder - increment of energy epsder*y(13) is used to calculate
   *                       df(i)/dy(13) numerically
   *-----------------------------------------------------------------------------*/
  int RHSFull(const Real y[NEQN], Real f[NEQN], Real jac[NEQN][NEQN], Real * rdata);

  //calculate Jacobian with numerical differentiation 
  void Jacobian_numerical(const Real t, const Real y[NSCALARS],
                          const Real ydot[NSCALARS], const Real ED,
                          AthenaArray<Real> &jac);
};

#endif // ALPHA13_HPP
