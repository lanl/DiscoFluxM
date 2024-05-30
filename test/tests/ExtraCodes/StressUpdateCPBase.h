//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "DelimitedFileReader.h"

class StressUpdateCPBase : public Material
{
public:
  static InputParameters validParams();

  StressUpdateCPBase(const InputParameters & parameters);

  /// Sets the value of the global variable _qp for inheriting classes
  void setQp(const unsigned int & qp);

  /// Sets the value of the _substep_dt for inheriting classes
  void setSubstepDt(const Real & substep_dt);

  ///@{ Retained as empty methods to avoid a warning from Material.C in framework. These methods are unused in all inheriting classes and should not be overwritten.
  virtual void resetQpProperties() final {}
  virtual void resetProperties() final {}
  ///@}

  virtual void initQpStatefulProperties() override;
  virtual void setMaterialVectorSize() ;

  /**
   * Slip-system information
   */
  virtual void getSlipSystems();
  virtual void storeDislocationMobilityInformation();
  void calculateSchmidTensor();

  /**
   * Computes the shear stess for each slip system
   */
  void calculateShearStress(const RankTwoTensor & pk2,
                            const RankTwoTensor & inverse_eigenstrain_deformation_grad,
                            const unsigned int & num_eigenstrains);


  virtual void calculateTotalPlasticDeformationGradientDerivative(
      RankFourTensor & dfpinvdpk2,
      const RankTwoTensor & inverse_plastic_deformation_grad_old,
      const RankTwoTensor & inverse_eigenstrain_deformation_grad_old);

  void calculateSchmidTensor(const unsigned int & number_dislocation_systems,
                             const std::vector<RealVectorValue> & plane_normal_vector,
                             const std::vector<RealVectorValue> & direction_vector,
                             std::vector<RankTwoTensor> & schmid_tensor);

  void sortCrossSlipFamilies();

  unsigned int identifyCrossSlipFamily(const unsigned int index);

  virtual void setInitialConstitutiveVariableValues() {}

  virtual void setSubstepConstitutiveVariableValues() {}

  virtual void updateSubstepConstitutiveVariableValues() {}

  virtual bool calculateSlipRate() = 0;

  virtual void calculateEquivalentSlipIncrement(RankTwoTensor & /*equivalent_slip_increment*/);

  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & /*dslip_dtau*/) = 0;

  virtual void cacheStateVariablesBeforeUpdate() {}

  virtual void calculateStateVariableEvolutionRateComponent() {}

  virtual bool updateStateVariables() = 0;

  virtual void calculateSlipResistance() {}

  virtual void coputeStateVariables() {}

  virtual bool areConstitutiveStateVariablesConverged() { return true; }


  virtual bool isConstitutiveStateVariableConverged(const std::vector<Real> & current_var,
                                                    const std::vector<Real> & var_before_update,
                                                    const std::vector<Real> & previous_substep_var,
                                                    const Real & tolerance);

protected:
  /// Base name prepended to all material property names to allow for
  /// multi-material systems
  const std::string _base_name;

  const enum class CrystalLatticeType { BCC, FCC, HCP } _crystal_lattice_type;

  const std::vector<Real> _unit_cell_dimension;

  ///Maximum number of active slip systems for the crystalline material being modeled
  const unsigned int _number_slip_systems = 12;

  /// File should contain slip plane normal and direction.
  std::string _slip_sys_file_name;

  /// @{Parameters to characterize the cross slip behavior of the crystal
  const Real _number_cross_slip_directions;
  const Real _number_cross_slip_planes;
  ///@}

  /// Internal variable update equation tolerance
  Real _rel_state_var_tol;
  /// Slip increment tolerance
  Real _slip_incr_tol;
  /// Tolerance for change in slip system resistance over an increment
  Real _resistance_tol;
  /// Residual tolerance when variable value is zero. Default 1e-12.
  Real _zero_tol;

  ///@{Slip system resistance
  MaterialProperty<std::vector<Real>> & _slip_resistance;
  const MaterialProperty<std::vector<Real>> & _slip_resistance_old;
  ///@}

  /// Current slip increment material property
  MaterialProperty<std::vector<Real>> & _slip_rate;
  MaterialProperty<std::vector<Real>> & _slip_increment;

  ///@{Slip system direction and normal and associated Schmid tensors
  std::vector<RealVectorValue> _slip_direction;
  std::vector<RealVectorValue> _slip_plane_normal;
  MaterialProperty<std::vector<RankTwoTensor>> & _flow_direction;
  ///@}

	MaterialProperty<std::vector<RealVectorValue>> & _slip_direction_edge;
	MaterialProperty<std::vector<RealVectorValue>> & _slip_direction_screw;
	MaterialProperty<std::vector<RealVectorValue>> & _slip_plane_normalboth;

  /// Resolved shear stress on each slip system
  MaterialProperty<std::vector<Real>> & _tau;
  const MaterialProperty<RankTwoTensor> & _crysrot;

  /// Flag to print to console warning messages on stress, constitutive model convergence
  const bool _print_convergence_message;

  /// Substepping time step value used within the inheriting constitutive models
  Real _substep_dt;

  /// Sorted slip system indices into cross slip family groups
  std::vector<std::vector<unsigned int>> _cross_slip_familes;

  /// Flag to run the cross slip calculations if cross slip numbers are specified
  bool _calculate_cross_slip;
  
  MaterialProperty<std::vector<Real>> & _dislocation_mobile;
  MaterialProperty<std::vector<Real>> & _dislocation_immobile;
  const MaterialProperty<std::vector<Real>> & _dislocation_immobile_old;
  MaterialProperty<std::vector<Real>> & _tau_b;
  MaterialProperty<std::vector<Real>> & _kappa;

// Discoflux related material parameters that are constant
//const      Real c1 =2.0;
const  Real mu = 76e+03;
const  Real nu =0.3;
const  Real rho_m = 8960;
const  Real B0 = 3.0e-11;
const  Real g0 = 0.87;
//const  Real q1 = 0.23;
//const  Real q2 = 1.96;
const  Real boltz = 1.38e-23; // Boltzman constant in Jule/Kelvin
//const  Real temp =300;
const  Real omega0 = 2.0e+2; //8.0e+11;

std::vector<Real> _hb;
std::vector<Real> _slip_resistance_increment;  
std::vector<Real> _dislocation_mobile_increment;
std::vector<Real> _dislocation_immobile_increment; 
std::vector<Real> _dv_dtau; 
std::vector<Real> _L_bar;  

//For dislocation velocity computation
Real small2 = 1.0e-10, exp_limit = 2.0e+2;
std::vector<Real> t_wait, t_run, vel_run, dislocation_density, tau_b, xi0, tau_eff, tau_effAbs, tau_effSign, slip_r; 
Real vcrit  = std::sqrt(mu*1.0e+06/rho_m)*1000 ; 
Real deltaG0, inner, deltaG, exp_arg, dtw_dtau, dtr_dtau;
};
