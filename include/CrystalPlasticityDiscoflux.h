// Subhendu Chakraborty, Abigail Hunter, Darby J. Luscher
// Funding: LDRD project XX9A, XXNA
/*----------
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for 
Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC 
for the U.S. Department of Energy/National Nuclear Security Administration. 
All rights in the program are reserved by Triad National Security, LLC, and the 
U.S. Department of Energy/National Nuclear Security Administration. 
The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, 
irrevocable worldwide license in this material to reproduce, prepare. derivative works, 
distribute copies to the public, perform publicly and display publicly, and to 
permit others to do so.
------------*/

#include "MooseVariableFE.h"
#include "MooseVariableScalar.h"
#include "MooseVariableInterface.h"
#include "StressUpdateCPBase.h"


class CrystalPlasticityDiscoflux;


class CrystalPlasticityDiscoflux : public StressUpdateCPBase
{
public:
  static InputParameters validParams();

  CrystalPlasticityDiscoflux(const InputParameters & parameters);

protected:
  // initializes the stateful properties
  virtual void initQpStatefulProperties() override;

  virtual void setInitialConstitutiveVariableValues() override;

  virtual void setSubstepConstitutiveVariableValues() override;

  /**
   * Stores the current value of the slip system resistance into a separate
   * material property in case substepping is needed.
   */
  virtual void updateSubstepConstitutiveVariableValues() override;

  virtual bool calculateSlipRate() override;

  virtual void calculateConstitutiveSlipDerivative(std::vector<Real> & dslip_dtau) override;

  // Cache the slip system value before the update for the diff in the convergence check
  virtual void cacheStateVariablesBeforeUpdate() override;

  /**
   * Following the Constitutive model for slip system resistance as given in:
   */
  virtual void calculateStateVariableEvolutionRateComponent() override;

  /*
   * Finalizes the values of the state variables and slip system resistance
   * for the current timestep after convergence has been reached.
   */
  virtual bool updateStateVariables() override;

  /*
   * Determines if the state variables, e.g. defect densities, have converged
   * by comparing the change in the values over the iteration period.
   */
  virtual bool areConstitutiveStateVariablesConverged() override;
  
  
  virtual void getDDIncrements();
  
    /**
  * This function
  * stores the dislocation velocity value
  * to couple with dislocation transport
  */
  virtual void getDisloVelocity();

  
  const Real _lattice_friction;
  const Real _burgers_vector_mag;
  const Real _dislo_density_initial;
  const Real _dislo_density_factor_CDT;
  const Real _C_multi, _C_trap, _C_m_ann, _C_im_ann;
  Real _Coeff_hardening, _q1, _q2, _c1, _temp;
  Real Lbar;
  
  //const VariableValue & _Rho_EdgePositive_01;
  const ArrayVariableValue & _DD_EdgePositive;
  const ArrayVariableValue & _DD_EdgeNegative;
  const ArrayVariableValue & _DD_ScrewPositive;
  const ArrayVariableValue & _DD_ScrewNegative;
  
  const ArrayVariableGradient & _DD_EdgePositive_Grad;
  const ArrayVariableGradient & _DD_EdgeNegative_Grad;
  const ArrayVariableGradient & _DD_ScrewPositive_Grad;
  const ArrayVariableGradient & _DD_ScrewNegative_Grad;
  
  MaterialProperty<std::vector<Real>> & _dislocation_immobile;
  const MaterialProperty<std::vector<Real>> & _dislocation_immobile_old;
  MaterialProperty<std::vector<Real>> & _dislo_velocity_edge;
  MaterialProperty<std::vector<Real>> & _dislo_velocity_screw;
  const MaterialProperty<std::vector<Real>> & _tau_old;
  const MaterialProperty<std::vector<Real>> & _GND_density;
  MaterialProperty<std::vector<Real>> & _kappa;
  
  // DDC related variables
  std::vector<RealVectorValue> _DD_grad;
  std::vector<Real> _tau_b_local;

  // Stores the values of the slip system resistance 
  std::vector<Real> _previous_substep_slip_resistance;
  std::vector<Real> _previous_substep_dislocation_mobile;
  std::vector<Real> _previous_substep_dislocation_immobile;
  std::vector<Real> _previous_substep_slip_accumulated;

  //Caches the value of the current slip system resistance
  std::vector<Real> _slip_resistance_before_update;
  std::vector<Real> _dislocation_mobile_before_update;
  std::vector<Real> _dislocation_immobile_before_update;

  RealVectorValue slip_direction_rotated, slip_plane_normal_rotated;
  RankTwoTensor Stress_internal;
};
