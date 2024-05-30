//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeFiniteStrainElasticStress.h"

#include "StressUpdateCPBase.h"
#include "ComputeCrystalPlasticityEigenstrainBase.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"

/**
 * Multiplicative decomposition based stress update algorithm for crystal plasticity
 */
class StressUpdateCP : public ComputeFiniteStrainElasticStress
{
public:
  static InputParameters validParams();

  StressUpdateCP(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void computeQpStress() override;

  /**
   * Updates the stress (PK2)
   */
  virtual void updateStress(RankTwoTensor & cauchy_stress, RankFourTensor & jacobian_mult);

  virtual void initQpStatefulProperties() override;

  /**
   * Residual and jacobian for stress update
   */
  void calculateResidualAndJacobian();


  void preSolveQp();
  void solveQp();
  void postSolveQp();

  /**
   * Update state variables
   */
  void solveStateVariables();

  void solveStress();

  /**
   * Stress residual 
   */
  void calculateResidual();

  /**
   * @brief 
   * 
   * @return * Stress 
   */
  void calculateJacobian();

  void calcTangentModuli(RankFourTensor & jacobian_mult);
  void elasticTangentModuli(RankFourTensor & jacobian_mult);
  void elastoPlasticTangentModuli(RankFourTensor & jacobian_mult);
  ///@}

  /// performs the line search update
  bool lineSearchUpdate(const Real & rnorm_prev, const RankTwoTensor & dpk2);

  // compute some state variables for postprocessing
	virtual void computeOtherQuantities();
  /**
   * Calculates the deformation gradient due to eigenstrain
   */
  void calculateEigenstrainDeformationGrad();

  /// number of plastic models
  const unsigned _num_models;

  /// The user supplied cyrstal plasticity consititutive models
  std::vector<StressUpdateCPBase *> _models;

  /// number of eigenstrains
  const unsigned _num_eigenstrains;

  /// The user supplied cyrstal plasticity eigenstrains
  std::vector<ComputeCrystalPlasticityEigenstrainBase *> _eigenstrains;

  /// optional parameter to define several mechanical systems on the same block, e.g. multiple phases
  const std::string _base_name;

  /// Elasticity tensor as defined by a separate class
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  /// Stress residual equation relative tolerance
  Real _rtol = 1e-2;
  /// Stress residual equation absolute tolerance
  Real _abs_tol = 1e-6;

  /// Residual tensor
  RankTwoTensor _residual_tensor;
  /// Jacobian tensor
  RankFourTensor _jacobian;

  /// iteration number
  unsigned int _maxiter=100;
  unsigned int _maxiterg = 100;

  /// Maximum number of substep iterations
  unsigned int _max_substep_iter;
  const enum class TangentModuliType { EXACT, NONE } _tan_mod_type;
  
  /// time step size during substepping
  Real _substep_dt;

  /// Flag to activate line serach
  bool _use_line_search;

  /// Minimum line search step size
  Real _min_line_search_step_size;

  /// Line search bisection method tolerance
  Real _line_search_tolerance;

  /// Line search bisection method maximum iteration number
  unsigned int _line_search_max_iterations;

  /// strain formulation
  const enum class LineSearchMethod { CutHalf, Bisection } _line_search_method;

  ///@{Plastic deformation gradient RankTwoTensor for the crystal
  MaterialProperty<RankTwoTensor> & _plastic_deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _plastic_deformation_gradient_old;
  ///@}

  ///@{ Generalized eigenstrain deformation gradient RankTwoTensor for the crystal
  MaterialProperty<RankTwoTensor> * _eigenstrain_deformation_gradient;
  const MaterialProperty<RankTwoTensor> * _eigenstrain_deformation_gradient_old;
  ///@}

  ///@{Total deformation gradient RankTwoTensor for the crystal
  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  ///@}

  ///@{Second Piola-Kirchoff stress measure
  MaterialProperty<RankTwoTensor> & _pk2;
  const MaterialProperty<RankTwoTensor> & _pk2_old;
  ///@}

  /// Lagrangian total strain measure for the entire crystal
  MaterialProperty<RankTwoTensor> & _total_lagrangian_strain;

  /**
   * Tracks the rotation of the crystal during deformation
   * Note: this rotation tensor is not applied to the crystal lattice
   */
  MaterialProperty<RankTwoTensor> & _updated_rotation;

  /**
   * Crystal rotation in the original, or reference, configuration as defined by
   * Euler angle arguments in the ComputeElasticityTensor classes
   */
  const MaterialProperty<RankTwoTensor> & _crysrot;

  ///@{Helper deformation gradient tensor variables used in iterative solve
  RankTwoTensor _temporary_deformation_gradient;
  RankTwoTensor _elastic_deformation_gradient;
  RankTwoTensor _inverse_plastic_deformation_grad;
  RankTwoTensor _inverse_plastic_deformation_grad_old;
  RankTwoTensor _inverse_eigenstrain_deformation_grad;
  ///@}

  /// Flag to print to console warning messages on stress, constitutive model convergence
  const bool _print_convergence_message;
  bool _print_debug_message = false;
  const unsigned int target_element = 0;

  /// Flag to check whether convergence is achieved or if substepping is needed
  bool _convergence_failed;

  ///@{ Used for substepping; Uniformly divides the increment in deformation gradient
  RankTwoTensor _delta_deformation_gradient;
  RankTwoTensor _temporary_deformation_gradient_old;
  ///@}

  /// Scales the substepping increment to obtain deformation gradient at a substep iteration
  Real _dfgrd_scale_factor;

  MaterialProperty<RankTwoTensor> & _plastic_strain;
  MaterialProperty<Real> & _plastic_strain_Eff;
  MaterialProperty<RankTwoTensor> & _stress_deviatoric;
  MaterialProperty<Real> & _stress_VM;
};
