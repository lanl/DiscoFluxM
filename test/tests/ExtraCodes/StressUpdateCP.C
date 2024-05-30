//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StressUpdateCP.h"

#include "StressUpdateCPBase.h"
#include "libmesh/utility.h"
#include "Conversion.h"
#include "MooseException.h"

registerMooseObject("TensorMechanicsApp", StressUpdateCP);

InputParameters
StressUpdateCP::validParams()
{
  InputParameters params = ComputeFiniteStrainElasticStress::validParams();

  params.addClassDescription("Crystal Plasticity class according to the DiscoFlux material model");
  params.addParam<std::string>("base_name","customize for multiple material system");
  params.addRequiredParam<std::vector<MaterialName>>("discoflux_model_name", "Name of the DiscoFlux material model.");
  params.addParam<std::vector<MaterialName>>("eigenstrain_names", {}, "The material objects to calculate eigenstrains.");
  params.addParam<unsigned int>("maximum_substep_iteration", 1, "Maximum number of substep iteration");
  params.addParam<MooseEnum>("tan_mod_type",MooseEnum("exact none", "exact"),"Tangent modulus");
  params.addParam<bool>("use_line_search", false, "Use line search in constitutive update");
  params.addParam<Real>("min_line_search_step_size", 0.01, "Minimum line search step size");
  params.addParam<Real>("line_search_tol", 0.5, "Line search bisection method tolerance");
  params.addParam<unsigned int>("line_search_maxiter", 20, "Line search bisection method maximum number of iteration");
  params.addParam<MooseEnum>("line_search_method",MooseEnum("CUT_HALF BISECTION", "CUT_HALF"),"The method used in line search");
  params.addParam<bool>(
      "print_state_variable_convergence_error_messages",
      false,
      "Whether or not to print warning messages from the crystal plasticity specific convergence "
      "checks on the stress measure and general constitutive model quantinties.");
  return params;
}

StressUpdateCP::StressUpdateCP(
    const InputParameters & parameters)
  : ComputeFiniteStrainElasticStress(parameters),
    _num_models(1),
    _num_eigenstrains(getParam<std::vector<MaterialName>>("eigenstrain_names").size()),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") : ""),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _max_substep_iter(getParam<unsigned int>("maximum_substep_iteration")),
    _tan_mod_type(getParam<MooseEnum>("tan_mod_type").getEnum<TangentModuliType>()),
    _use_line_search(getParam<bool>("use_line_search")),
    _min_line_search_step_size(getParam<Real>("min_line_search_step_size")),
    _line_search_tolerance(getParam<Real>("line_search_tol")),
    _line_search_max_iterations(getParam<unsigned int>("line_search_maxiter")),
    _line_search_method(getParam<MooseEnum>("line_search_method").getEnum<LineSearchMethod>()),
    _plastic_deformation_gradient(declareProperty<RankTwoTensor>("plastic_deformation_gradient")),
    _plastic_deformation_gradient_old(
        getMaterialPropertyOld<RankTwoTensor>("plastic_deformation_gradient")),
    _eigenstrain_deformation_gradient(
        _num_eigenstrains ? &declareProperty<RankTwoTensor>("eigenstrain_deformation_gradient")
                          : nullptr),
    _eigenstrain_deformation_gradient_old(
        _num_eigenstrains
            ? &getMaterialPropertyOld<RankTwoTensor>("eigenstrain_deformation_gradient")
            : nullptr),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>(_base_name + "deformation_gradient")),
    _deformation_gradient_old(
        getMaterialPropertyOld<RankTwoTensor>(_base_name + "deformation_gradient")),
    _pk2(declareProperty<RankTwoTensor>("second_piola_kirchhoff_stress")),
    _pk2_old(getMaterialPropertyOld<RankTwoTensor>("second_piola_kirchhoff_stress")),
    _total_lagrangian_strain(
        declareProperty<RankTwoTensor>("total_lagrangian_strain")),
    _updated_rotation(declareProperty<RankTwoTensor>("updated_rotation")),
    _crysrot(getMaterialProperty<RankTwoTensor>(_base_name + "crysrot")), 
    _print_convergence_message(getParam<bool>("print_state_variable_convergence_error_messages")),


  _plastic_strain(declareProperty<RankTwoTensor>(_base_name + "plastic_strain")),
	_plastic_strain_Eff(declareProperty<Real>(_base_name + "plastic_strain_Eff")),
	_stress_deviatoric(declareProperty<RankTwoTensor>(_base_name + "stress_deviatoric")),
	_stress_VM(declareProperty<Real>(_base_name + "stress_VM"))
{
  _convergence_failed = false;
}

void
StressUpdateCP::initQpStatefulProperties()
{
  _plastic_deformation_gradient[_qp].zero();
  _plastic_deformation_gradient[_qp].addIa(1.0);

  if (_num_eigenstrains)
  {
    (*_eigenstrain_deformation_gradient)[_qp].zero();
    (*_eigenstrain_deformation_gradient)[_qp].addIa(1.0);
  }
  else
  {
    // set to identity if no eigenstrain is added
    _inverse_eigenstrain_deformation_grad.zero();
    _inverse_eigenstrain_deformation_grad.addIa(1.0);
  }

  _pk2[_qp].zero();

  _total_lagrangian_strain[_qp].zero();

  _updated_rotation[_qp].zero();
  _updated_rotation[_qp].addIa(1.0);

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->setQp(_qp);
    _models[i]->setMaterialVectorSize();
    _models[i]->initQpStatefulProperties();
  }

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    _eigenstrains[i]->setQp(_qp);
    _eigenstrains[i]->initQpStatefulProperties();
  }

  _plastic_strain[_qp].zero();
  _plastic_strain_Eff[_qp] = 0.00;
  _stress_deviatoric[_qp].zero();
  _stress_VM[_qp] = 0.00;
}

void
StressUpdateCP::initialSetup()
{
  // get crystal plasticity models
  std::vector<MaterialName> model_names =
      getParam<std::vector<MaterialName>>("discoflux_model_name");

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    StressUpdateCPBase * model =
        dynamic_cast<StressUpdateCPBase *>(&getMaterialByName(model_names[i]));

      _models.push_back(model);
  }

  // get crystal plasticity eigenstrains
  std::vector<MaterialName> eigenstrain_names =
      getParam<std::vector<MaterialName>>("eigenstrain_names");

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    ComputeCrystalPlasticityEigenstrainBase * eigenstrain =
        dynamic_cast<ComputeCrystalPlasticityEigenstrainBase *>(
            &getMaterialByName(eigenstrain_names[i]));

    if (eigenstrain)
      _eigenstrains.push_back(eigenstrain);
    else
      mooseError("Eigenstrain" + eigenstrain_names[i] +
                 " is not compatible with StressUpdateCP");
  }
}

void
StressUpdateCP::computeQpStress()
{
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->setQp(_qp);
  }

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
    _eigenstrains[i]->setQp(_qp);
  updateStress(_stress[_qp], _Jacobian_mult[_qp]); 
}

void
StressUpdateCP::updateStress(RankTwoTensor & cauchy_stress,
                                                     RankFourTensor & jacobian_mult)
{

  // Initialize substepping variables
  unsigned int substep_iter = 1;
  unsigned int num_substep = 1;

  _temporary_deformation_gradient_old = _deformation_gradient_old[_qp];
  if (_temporary_deformation_gradient_old.det() == 0)
    _temporary_deformation_gradient_old.addIa(1.0);

  _delta_deformation_gradient = _deformation_gradient[_qp] - _temporary_deformation_gradient_old;

  for (unsigned int i = 0; i < _num_models; ++i)
    _models[i]->calculateSchmidTensor();

  do
  {
    _convergence_failed = false;
    preSolveQp();

    _substep_dt = _dt / num_substep;

    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->setSubstepDt(_substep_dt);

    // calculate F^{eigen} only when we have eigenstrain
    if (_num_eigenstrains)
      calculateEigenstrainDeformationGrad();

    for (unsigned int istep = 0; istep < num_substep; ++istep)
    {
      _temporary_deformation_gradient =
          (static_cast<Real>(istep) + 1) / num_substep * _delta_deformation_gradient;
      _temporary_deformation_gradient += _temporary_deformation_gradient_old;

      solveQp();

      if (_convergence_failed)
      {
        substep_iter++;
        num_substep *= 2;
        break;
      }
    }

    if (substep_iter > _max_substep_iter && _convergence_failed)
    {
      mooseException("StressUpdateCP: Constitutive failure");
    }
  } while (_convergence_failed);


  cauchy_stress = _elastic_deformation_gradient * _pk2[_qp] *
                  _elastic_deformation_gradient.transpose() / _elastic_deformation_gradient.det();

  calcTangentModuli(jacobian_mult);

  _total_lagrangian_strain[_qp] =
      _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] -
      RankTwoTensor::Identity();
  _total_lagrangian_strain[_qp] = _total_lagrangian_strain[_qp] * 0.5;


  postSolveQp();
}

void
StressUpdateCP::preSolveQp()
{
  for (unsigned int i = 0; i < _num_models; ++i)
    _models[i]->setInitialConstitutiveVariableValues();

  _pk2[_qp] = _pk2_old[_qp];
  _inverse_plastic_deformation_grad_old = _plastic_deformation_gradient_old[_qp].inverse();
}

void
StressUpdateCP::solveQp()
{
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->setSubstepConstitutiveVariableValues();
    _models[i]->calculateSlipResistance();
    _models[i]->coputeStateVariables();
  }

  _inverse_plastic_deformation_grad = _inverse_plastic_deformation_grad_old;

  solveStateVariables();
  if (_convergence_failed)
    return; // pop back up and take a smaller substep

  for (unsigned int i = 0; i < _num_models; ++i)
    _models[i]->updateSubstepConstitutiveVariableValues();

  // save off the old F^{p} inverse now that have converged on the stress and state variables
  _inverse_plastic_deformation_grad_old = _inverse_plastic_deformation_grad;
}

void
StressUpdateCP::postSolveQp()
{
  // Calculate crystal rotation to track separately
  RankTwoTensor rot;
  _elastic_deformation_gradient.getRUDecompositionRotation(rot);
  _updated_rotation[_qp] =  _crysrot[_qp] * rot;

// compute the other quantities
  computeOtherQuantities();
}

void
StressUpdateCP::solveStateVariables()
{
  unsigned int iteration;
  bool iter_flag = true;

  iteration = 0;
  // Check for slip system resistance update tolerance
  do
  {
    solveStress();

    if (_convergence_failed)
      return;

    _plastic_deformation_gradient[_qp] =
        _inverse_plastic_deformation_grad.inverse(); // the postSoveStress

    // Update slip system resistance and state variable after the stress has been finalized
    // We loop through all the models for each calculation
    // in order to make sure that when coupling appears, the state variables are updated based on
    // the same components
    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->cacheStateVariablesBeforeUpdate();

    for (unsigned int i = 0; i < _num_models; ++i)
      _models[i]->calculateStateVariableEvolutionRateComponent();

    for (unsigned int i = 0; i < _num_models; ++i)
      if (!_models[i]->updateStateVariables())
        _convergence_failed = true;

    if (_convergence_failed)
      return;

    for (unsigned int i = 0; i < _num_models; ++i)
    {
      // iter_flag = false, stop iteration if all values are converged and good to go
      // iter_flag = true, continue iteration if any value is not converged
      if (!_models[i]->areConstitutiveStateVariablesConverged())
      {
        iter_flag = true;
        break;
      }
      else
        iter_flag = false; // iter_flag = false, stop iteration only when all models returns true
    }
    iteration++;
  } while (iter_flag && iteration < _maxiterg);

  if (iteration == _maxiterg)
  {
    if (_print_convergence_message)
      mooseWarning(" Convergence failed.... \n");

    _convergence_failed = true;
  }
}

void
StressUpdateCP::solveStress()
{
  unsigned int iteration = 0;
  RankTwoTensor dpk2;
  Real rnorm, rnorm0, rnorm_prev;

  // Calculate stress residual
  calculateResidualAndJacobian();

  if (_convergence_failed)
  {
    return;
  }

  rnorm = _residual_tensor.L2norm();
  rnorm0 = rnorm;
  
  // Check for stress residual tolerance; different from user object version which
  // compares the absolute tolerance of only the original rnorm value
  while (rnorm > _rtol * rnorm0 && rnorm > _abs_tol && iteration < _maxiter)
  {
    
    // Calculate stress increment
    dpk2 = -_jacobian.invSymm() * _residual_tensor;
    _pk2[_qp] = _pk2[_qp] + dpk2;

    calculateResidualAndJacobian();
    if (_convergence_failed)
    {
      if (_print_convergence_message)
        mooseWarning(" Convergence failed.... \n");
      return;
    }

    rnorm_prev = rnorm;
    rnorm = _residual_tensor.L2norm();
    if (_use_line_search && rnorm > rnorm_prev && !lineSearchUpdate(rnorm_prev, dpk2))
    {
      if (_print_convergence_message)
        mooseWarning("Failed with line search.... \n");

      _convergence_failed = true;
      return;
    }

    if (_use_line_search)
      rnorm = _residual_tensor.L2norm();

    iteration++;
  }

  if (iteration >= _maxiter)
  {
    if (_print_convergence_message)
      mooseWarning(" Convergence failed.... \n");

    _convergence_failed = true;
  }
}

// Calculates stress residual equation and jacobian
void
StressUpdateCP::calculateResidualAndJacobian()
{
  calculateResidual();
  if (_convergence_failed)
    return;
  calculateJacobian();
}

void
StressUpdateCP::calculateResidual()
{
  RankTwoTensor ce, elastic_strain, ce_pk2, equivalent_slip_increment_per_model,
      equivalent_slip_increment, pk2_new;

  equivalent_slip_increment.zero();

  // calculate slip rate in order to compute F^{p-1}
  for (unsigned int i = 0; i < _num_models; ++i)
  {
    equivalent_slip_increment_per_model.zero();

    // calculate shear stress with consideration of contribution from other physics
    _models[i]->calculateShearStress(
        _pk2[_qp], _inverse_eigenstrain_deformation_grad, _num_eigenstrains);

    _convergence_failed = !_models[i]->calculateSlipRate();

    if (_convergence_failed)
    {
      return;
    }

    _models[i]->calculateEquivalentSlipIncrement(equivalent_slip_increment_per_model);
    equivalent_slip_increment += equivalent_slip_increment_per_model;
  }

  RankTwoTensor residual_equivalent_slip_increment =
      RankTwoTensor::Identity() - equivalent_slip_increment;
  _inverse_plastic_deformation_grad =
      _inverse_plastic_deformation_grad_old * residual_equivalent_slip_increment;

  _elastic_deformation_gradient = _temporary_deformation_gradient *
                                  _inverse_eigenstrain_deformation_grad *
                                  _inverse_plastic_deformation_grad;

  ce = _elastic_deformation_gradient.transpose() * _elastic_deformation_gradient;

  elastic_strain = ce - RankTwoTensor::Identity();
  elastic_strain *= 0.5;

  pk2_new = _elasticity_tensor[_qp] * elastic_strain;
  _residual_tensor = _pk2[_qp] - pk2_new;
}

void
StressUpdateCP::calculateJacobian()
{
  // may not need to cache the dfpinvdpk2 here. need to double check
  RankFourTensor dfedfpinv, deedfe, dfpinvdpk2, dfpinvdpk2_per_model;

  RankTwoTensor ffeiginv = _temporary_deformation_gradient * _inverse_eigenstrain_deformation_grad;

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
        dfedfpinv(i, j, k, j) = ffeiginv(i, k);

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _elastic_deformation_gradient(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _elastic_deformation_gradient(k, i) * 0.5;
      }

  for (unsigned int i = 0; i < _num_models; ++i)
  {
    _models[i]->calculateTotalPlasticDeformationGradientDerivative(
        dfpinvdpk2_per_model,
        _inverse_plastic_deformation_grad_old,
        _inverse_eigenstrain_deformation_grad);
    dfpinvdpk2 += dfpinvdpk2_per_model;
  }

  _jacobian =
      RankFourTensor::IdentityFour() - (_elasticity_tensor[_qp] * deedfe * dfedfpinv * dfpinvdpk2);
}

void
StressUpdateCP::elasticTangentModuli(RankFourTensor & jacobian_mult)
{
  // update jacobian_mult
  jacobian_mult = _elasticity_tensor[_qp];
}

void
StressUpdateCP::calcTangentModuli(RankFourTensor & jacobian_mult)
{
  switch (_tan_mod_type)
  {
    case TangentModuliType::EXACT:
      elastoPlasticTangentModuli(jacobian_mult); // has some convergence issue
      break;
    default:
      elasticTangentModuli(jacobian_mult);
  }
}

void
StressUpdateCP::elastoPlasticTangentModuli(RankFourTensor & jacobian_mult)
{
  RankFourTensor tan_mod;
  RankTwoTensor pk2fet, fepk2, feiginvfpinv;
  RankFourTensor deedfe, dsigdpk2dfe, dfedf;

  // Fill in the matrix stiffness material property
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _elastic_deformation_gradient(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _elastic_deformation_gradient(k, i) * 0.5;
      }

  usingTensorIndices(i_, j_, k_, l_);
  dsigdpk2dfe = _elastic_deformation_gradient.times<i_, k_, j_, l_>(_elastic_deformation_gradient) *
                _elasticity_tensor[_qp] * deedfe;

  pk2fet = _pk2[_qp] * _elastic_deformation_gradient.transpose();
  fepk2 = _elastic_deformation_gradient * _pk2[_qp];

  for (const auto l : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto i : make_range(Moose::dim))
      {
        tan_mod(i, j, i, l) += pk2fet(l, i);
        tan_mod(i, j, j, l) += fepk2(j, l);
      }

  tan_mod += dsigdpk2dfe;

  const auto je = _elastic_deformation_gradient.det();
  if (je > 0.0)
    tan_mod /= je;

  feiginvfpinv = _inverse_eigenstrain_deformation_grad * _inverse_plastic_deformation_grad;
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto l : make_range(Moose::dim))
        dfedf(i, j, i, l) = feiginvfpinv(l, j);

  jacobian_mult = tan_mod * dfedf;
}

bool
StressUpdateCP::lineSearchUpdate(const Real & rnorm_prev,const RankTwoTensor & dpk2)
{
  if (_line_search_method == LineSearchMethod::CutHalf)
  {
    Real rnorm;
    Real step = 1.0;

    do
    {
      _pk2[_qp] = _pk2[_qp] - step * dpk2;
      step /= 2.0;
      _pk2[_qp] = _pk2[_qp] + step * dpk2;

      calculateResidual();
      rnorm = _residual_tensor.L2norm();
    } while (rnorm > rnorm_prev && step > _min_line_search_step_size);

    // has norm improved or is the step still above minumum search step size?
    return (rnorm <= rnorm_prev || step > _min_line_search_step_size);
  }
  else if (_line_search_method == LineSearchMethod::Bisection)
  {
    unsigned int count = 0;
    Real step_a = 0.0;
    Real step_b = 1.0;
    Real step = 1.0;
    Real s_m = 1000.0;
    Real rnorm = 1000.0;

    calculateResidual();
    auto s_b = _residual_tensor.doubleContraction(dpk2);
    const auto rnorm1 = _residual_tensor.L2norm();
    _pk2[_qp] = _pk2[_qp] - dpk2;
    calculateResidual();
    auto s_a = _residual_tensor.doubleContraction(dpk2);
    const auto rnorm0 = _residual_tensor.L2norm();
    _pk2[_qp] = _pk2[_qp] + dpk2;

    if ((rnorm1 / rnorm0) < _line_search_tolerance || s_a * s_b > 0)
    {
      calculateResidual();
      return true;
    }

    while ((rnorm / rnorm0) > _line_search_tolerance && count < _line_search_max_iterations)
    {
      _pk2[_qp] = _pk2[_qp] - step * dpk2;
      step = 0.5 * (step_b + step_a);
      _pk2[_qp] = _pk2[_qp] + step * dpk2;
      calculateResidual();
      s_m = _residual_tensor.doubleContraction(dpk2);
      rnorm = _residual_tensor.L2norm();

      if (s_m * s_a < 0.0)
      {
        step_b = step;
        s_b = s_m;
      }
      if (s_m * s_b < 0.0)
      {
        step_a = step;
        s_a = s_m;
      }
      count++;
    }

    // below tolerance and max iterations?
    return ((rnorm / rnorm0) < _line_search_tolerance && count < _line_search_max_iterations);
  }
  else
    mooseError("Line search method is not provided.");
}

void
StressUpdateCP::computeOtherQuantities()
{
  // This is to compute few extra state variables once solution is converged

  // Calculate equivalent plastic strain
    _plastic_strain[_qp] = 
      (0.5)*((_plastic_deformation_gradient[_qp].transpose() * _plastic_deformation_gradient[_qp]) - RankTwoTensor::Identity());
  _plastic_strain_Eff[_qp] = _plastic_strain[_qp].L2norm() * std::sqrt(2.0/3); 

  // Compute Von-Misses stress
  // Need to use RankTwoScalarTools::vonMisesStress(stress)
  _stress_deviatoric[_qp] = _stress[_qp].deviatoric(); // change it to local variable of rank2tensor
  _stress_VM[_qp] =  _stress_deviatoric[_qp].L2norm() * std::sqrt(3.0/2); 
  
}


void
StressUpdateCP::calculateEigenstrainDeformationGrad()
{
  _inverse_eigenstrain_deformation_grad.zero();
  _inverse_eigenstrain_deformation_grad.addIa(1.0);

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    _eigenstrains[i]->setSubstepDt(_substep_dt);
    _eigenstrains[i]->computeQpProperties();
    _inverse_eigenstrain_deformation_grad *= _eigenstrains[i]->getDeformationGradientInverse();
  }
  (*_eigenstrain_deformation_gradient)[_qp] = _inverse_eigenstrain_deformation_grad.inverse();
}
