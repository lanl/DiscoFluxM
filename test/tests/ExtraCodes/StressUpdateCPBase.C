//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StressUpdateCPBase.h"

#include "libmesh/utility.h"
#include "libmesh/int_range.h"
#include "Conversion.h"
#include "MooseException.h"

InputParameters
StressUpdateCPBase::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<std::string>(
      "base_name",
      "Optional parameter that allows the user to define multiple crystal plasticity mechanisms");
  params.addClassDescription(
      "Crystal Plasticity base class: handles the Newton iteration over the stress residual and "
      "calculates the Jacobian based on constitutive laws provided by inheriting classes");

  // The return stress increment classes are intended to be iterative materials, so must set compute
  // = false for all inheriting classes
  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");

  params.addParam<MooseEnum>(
      "crystal_lattice_type",
      MooseEnum("BCC FCC HCP", "FCC"),
      "Crystal lattice type or representative unit cell, i.e., BCC, FCC, HCP, etc.");

  params.addRangeCheckedParam<std::vector<Real>>(
      "unit_cell_dimension",
      std::vector<Real>{1.0, 1.0, 1.0},
      "unit_cell_dimension_size = 3",
      "The dimension of the unit cell along three directions, where a cubic unit cell is assumed "
      "for cubic crystals and a hexagonal unit cell (a, a, c) is assumed for HCP crystals. These "
      "dimensions will be taken into account while computing the slip systems."
      " Default size is 1.0 along all three directions.");

  params.addRequiredCoupledVar("DD_EdgePositive", 1.0, "Coupled dislocation density, EdgePositive");
  params.addRequiredCoupledVar("DD_EdgeNegative", 1.0, "Coupled dislocation density, EdgeNegative");
  params.addRequiredCoupledVar("DD_ScrewPositive", 1.0, "Coupled dislocation density, ScrewPositive");
  params.addRequiredCoupledVar("DD_ScrewNegative", 1.0, "Coupled dislocation density, ScrewNegative");

  params.addRequiredParam<unsigned int>(
      "number_slip_systems",
      "The total number of possible active slip systems for the crystalline material");
  params.addRequiredParam<FileName>(
      "slip_sys_file_name",
      "Name of the file containing the slip systems, one slip system per row, with the slip plane "
      "normal given before the slip plane direction.");
  params.addParam<Real>("number_cross_slip_directions",
                        0,
                        "Quanity of unique slip directions, used to determine cross slip familes");
  params.addParam<Real>("number_cross_slip_planes",
                        0,
                        "Quanity of slip planes belonging to a single cross slip direction; used "
                        "to determine cross slip families");
  params.addParam<Real>(
      "slip_increment_tolerance",
      2e-2,
      "Maximum allowable slip in an increment for each individual constitutive model");
  params.addParam<Real>(
      "stol", 1e-2, "Constitutive internal state variable relative change tolerance");
  params.addParam<Real>("resistance_tol",
                        1.0e-2,
                        "Constitutive slip system resistance relative residual tolerance for each "
                        "individual constitutive model");
  params.addParam<Real>("zero_tol",
                        1e-12,
                        "Tolerance for residual check when variable value is zero for each "
                        "individual constitutive model");

  params.addCoupledVar("mobile_DD", 1.0, "mobile dislocations"), // for density based Crystal Plasticity
  params.addCoupledVar("immobile_DD", 1.0, "immobile dislocations"),
  params.addParam<bool>(
      "print_state_variable_convergence_error_messages",
      false,
      "Whether or not to print warning messages from the crystal plasticity specific convergence "
      "checks on both the constiutive model internal state variables.");
  return params;
}

StressUpdateCPBase::StressUpdateCPBase(
    const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") : ""),
    _crystal_lattice_type(
        getParam<MooseEnum>("crystal_lattice_type").getEnum<CrystalLatticeType>()),
    _unit_cell_dimension(getParam<std::vector<Real>>("unit_cell_dimension")),
    _number_slip_systems(getParam<unsigned int>("number_slip_systems")),
    _slip_sys_file_name(getParam<FileName>("slip_sys_file_name")),
    _number_cross_slip_directions(getParam<Real>("number_cross_slip_directions")),
    _number_cross_slip_planes(getParam<Real>("number_cross_slip_planes")),

    _rel_state_var_tol(getParam<Real>("stol")),
    _slip_incr_tol(getParam<Real>("slip_increment_tolerance")),
    _resistance_tol(getParam<Real>("resistance_tol")),
    _zero_tol(getParam<Real>("zero_tol")),

    _slip_resistance(declareProperty<std::vector<Real>>(_base_name + "slip_resistance")),
    _slip_resistance_old(getMaterialPropertyOld<std::vector<Real>>(_base_name + "slip_resistance")),
    _slip_rate(declareProperty<std::vector<Real>>(_base_name + "slip_rate")),
    _slip_increment(declareProperty<std::vector<Real>>(_base_name + "slip_increment")),

    _slip_direction(_number_slip_systems),
    _slip_plane_normal(_number_slip_systems),
    _flow_direction(declareProperty<std::vector<RankTwoTensor>>(_base_name + "flow_direction")),
	_slip_direction_edge(declareProperty<std::vector<RealVectorValue>>("slip_direction_edge")),
	_slip_direction_screw(declareProperty<std::vector<RealVectorValue>>("slip_direction_screw")),
	_slip_plane_normalboth(declareProperty<std::vector<RealVectorValue>>("slip_plane_normalboth")),

    _tau(declareProperty<std::vector<Real>>(_base_name + "applied_shear_stress")),
    _crysrot(getMaterialProperty<RankTwoTensor>(_base_name + "crysrot")),
    _print_convergence_message(getParam<bool>("print_state_variable_convergence_error_messages")),
	
	_dislocation_mobile(declareProperty<std::vector<Real>>("dislocation_mobile")),
	_dislocation_immobile(declareProperty<std::vector<Real>>("dislocation_immobile")),
	_dislocation_immobile_old(getMaterialPropertyOld<std::vector<Real>>("dislocation_immobile")),
  _tau_b(declareProperty<std::vector<Real>>("back_stress")),
  _kappa(declareProperty<std::vector<Real>>("GND_density")),

    // resize vectors used in the consititutive slip hardening
  _hb(_number_slip_systems, 0.00),
  _slip_resistance_increment(_number_slip_systems, 0.00),
  _dislocation_mobile_increment(_number_slip_systems, 0.00),
  _dislocation_immobile_increment(_number_slip_systems, 0.00),
  _dv_dtau(_number_slip_systems, 0.00),
  _L_bar(_number_slip_systems, 0.00),

  // resize local variables realted to dislocation velocity
  t_wait(_number_slip_systems,0.00),
  t_run(_number_slip_systems,0.00),
  vel_run(_number_slip_systems,0.00),
  dislocation_density(_number_slip_systems,0.00),
  tau_b(_number_slip_systems,0.00),
  xi0(_number_slip_systems,0.00),
  tau_eff(_number_slip_systems,0.00),
  tau_effAbs(_number_slip_systems,0.00),
  tau_effSign(_number_slip_systems,0.00),
  slip_r(_number_slip_systems,0.00)
{
  getSlipSystems();
  sortCrossSlipFamilies();

  if (parameters.isParamSetByUser("number_cross_slip_directions"))
    _calculate_cross_slip = true;
  else
    _calculate_cross_slip = false;
}

void
StressUpdateCPBase::setMaterialVectorSize()
{

  _slip_direction_edge[_qp].resize(_slip_direction.size());
  _slip_direction_screw[_qp].resize(_slip_direction.size());
  _slip_plane_normalboth[_qp].resize(_slip_direction.size());

  _tau[_qp].resize(_number_slip_systems);

  _flow_direction[_qp].resize(_number_slip_systems);

  _slip_resistance[_qp].resize(_number_slip_systems);
  _slip_rate[_qp].resize(_number_slip_systems);
  
  _dislocation_mobile[_qp].resize(_slip_direction_edge[_qp].size());
  _dislocation_immobile[_qp].resize(_slip_direction_edge[_qp].size());
  _tau_b[_qp].resize(_number_slip_systems);


}

void
StressUpdateCPBase::initQpStatefulProperties()
{
  setMaterialVectorSize();
  for (const auto i : make_range(_number_slip_systems))
  {
    _flow_direction[_qp][i].zero();
    _tau[_qp][i] = 0.0;
  }

}

void
StressUpdateCPBase::getSlipSystems()
{
  bool orthonormal_error = false;

  // read in the slip system data from auxiliary text file
  MooseUtils::DelimitedFileReader _reader(_slip_sys_file_name);
  _reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::ROWS);
  _reader.read();

  for (const auto i : make_range(_number_slip_systems))
  {
    // initialize to zero
    _slip_direction[i].zero();
    _slip_plane_normal[i].zero();
  }

  if (_crystal_lattice_type == CrystalLatticeType::HCP)
    mooseError(" HCP slip-system is not compatible with this material class");
  else if (_crystal_lattice_type == CrystalLatticeType::BCC ||
           _crystal_lattice_type == CrystalLatticeType::FCC)
  {
    for (const auto i : make_range(_number_slip_systems))
    {
      // directly grab the raw data and scale it by the unit cell dimension
      for (const auto j : index_range(_reader.getData(i)))
      {
        if (j < LIBMESH_DIM)
          _slip_plane_normal[i](j) = _reader.getData(i)[j] / _unit_cell_dimension[j];
        else
          _slip_direction[i](j - LIBMESH_DIM) =
              _reader.getData(i)[j] * _unit_cell_dimension[j - LIBMESH_DIM];
      }
    }
  }

  for (const auto i : make_range(_number_slip_systems))
  {
    // normalize
    _slip_plane_normal[i] /= _slip_plane_normal[i].norm();
    _slip_direction[i] /= _slip_direction[i].norm();

    if (_crystal_lattice_type != CrystalLatticeType::HCP)
    {
      const auto magnitude = _slip_plane_normal[i] * _slip_direction[i];
      if (std::abs(magnitude) > libMesh::TOLERANCE)
      {
        orthonormal_error = true;
        break;
      }
    }
  }

  if (orthonormal_error)
    mooseError("StressUpdateCPBase Error: The slip system file contains a slip "
               "direction and plane normal pair that are not orthonormal in the Cartesian "
               "coordinate system.");
	
}

void
StressUpdateCPBase::storeDislocationMobilityInformation()
{
  for (const auto i : make_range(_number_slip_systems))
  {
	_slip_direction_edge[_qp][i] = _slip_direction[i];
	_slip_direction_edge[_qp][i] = _crysrot[_qp] * _slip_direction_edge[_qp][i]; 

	_slip_plane_normalboth[_qp][i] = _slip_plane_normal[i];
	_slip_plane_normalboth[_qp][i] = _crysrot[_qp] * _slip_plane_normalboth[_qp][i];
  }
}

void
StressUpdateCPBase::sortCrossSlipFamilies()
{
  if (_number_cross_slip_directions == 0)
  {
    _cross_slip_familes.resize(0);
    return;
  }

  // If cross slip does occur, then set up the system of vectors for the families
  _cross_slip_familes.resize(_number_cross_slip_directions);
  // and set the first index of each inner vector
  for (unsigned int i = 0; i < _number_cross_slip_directions; ++i)
    _cross_slip_familes[i].resize(1);

  // Sort the index of the slip system based vectors into separte families
  unsigned int family_counter = 1;
  _cross_slip_familes[0][0] = 0;

  for (unsigned int i = 1; i < _number_slip_systems; ++i)
  {
    for (unsigned int j = 0; j < family_counter; ++j)
    {
      // check to see if the slip system direction i matches any of the existing slip directions
      // First calculate the dot product
      Real dot_product = 0.0;
      for (const auto k : make_range(Moose::dim))
      {
        unsigned int check_family_index = _cross_slip_familes[j][0];
        dot_product += std::abs(_slip_direction[check_family_index](k) - _slip_direction[i](k));
      }
      // Then check if the dot product is one, if yes, add to family and break
      if (MooseUtils::absoluteFuzzyEqual(dot_product, 0.0))
      {
        _cross_slip_familes[j].push_back(i);
        if (_cross_slip_familes[j].size() > _number_cross_slip_planes)
          mooseError(
              "Exceeded the number of cross slip planes allowed in a single cross slip family");

        break; // exit the loop over the exisiting cross slip families and move to the next slip
               // direction
      }
      // The slip system in question does not belong to an existing family
      else if (j == (family_counter - 1) && !MooseUtils::absoluteFuzzyEqual(dot_product, 0.0))
      {
        if (family_counter > _number_cross_slip_directions)
          mooseError("Exceeds the number of cross slip directions specified for this material");

        _cross_slip_familes[family_counter][0] = i;
        family_counter++;
        break;
      }
    }
  }

  if (_print_convergence_message)
  {
    mooseWarning("Checking the slip system ordering now:");
    for (unsigned int i = 0; i < _number_cross_slip_directions; ++i)
    {
      Moose::out << "In cross slip family " << i << std::endl;
      for (unsigned int j = 0; j < _number_cross_slip_planes; ++j)
        Moose::out << " is the slip direction number " << _cross_slip_familes[i][j] << std::endl;
    }
  }
}

unsigned int
StressUpdateCPBase::identifyCrossSlipFamily(const unsigned int index)
{
  for (unsigned int i = 0; i < _number_cross_slip_directions; ++i)
    for (unsigned int j = 0; j < _number_cross_slip_planes; ++j)
      if (_cross_slip_familes[i][j] == index)
        return i;

  // Should never reach this statement
  mooseError("The supplied slip system index is not among the slip system families sorted.");
}

void
StressUpdateCPBase::calculateTotalPlasticDeformationGradientDerivative(
    RankFourTensor & dfpinvdpk2,
    const RankTwoTensor & inverse_plastic_deformation_grad_old,
    const RankTwoTensor & inverse_eigenstrain_deformation_grad_old)
{
  std::vector<Real> dslip_dtau(12, 0.0);
  std::vector<RankTwoTensor> dtaudpk2(12);
  std::vector<RankTwoTensor> dfpinvdslip(12);

  calculateConstitutiveSlipDerivative(dslip_dtau);

  for (const auto j : make_range(12))
  {
    dtaudpk2[j] = _flow_direction[_qp][j];
    dfpinvdslip[j] = -inverse_plastic_deformation_grad_old * _flow_direction[_qp][j];
    dfpinvdpk2 += (dfpinvdslip[j] * dslip_dtau[j] * _substep_dt).outerProduct(dtaudpk2[j]);
  }
}

void
StressUpdateCPBase::calculateSchmidTensor()
{
  calculateSchmidTensor(
      _number_slip_systems, _slip_plane_normal, _slip_direction, _flow_direction[_qp]);
}

void
StressUpdateCPBase::calculateSchmidTensor(
    const unsigned int & number_slip_systems,
    const std::vector<RealVectorValue> & plane_normal_vector,
    const std::vector<RealVectorValue> & direction_vector,
    std::vector<RankTwoTensor> & schmid_tensor)
{
  std::vector<RealVectorValue> local_direction_vector, local_plane_normal;
  local_direction_vector.resize(number_slip_systems);
  local_plane_normal.resize(number_slip_systems);

  // Update slip direction and normal with crystal orientation
  for (const auto i : make_range(_number_slip_systems))
  {
    local_direction_vector[i].zero();
    local_plane_normal[i].zero();

    for (const auto j : make_range(LIBMESH_DIM))
      for (const auto k : make_range(LIBMESH_DIM))
      {
        local_direction_vector[i](j) =
            local_direction_vector[i](j) + _crysrot[_qp](j, k) * direction_vector[i](k);

        local_plane_normal[i](j) =
            local_plane_normal[i](j) + _crysrot[_qp](j, k) * plane_normal_vector[i](k);
      }

    // Calculate Schmid tensor
    for (const auto j : make_range(LIBMESH_DIM))
      for (const auto k : make_range(LIBMESH_DIM))
      {
        schmid_tensor[i](j, k) = local_direction_vector[i](j) * local_plane_normal[i](k);
      }
  }
}

void
StressUpdateCPBase::calculateShearStress(
    const RankTwoTensor & pk2,
    const RankTwoTensor & inverse_eigenstrain_deformation_grad,
    const unsigned int & num_eigenstrains)
{
  if (!num_eigenstrains)
  {
    for (const auto i : make_range(_number_slip_systems))
      _tau[_qp][i] = pk2.doubleContraction(_flow_direction[_qp][i]);

    return;
  }

  RankTwoTensor eigenstrain_deformation_grad = inverse_eigenstrain_deformation_grad.inverse();
  for (const auto i : make_range(_number_slip_systems))
  {
    // compute PK2_hat using deformation gradient
    RankTwoTensor pk2_hat = eigenstrain_deformation_grad.det() *
                            eigenstrain_deformation_grad.transpose() * pk2 *
                            inverse_eigenstrain_deformation_grad.transpose();
    _tau[_qp][i] = pk2_hat.doubleContraction(_flow_direction[_qp][i]);
  }
}

void
StressUpdateCPBase::calculateEquivalentSlipIncrement(
    RankTwoTensor & equivalent_slip_increment)
{
  // Sum up the slip increments to find the equivalent plastic strain due to slip
  for (const auto i : make_range(_number_slip_systems))
    equivalent_slip_increment += _flow_direction[_qp][i] * _slip_rate[_qp][i] * _substep_dt;
}

void
StressUpdateCPBase::setQp(const unsigned int & qp)
{
  _qp = qp;
}

void
StressUpdateCPBase::setSubstepDt(const Real & substep_dt)
{
  _substep_dt = substep_dt;
}

bool
StressUpdateCPBase::isConstitutiveStateVariableConverged(
    const std::vector<Real> & current_var,
    const std::vector<Real> & var_before_update,
    const std::vector<Real> & previous_substep_var,
    const Real & tolerance)
{
  // sometimes the state variable size may not equal to the number of slip systems
  unsigned int sz = current_var.size();
  mooseAssert(current_var.size() == sz, "Current variable size does not match");
  mooseAssert(var_before_update.size() == sz, "Variable before update size does not match");
  mooseAssert(previous_substep_var.size() == sz, "Previous substep variable size does not match");

  bool is_converged = true;

  Real diff_val = 0.0;
  Real abs_prev_substep_val = 0.0;
  for (const auto i : make_range(sz))
  {
    diff_val = std::abs(var_before_update[i] - current_var[i]);
    abs_prev_substep_val = std::abs(previous_substep_var[i]);

    // set to false if the state variable is not converged
    if (abs_prev_substep_val < _zero_tol && diff_val > _zero_tol)
      is_converged = false;
    else if (abs_prev_substep_val > _zero_tol && diff_val > tolerance * abs_prev_substep_val)
      is_converged = false;
  }
  return is_converged;
}
