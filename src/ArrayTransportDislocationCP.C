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
#include "ArrayTransportDislocationCP.h"

registerMooseObject("TensorMechanicsApp", ArrayTransportDislocationCP);

InputParameters
ArrayTransportDislocationCP::validParams()
{
  InputParameters params = ArrayKernel::validParams();
  params.addClassDescription("Continuum transport of dislocations(array variable) is modeled using advection model.");
  MooseEnum stabilization_method("none full", "full");
  params.addParam<MooseEnum>("stabilization_method",
                             stabilization_method,
                             "Stabilization method used for the transport term");
  MooseEnum dislocation_character("edge screw");
  params.addRequiredParam<MooseEnum>("dislocation_character", dislocation_character, "Character of dislocation");
  MooseEnum dislocation_sign("positive negative");
  params.addRequiredParam<MooseEnum>("dislocation_sign", dislocation_sign, "Sign of dislocation");
  params.addClassDescription(
      "The array Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
      "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  return params;
}

ArrayTransportDislocationCP::ArrayTransportDislocationCP(const InputParameters & parameters)
  : ArrayKernel(parameters),
	_dislo_velocity_CP_edge(getMaterialPropertyOld<std::vector<Real>>("dislo_velocity_edge")),
  _Dislo_Velocity_atQP (_count),

  _stabilization_method(getParam<MooseEnum>("stabilization_method").getEnum<UpwindingType>()),
	_u_nodal(_var.dofValuesOld()),
  _upwind_node(0),
  _dtotal_mass_out(0),
  work_vector(_count),
	_dislocationcharacter(getParam<MooseEnum>("dislocation_character").getEnum<DislocationCharacter>()),
	_dislocationsign(getParam<MooseEnum>("dislocation_sign").getEnum<DislocationSign>()), 
	_slip_direction_edge(getMaterialProperty<std::vector<RealVectorValue>>("slip_direction_edge")),
	_slip_plane_normalboth(getMaterialProperty<std::vector<RealVectorValue>>("slip_plane_normalboth"))
	
{
}


void
ArrayTransportDislocationCP::computeQpResidual(RealEigenVector & residual)
{
  /*
  Real _dislocation_mobile_increment_mult = 0.00;
  Real _dislocation_mobile_increment_ann = 0.00;
  Real _dislocation_mobile_increment_trap = 0.00;

  RealEigenVector slip_rate, dislocation_mobile_increment;
  dislocation_mobile_increment.resize(_number_slip_systems);
  slip_rate.resize(_number_slip_systems);
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {  
	dislocation_mobile_increment[i] = 0.00;

  slip_rate[i] = _u[_qp][i] *_dislo_density_factor_CDT * _burgers_vector_mag * _Dislo_Velocity_atQP[i];
	_dislocation_mobile_increment_mult = (_C_multi / std::pow(_burgers_vector_mag,2))* std::abs(slip_rate[i]);
	_dislocation_mobile_increment_trap = (_C_trap/_burgers_vector_mag)*std::pow(_C_m_ann*_dislocation_immobile[_qp][i],0.5)* std::abs(slip_rate[i]);
	_dislocation_mobile_increment_ann = (_C_m_ann*_dislocation_immobile[_qp][i])* std::abs(slip_rate[i]);

    dislocation_mobile_increment[i] =  (_dislocation_mobile_increment_mult - _dislocation_mobile_increment_trap - _dislocation_mobile_increment_ann);
	dislocation_mobile_increment[i] *= 1.0/_dislo_density_factor_CDT;  
	
	//residual[i] += -_test[_i][_qp] * dislocation_mobile_increment[i];
  }  
*/


  residual += (_u[_qp].asDiagonal() * negativeVelocityGradTestQp());
}

RealEigenVector
ArrayTransportDislocationCP::computeQpJacobian()
{	
  RealEigenVector jac ;
  jac.resize(_count);
  jac.setZero();

	jac = _phi[_j][_qp] * negativeVelocityGradTestQp();

	return jac;
}



void
ArrayTransportDislocationCP::computeJacobian()
{
  switch (_stabilization_method)
  {
    case UpwindingType::none:
      ArrayKernel::computeJacobian();
      break;
    case UpwindingType::full:
      ArrayKernel::computeJacobian();
      break;
  }
}

// Stabilization of the advection kernel 
void
ArrayTransportDislocationCP::stabilizedScheme(JacRes res_or_jac)
{
  const unsigned int num_nodes = _test.size();

  prepareVectorTag(_assembly, _var.number());

  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
    prepareMatrixTag(_assembly, _var.number(), _var.number());


  std::vector<Real> total_mass_out;
  std::vector<Real> total_in;
  std::vector<std::vector<bool>> _upwind_node;
  std::vector<std::vector<Real>> _dtotal_mass_out;

  _upwind_node.resize(_count,std::vector<bool>(num_nodes));
  _dtotal_mass_out.resize(_count,std::vector<Real>(num_nodes));
  RealEigenMatrix work_vector_res;
  work_vector_res.resize(_count,num_nodes); work_vector_res.setZero();
  work_vector.setZero(_count);
  RealEigenVector SpeedQp(_count); SpeedQp.setZero();
  RealEigenVector SpeedQpJac(_count); SpeedQpJac.setZero();

for (_i = 0; _i < num_nodes; _i++)
  {

   for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    initQpResidual();
      SpeedQp =negativeVelocityGradTestQp(); 
      mooseAssert(SpeedQp.size() == _count,
                  "Size of local residual is not equal to the number of array variable compoments");
      
      for (unsigned int i_var=0; i_var<_count; i_var++)
      {
      work_vector_res(i_var,_i) += _JxW[_qp] * _coord[_qp] * SpeedQp[i_var];
      }
  }

  for (unsigned int i_var=0; i_var<_count; i_var++)
  _upwind_node[i_var][_i] = (work_vector_res(i_var,_i) >= 0.0);
  }


  total_mass_out.assign(_count, 0.0);
  total_in.assign(_count, 0.0);

  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
    for (unsigned int i_var=0; i_var<_count; i_var++)
    _dtotal_mass_out[i_var].assign(num_nodes, 0.0);

for (unsigned int i_var = 0; i_var < _count; i_var++)
{
  for (unsigned int n = 0; n < num_nodes; n++)
  {
    if (_upwind_node[i_var][n])
    {
      work_vector_res(i_var,n) *= _u_nodal[n][i_var];
      total_mass_out[i_var] += work_vector_res(i_var,n); 
    }
    else                   
      total_in[i_var] -= work_vector_res(i_var,n); 
  }
  total_in[i_var] /= num_nodes;
}

  if (res_or_jac == JacRes::CALCULATE_RESIDUAL)
  {
    for (_i = 0; _i < _test.size(); _i++)
    {
      work_vector.setZero(_count);
      work_vector02.setZero(_count);
      prepareVectorTag(_assembly, _var.number());
      for (unsigned int i_var = 0; i_var < _count; i_var++)
       _local_re(_i) += _JxW[_qp] * _coord[_qp] * work_vector_res(i_var,_i);

    }
  accumulateTaggedLocalResidual();

  ArrayKernel::computeResidual();
  }

  else
  return;
  
  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
  {
  //ArrayKernel::computeOffDiagJacobian(jvar_num); return;
  //work_matrix = computeQpOffDiagJacobian(jvar) * _JxW[_qp] * _coord[_qp];
  }

}

RealEigenVector
ArrayTransportDislocationCP::negativeVelocityGradTestQp()
{    
  RealEigenVector negativeVelocityGradTest;
  negativeVelocityGradTest.resize(_count);
  negativeVelocityGradTest.setZero();
  computeDisloVelocity();
  rateSourceTerm();

  for (unsigned int i = 0; i < _var.count(); ++i) 
	{
    negativeVelocityGradTest[i] = (-1) * _Dislo_Velocity_atQP[i] * _grad_test[_i][_qp];
  }
	
  return negativeVelocityGradTest;
}

void
ArrayTransportDislocationCP::computeDisloVelocity()
{
	RealVectorValue _slip_direction_rotated;
	unsigned int _number_slip_systems = _count; 
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {

    _slip_direction_rotated = _slip_direction_edge[_qp][i]; // already rotated

		switch (_dislocationsign)
		{
			case DislocationSign::negative:
			_slip_direction_rotated *= (-1);
			break;
		}

    _Dislo_Velocity_atQP[i] = _dislo_velocity_CP_edge[_qp][i]*_slip_direction_rotated;

  }
}

void
ArrayTransportDislocationCP::rateSourceTerm()
{
	unsigned int _number_slip_systems = _count; 

  dislo_velocity.resize(_number_slip_systems);
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
	  dislo_velocity[i] = 0.00;
	  dislo_velocity[i] = std::abs(_dislo_velocity_CP_edge[_qp][i]);

  }
}


RealEigenMatrix
ArrayTransportDislocationCP::computeQpOffDiagJacobian(const MooseVariableFEBase & jvar)
{
    return ArrayKernel::computeQpOffDiagJacobian(jvar);
}


void
ArrayTransportDislocationCP::computeResidual()
{
  switch (_stabilization_method)
  {
    case UpwindingType::none:
      ArrayKernel::computeResidual();
      break;
    case UpwindingType::full:
      stabilizedScheme(JacRes::CALCULATE_RESIDUAL);
      break;
  }
}

void
ArrayTransportDislocationCP::computeOffDiagJacobian(const unsigned int jvar_num)
{
  switch (_stabilization_method)
  {
    case UpwindingType::none:
      ArrayKernel::computeOffDiagJacobian(jvar_num);
      break;
    case UpwindingType::full:
      jvar_num02 = jvar_num;
      stabilizedScheme(JacRes::CALCULATE_JACOBIAN);
      break;
  }
}