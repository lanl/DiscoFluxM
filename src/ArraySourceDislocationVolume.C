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
#include "ArraySourceDislocationVolume.h"

registerMooseObject("TensorMechanicsApp", ArraySourceDislocationVolume);

InputParameters
ArraySourceDislocationVolume::validParams()
{
  InputParameters params = ArrayKernel::validParams();
  params.addClassDescription("Updates the dislocation density based on local evolution law of DiscoFlux crystal plasticity material model.");
  // parameters for local evolution of dislocation density.
  params.addParam<Real>("dislo_density_factor_CDT",1.0,"factor to convert the dislocation density from CDT to CP, this is for scaling of solution variable");
  params.addParam<Real>("C_multi", 8.96e-05, "parameter for dislocation multiplication");
  params.addParam<Real>("C_trap", 3.0e-03, "parameter for dislocation trapping");
  params.addParam<Real>("C_m_ann", 0.5, "parameter for dislocation mobile annihilation");
  params.addParam<Real>("C_im_ann", 0.5, "parameter for dislocation immobile annihilation");
  params.addParam<Real>("burgers_vector_mag", 3.52e-7, "magnitude of the Burger Vector in mm");
  
  MooseEnum dislocation_character("edge screw");
  params.addRequiredParam<MooseEnum>("dislocation_character", dislocation_character, "Character of dislocation");
  MooseEnum dislocation_sign("positive negative");
  params.addRequiredParam<MooseEnum>("dislocation_sign", dislocation_sign, "Sign of dislocation");

  params.addClassDescription(
      "The array Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
      "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  return params;
}

ArraySourceDislocationVolume::ArraySourceDislocationVolume(const InputParameters & parameters)
  : ArrayKernel(parameters),

	_dislo_density_factor_CDT(getParam<Real>("dislo_density_factor_CDT")),
 	_C_multi(getParam<Real>("C_multi")),
	_C_trap(getParam<Real>("C_trap")),
	_C_m_ann(getParam<Real>("C_m_ann")),
	_C_im_ann(getParam<Real>("C_im_ann")), 
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")), 
	
	_dislo_velocity_CP_edge(getMaterialProperty<std::vector<Real>>("dislo_velocity_edge")),

	_slip_direction_edge(getMaterialProperty<std::vector<RealVectorValue>>("slip_direction_edge")),
	_slip_plane_normalboth(getMaterialProperty<std::vector<RealVectorValue>>("slip_plane_normalboth")),

	_dislocationcharacter(getParam<MooseEnum>("dislocation_character").getEnum<DislocationCharacter>()),
	_dislocationsign(getParam<MooseEnum>("dislocation_sign").getEnum<DislocationSign>()), 

	_dislocation_mobile(getMaterialProperty<std::vector<Real>>("dislocation_mobile")),
  _dislocation_immobile(getMaterialProperty<std::vector<Real>>("dislocation_immobile")),
	_slip_rate(getMaterialProperty<std::vector<Real>>("slip_rate")),
	
	_dt(_fe_problem.dt())
	
{
}

void
ArraySourceDislocationVolume::computeQpResidual(RealEigenVector & residual)
{
	unsigned int _number_slip_systems = _count;
  Real A_f_ij, dislocation_forest;

  Real _dislocation_mobile_increment_mult = 0.00;
  Real _dislocation_mobile_increment_ann = 0.00;
  Real _dislocation_mobile_increment_trap = 0.00;

  RealEigenVector slip_rate, dislocation_mobile_increment;
  dislocation_mobile_increment.resize(_number_slip_systems);
  slip_rate.resize(_number_slip_systems);

  rateSourceTerm();
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {  
	dislocation_mobile_increment[i] = 0.00;

  dislocation_forest = 0.00;
	for (unsigned int j = 0; j < _number_slip_systems; ++j)
	{    
		A_f_ij = 0.5 * std::abs(_slip_plane_normalboth[_qp][i] * (_slip_direction_edge[_qp][j])) + 
				 0.5 * std::abs(_slip_plane_normalboth[_qp][i] * (_slip_plane_normalboth[_qp][i].cross(_slip_direction_edge[_qp][j])));
		dislocation_forest += A_f_ij *(_dislocation_mobile[_qp][j] + _dislocation_immobile[_qp][j]);
	}

  slip_rate[i] = _u[_qp][i] *_dislo_density_factor_CDT * _burgers_vector_mag * dislo_velocity[i];
	_dislocation_mobile_increment_mult = (_C_multi / std::pow(_burgers_vector_mag,2))* std::abs(slip_rate[i]);
	_dislocation_mobile_increment_trap = (_C_trap/_burgers_vector_mag)*std::pow(dislocation_forest,0.5)* std::abs(slip_rate[i]);
	_dislocation_mobile_increment_ann = (_C_m_ann*_dislocation_immobile[_qp][i])* std::abs(slip_rate[i]);

    dislocation_mobile_increment[i] =  (_dislocation_mobile_increment_mult - _dislocation_mobile_increment_trap - _dislocation_mobile_increment_ann);
	dislocation_mobile_increment[i] *= 1.0/_dislo_density_factor_CDT;  
	
	residual[i] = -_test[_i][_qp] * dislocation_mobile_increment[i];
  }  
}

RealEigenVector
ArraySourceDislocationVolume::computeQpJacobian()
{
  RealEigenVector jac ;
  jac.resize(_count);
  jac.setZero();
  
  rateSourceTerm();
  for (unsigned int i = 0; i < _var.count(); i++) 
  	jac[i] = -_test[_i][_qp] * ((0.00) + (_C_multi / std::pow(_burgers_vector_mag,2)) - (_C_trap/_burgers_vector_mag)*std::pow(_u[_qp][i]*_dislo_density_factor_CDT,0.5) - 
              (_C_m_ann*_u[_qp][i]*_dislo_density_factor_CDT) ) * (1/_dislo_density_factor_CDT) * (_phi[_j][_qp] * _burgers_vector_mag * std::abs(dislo_velocity[i]));

  return jac;

}

void
ArraySourceDislocationVolume::rateSourceTerm()
{
	unsigned int _number_slip_systems = _count; 

  dislo_velocity.resize(_number_slip_systems); // declared in the header file
  SourceRateTerm.resize(_number_slip_systems);
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
	  dislo_velocity[i] = 0.00;
	  dislo_velocity[i] = std::abs(_dislo_velocity_CP_edge[_qp][i]);
    SourceRateTerm[i] = _u[_qp][i] * dislo_velocity[i];
  }

}