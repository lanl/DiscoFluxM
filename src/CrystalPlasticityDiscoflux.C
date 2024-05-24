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
#include "CrystalPlasticityDiscoflux.h"
#include "RankTwoTensor.h"

#include "SystemBase.h"
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

#include "Assembly.h" 
#include "MooseMesh.h"

#include "libmesh/node.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/remote_elem.h"

registerMooseObject("TensorMechanicsApp", CrystalPlasticityDiscoflux);

InputParameters
CrystalPlasticityDiscoflux::validParams()
{
  InputParameters params = StressUpdateCPBase::validParams();
  
  params.addClassDescription("Calculates the plastic slip based on DiscoFlux crystal plasticity material model.");
  
  // Parameter description in the DiscoFlux paper(International Journal of Plasticity 76 (2016) 111e129)
  params.addParam<Real>("lattice_friction", 10, "initial lattice friction strength of the material");  
  params.addParam<Real>("burgers_vector_mag",1.0e-07,"Magnitude of the Burgers vector in mm");
  params.addParam<Real>("dislo_density_initial",1.0e+05,"Initial dislocation density");
  
  params.addParam<Real>("C_multi", 8.5e-06, "parameter for dislocation multiplication");
  params.addParam<Real>("C_trap", 5.5e-03, "parameter for dislocation trapping");
  params.addParam<Real>("C_m_ann", 0.5, "parameter for dislocation mobile annihilation");
  params.addParam<Real>("C_im_ann", 0.5, "parameter for dislocation immobile annihilation");
  params.addParam<Real>("Coeff_hardening", 0.75, "parameter to control the material hardening");

  params.addParam<Real>("q1", 0.1, "material parameter");
  params.addParam<Real>("q2", 1.9, "material parameter");
  params.addParam<Real>("c1", 2.0, "material parameter");
  params.addParam<Real>("temp", 300, "Temperature");
  params.addParam<Real>("dislo_density_factor_CDT",1.0e+06,"factor to convert the dislocation density from CDT to CP");

  return params;
}

CrystalPlasticityDiscoflux::CrystalPlasticityDiscoflux(
    const InputParameters & parameters)
  : StressUpdateCPBase(parameters),
	_lattice_friction(getParam<Real>("lattice_friction")),
	
	_burgers_vector_mag(getParam<Real>("burgers_vector_mag")),
	_dislo_density_initial(getParam<Real>("dislo_density_initial")),
	_dislo_density_factor_CDT(getParam<Real>("dislo_density_factor_CDT")),
	_C_multi(getParam<Real>("C_multi")),
	_C_trap(getParam<Real>("C_trap")),
	_C_m_ann(getParam<Real>("C_m_ann")),
	_C_im_ann(getParam<Real>("C_im_ann")),
	_Coeff_hardening(getParam<Real>("Coeff_hardening")),

	_q1(getParam<Real>("q1")),
	_q2(getParam<Real>("q2")),
	_c1(getParam<Real>("c1")),
	_temp(getParam<Real>("temp")),
	
	//_Rho_EdgePositive_01(coupledValue("Rho_EdgePositive_01")),
	_DD_EdgePositive(coupledArrayValue("DD_EdgePositive")), 
	_DD_EdgeNegative(coupledArrayValue("DD_EdgeNegative")), 
	_DD_ScrewPositive(coupledArrayValue("DD_ScrewPositive")), 
	_DD_ScrewNegative(coupledArrayValue("DD_ScrewNegative")),
	
	_DD_EdgePositive_Grad(coupledArrayGradient("DD_EdgePositive")),
	_DD_EdgeNegative_Grad(coupledArrayGradient("DD_EdgeNegative")),
	_DD_ScrewPositive_Grad(coupledArrayGradient("DD_ScrewPositive")),
	_DD_ScrewNegative_Grad(coupledArrayGradient("DD_ScrewNegative")),
	
	_dislocation_immobile(declareProperty<std::vector<Real>>("dislocation_immobile")),
	_dislocation_immobile_old(getMaterialPropertyOld<std::vector<Real>>("dislocation_immobile")),
	
	_dislo_velocity_edge(declareProperty<std::vector<Real>>("dislo_velocity_edge")),
	_dislo_velocity_screw(declareProperty<std::vector<Real>>("dislo_velocity_screw")),
	_tau_old(getMaterialPropertyOld<std::vector<Real>>("applied_shear_stress")),
	_GND_density(getMaterialProperty<std::vector<Real>>("GND_density")),

	// DDC related variables
	_kappa(declareProperty<std::vector<Real>>("kappa")),
	_DD_grad(_number_slip_systems, 0.00),
	_tau_b_local(_number_slip_systems, 0.00),

  // resize local caching vectors used for substepping
  _previous_substep_slip_resistance(_number_slip_systems, 0.00),
  _previous_substep_dislocation_mobile(_number_slip_systems, 0.00),
  _previous_substep_dislocation_immobile(_number_slip_systems, 0.00),
  
  _slip_resistance_before_update(_number_slip_systems, 0.00),
  _dislocation_mobile_before_update(_number_slip_systems, 0.00),
  _dislocation_immobile_before_update(_number_slip_systems, 0.00)

{
  
}

void
CrystalPlasticityDiscoflux::initQpStatefulProperties()
{
  StressUpdateCPBase::initQpStatefulProperties();
    
  _slip_direction_edge[_qp].resize(_number_slip_systems);
  _slip_direction_screw[_qp].resize(_number_slip_systems);
  _slip_plane_normalboth[_qp].resize(_number_slip_systems);


  _dislocation_immobile[_qp].resize(_number_slip_systems);
  _dislo_velocity_edge[_qp].resize(_number_slip_systems);
  _dislo_velocity_screw[_qp].resize(_number_slip_systems);
  _kappa[_qp].resize(_number_slip_systems);
  
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {  
    _slip_direction_edge[_qp][i].zero();
	_slip_direction_screw[_qp][i].zero();
	_slip_plane_normalboth[_qp][i].zero();

    _slip_resistance[_qp][i] = 5*_lattice_friction; // approximate initial value
    _slip_rate[_qp][i] = 0.0;
	_dislocation_mobile[_qp][i] = 4 * _dislo_density_initial;
	_dislocation_immobile[_qp][i] = 4.0 * _dislo_density_initial;
	_dislo_velocity_edge[_qp][i] = 0.00;
	_dislo_velocity_screw[_qp][i] = 0.00;
	_kappa[_qp][i] = 0.0;
	
  }
  
}

void
CrystalPlasticityDiscoflux::setInitialConstitutiveVariableValues()
{
  storeDislocationMobilityInformation();
  _slip_resistance[_qp] = _slip_resistance_old[_qp];
  _previous_substep_slip_resistance = _slip_resistance_old[_qp];
  
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
	_dislocation_mobile[_qp][i] = (_DD_EdgePositive[_qp][i] + _DD_EdgeNegative[_qp][i] + _DD_ScrewPositive[_qp][i] + _DD_ScrewNegative[_qp][i]) * _dislo_density_factor_CDT;
	//_dislocation_mobile[_qp][i] = 4.0 * _dislo_density_factor_CDT;
	_previous_substep_dislocation_mobile[i] = _dislocation_mobile[_qp][i] ;
  }
	_dislocation_immobile[_qp] = _dislocation_immobile_old[_qp];
    _previous_substep_dislocation_immobile = _dislocation_immobile_old[_qp];
}

void
CrystalPlasticityDiscoflux::setSubstepConstitutiveVariableValues()
{
  _slip_resistance[_qp] = _previous_substep_slip_resistance;
  _dislocation_immobile[_qp] = _previous_substep_dislocation_immobile;
}

bool
CrystalPlasticityDiscoflux::calculateSlipRate()
{	
// compute dislocation velocity according to DiscoFlux material model
getDisloVelocity();

for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
	//_dislocation_mobile[_qp][i] = ((_DD_EdgePositive[_qp][i] + _DD_EdgeNegative[_qp][i]) * _dislo_velocity_edge[_qp][i]
  	//						 + (_DD_ScrewPositive[_qp][i] + _DD_ScrewNegative[_qp][i]) * _dislo_velocity_screw[_qp][i]) * _dislo_density_factor_CDT * _burgers_vector_mag;
	_slip_rate[_qp][i] = _dislocation_mobile[_qp][i] * _burgers_vector_mag * _dislo_velocity_edge[_qp][i]; // For edge only
  }
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
	_tau_b[_qp][i] = _tau_b_local[i];
	_kappa[_qp][i] = (_DD_EdgePositive[_qp][i] -_DD_EdgeNegative[_qp][i])*_dislo_density_factor_CDT; // (_DD_ScrewPositive[_qp][i] - _DD_ScrewNegative[_qp][i]);
  }
  return true;
}

void
CrystalPlasticityDiscoflux::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
    if (MooseUtils::absoluteFuzzyEqual(_tau[_qp][i], 0.0))
      dslip_dtau[i] = 0.0;
    else
	{
	    dslip_dtau[i] = (_DD_EdgePositive[_qp][i] + _DD_EdgeNegative[_qp][i]) * _dislo_density_factor_CDT * _burgers_vector_mag * _dv_dtau[i];
	}	
 }
}

void
CrystalPlasticityDiscoflux::cacheStateVariablesBeforeUpdate()
{
  _slip_resistance_before_update = _slip_resistance[_qp];
  _dislocation_immobile_before_update = _dislocation_immobile[_qp];
}

void
CrystalPlasticityDiscoflux::calculateStateVariableEvolutionRateComponent()
{ 
  // calculate dislocation density increment
  getDDIncrements();
}


// Calculate Dislocation Density increment
void
CrystalPlasticityDiscoflux::getDDIncrements()
{  
  Real small2 = 1.0e-5;
  Real A_f_ij, dislocation_forest;

  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {  
	_dislocation_mobile_increment[i] = 0.00;
	_dislocation_immobile_increment[i] = 0.00;
	dislocation_forest = 0.00;
	for (unsigned int j = 0; j < _number_slip_systems; ++j)
	{    
		A_f_ij = 0.5 * std::abs(_slip_plane_normalboth[_qp][i] * (_slip_direction_edge[_qp][j])) + 
				 0.5 * std::abs(_slip_plane_normalboth[_qp][i] * (_slip_plane_normalboth[_qp][i].cross(_slip_direction_edge[_qp][j])));
		dislocation_forest += A_f_ij *(_dislocation_mobile[_qp][j] + _dislocation_immobile[_qp][j]);
	}

	if (std::abs(_slip_rate[_qp][i]) > small2) {
	_dislocation_immobile_increment[i] = (_C_trap/_burgers_vector_mag)*std::pow(dislocation_forest,0.5)* std::abs(_slip_rate[_qp][i])
		- (_C_im_ann*_dislocation_immobile[_qp][i])* std::abs(_slip_rate[_qp][i]);	
	 
	} 
	else
	 {
	  _dislocation_mobile_increment[i] = 0.0;
	  _dislocation_immobile_increment[i] = 0.0;
	  
	 }
 }
}

bool
CrystalPlasticityDiscoflux::updateStateVariables()
{
	Real Hij, eff_dislocation_density = 0.00;

     for (unsigned int i = 0; i < _number_slip_systems; ++i)
      {  
	    _dislocation_immobile[_qp][i] = _previous_substep_dislocation_immobile[i] + _dislocation_immobile_increment[i] * _substep_dt;
	  
      }

	  for (unsigned int i = 0; i < _number_slip_systems; ++i)
	  {
		  eff_dislocation_density = 0.00;
		  for (unsigned int j = 0; j < _number_slip_systems; ++j)
		  {
			if(i==j ? Hij=1.0: Hij=0.1) 
		    eff_dislocation_density += Hij*2.0*_dislocation_immobile[_qp][j]; // similar proportion of mobile and immobile dislocations
		  }
	  _slip_resistance[_qp][i] = _lattice_friction + _Coeff_hardening*mu*_burgers_vector_mag*std::sqrt(eff_dislocation_density);
	  }

  return true;
}

bool
CrystalPlasticityDiscoflux::areConstitutiveStateVariablesConverged()
{
	bool flagSlipResistanceConverged;

// Check convergence of _slip_resistance
   flagSlipResistanceConverged = isConstitutiveStateVariableConverged(_slip_resistance[_qp],
                                              _slip_resistance_before_update,
                                              _previous_substep_slip_resistance,
                                              _resistance_tol);

	return flagSlipResistanceConverged;
}

void
CrystalPlasticityDiscoflux::updateSubstepConstitutiveVariableValues()
{
  _previous_substep_slip_resistance = _slip_resistance[_qp];
  _previous_substep_dislocation_immobile = _dislocation_immobile[_qp];
}


//----Compute the DDC informed dislocation delocity

void
CrystalPlasticityDiscoflux::getDisloVelocity()
{ 
Stress_internal.zero();
Real dislocationsPositive, dislocationsNegative;
  
for (unsigned int i = 0; i < _number_slip_systems; ++i)
{	
	slip_direction_rotated = _slip_direction[i];
	slip_plane_normal_rotated = _slip_plane_normal[i];
	_L_bar[i] = std::pow((_dislocation_mobile[_qp][i] + _dislocation_immobile[_qp][i]),-0.5); 
	dislocationsPositive = _DD_EdgePositive[_qp][i];
	dislocationsNegative = _DD_EdgeNegative[_qp][i];
	_DD_grad[i](0) = (_DD_EdgePositive_Grad[_qp](i) - _DD_EdgeNegative_Grad[_qp](i))*_dislo_density_factor_CDT;
	_DD_grad[i](1) = (_DD_EdgePositive_Grad[_qp](i+_number_slip_systems) - _DD_EdgeNegative_Grad[_qp](i+_number_slip_systems))*_dislo_density_factor_CDT;
	_DD_grad[i](2) = (_DD_EdgePositive_Grad[_qp](i+2*_number_slip_systems) - _DD_EdgeNegative_Grad[_qp](i+2*_number_slip_systems))*_dislo_density_factor_CDT;
	_tau_b_local[i] = 0.2*(( mu * std::pow(_L_bar[i],1))/(2*3.141*(1-nu)))*_burgers_vector_mag * (_DD_grad[i]*slip_direction_rotated);
	Stress_internal += _tau_b_local[i]*(libMesh::outer_product(slip_direction_rotated, slip_plane_normal_rotated) + 
						libMesh::outer_product(slip_plane_normal_rotated, slip_direction_rotated));
}
for (unsigned int i = 0; i < _number_slip_systems; ++i) 
{
	slip_direction_rotated = _slip_direction[i];
	slip_plane_normal_rotated = _slip_plane_normal[i];
	_tau_b_local[i] = Stress_internal.contract(libMesh::outer_product(slip_direction_rotated, slip_plane_normal_rotated));
}


	//  compute velocity for each slip system
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
   {
	   tau_eff[i] = 0.00;
		t_wait[i]  = 0.00;
		t_run[i]   = 0.00;
		vel_run[i] = 0.00;
		_dislo_velocity_edge[_qp][i] =0.00;
		_dv_dtau[i] = 0.00;
   }
   
   for (unsigned int i = 0; i < _number_slip_systems; ++i)
   {
	   slip_r[i]  = _slip_resistance[_qp][i];

	   tau_eff[i] = (_tau[_qp][i] - _tau_b_local[i]);
	   tau_effAbs[i] = std::abs(tau_eff[i]);
	   tau_effSign[i] = std::copysign(1.0, tau_eff[i]);
   }

   deltaG0 = g0*mu*std::pow(_burgers_vector_mag,3)*1.0e-3;
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {	   
	  inner = 1.0 - std::pow((tau_effAbs[i] / slip_r[i] ),_q1);
	  if (inner > 0.00)
	  {
		  deltaG  = deltaG0*( std::pow(inner,_q2) );
		  exp_arg = deltaG / (boltz*_temp);
		  t_wait[i] = (exp(exp_arg) - 1.0) / omega0;
	  }
	  else
		 t_wait[i] = 0.00;
  }

  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  { 
	if (tau_effAbs[i]>small2)
	   {
	     xi0[i] = B0*vcrit / (2*_burgers_vector_mag*tau_effAbs[i]);
	     vel_run[i] = vcrit*(std::pow((std::pow(xi0[i],2)+1),0.5) - xi0[i]);
	    }
  }

  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
	if (vel_run[i] > small2)
	   {  
		t_run[i] = _L_bar[i] / vel_run[i];
		_dislo_velocity_edge[_qp][i] = tau_effSign[i]*_L_bar[i] / (t_wait[i] + t_run[i]);
		_dv_dtau[i] = 0.00;

		inner = 1.0 - std::pow((tau_effAbs[i] / slip_r[i] ),_q1);
		deltaG  = deltaG0*( std::pow(inner,_q2) );
		exp_arg = deltaG / (boltz*_temp);
		_dv_dtau[i] = 0.00; 
	
		}
	else
	  {
	  _dislo_velocity_edge[_qp][i] = 0.00;
	  _dv_dtau[i] = 0.00;
	  }
  }

// compute screw dislocation velocity.
//Stress_internal.zero();
for (unsigned int i = 0; i < _number_slip_systems; ++i)
{	
	slip_direction_rotated = _slip_direction[i];
	slip_plane_normal_rotated = _slip_plane_normal[i];
	_L_bar[i] = std::pow((_dislocation_mobile[_qp][i] + _dislocation_immobile[_qp][i]),-0.5); 
	dislocationsPositive = _DD_EdgePositive[_qp][i];
	dislocationsNegative = _DD_EdgeNegative[_qp][i];
	_DD_grad[i](0) = (_DD_ScrewPositive_Grad[_qp](i) - _DD_ScrewNegative_Grad[_qp](i))*_dislo_density_factor_CDT;
	_DD_grad[i](1) = (_DD_ScrewPositive_Grad[_qp](i+_number_slip_systems) - _DD_ScrewNegative_Grad[_qp](i+_number_slip_systems))*_dislo_density_factor_CDT;
	_DD_grad[i](2) = (_DD_ScrewPositive_Grad[_qp](i+2*_number_slip_systems) - _DD_ScrewNegative_Grad[_qp](i+2*_number_slip_systems))*_dislo_density_factor_CDT;
	_tau_b_local[i] = 0.2*(( mu * std::pow(_L_bar[i],1))/(2*3.141*(1-nu)))*_burgers_vector_mag * (_DD_grad[i]*slip_direction_rotated);
	Stress_internal += _tau_b_local[i]*(libMesh::outer_product(slip_direction_rotated, slip_plane_normal_rotated) + 
						libMesh::outer_product(slip_plane_normal_rotated, slip_direction_rotated));
}

for (unsigned int i = 0; i < _number_slip_systems; ++i) 
{
	slip_direction_rotated = _slip_direction[i];
	slip_plane_normal_rotated = _slip_plane_normal[i];
	_tau_b_local[i] = Stress_internal.contract(libMesh::outer_product(slip_direction_rotated, slip_plane_normal_rotated));
}


  for (unsigned int i = 0; i < _number_slip_systems; ++i)
   {
	   tau_eff[i] = 0.00;
		t_wait[i]  = 0.00;
		t_run[i]   = 0.00;
		vel_run[i] = 0.00;
		_dislo_velocity_screw[_qp][i] =0.00;
		_dv_dtau[i] = 0.00;
   }
   
   for (unsigned int i = 0; i < _number_slip_systems; ++i)
   {
	   slip_r[i]  = _slip_resistance[_qp][i];

	   tau_eff[i] = (_tau[_qp][i] - _tau_b_local[i]);
	   tau_effAbs[i] = std::abs(tau_eff[i]);
	   tau_effSign[i] = std::copysign(1.0, tau_eff[i]);
   }

   deltaG0 = g0*mu*std::pow(_burgers_vector_mag,3)*1.0e-3;
  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {	   
	  inner = 1.0 - std::pow((tau_effAbs[i] / slip_r[i] ),_q1);
	  if (inner > 0.00)
	  {
		  deltaG  = deltaG0*( std::pow(inner,_q2) );
		  exp_arg = deltaG / (boltz*_temp);
		  t_wait[i] = (exp(exp_arg) - 1.0) / omega0;
	  }
	  else
		 t_wait[i] = 0.00;
  }

  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  { 
	if (tau_effAbs[i]>small2)
	   {
	     xi0[i] = B0*vcrit / (2*_burgers_vector_mag*tau_effAbs[i]);
	     vel_run[i] = vcrit*(std::pow((std::pow(xi0[i],2)+1),0.5) - xi0[i]);
	    }
  }

  for (unsigned int i = 0; i < _number_slip_systems; ++i)
  {
	if (vel_run[i] > small2)
	   {  
		t_run[i] = _L_bar[i] / vel_run[i];
		_dislo_velocity_screw[_qp][i] = tau_effSign[i]*_L_bar[i] / (t_wait[i] + t_run[i]);
		_dv_dtau[i] = 0.00;

		inner = 1.0 - std::pow((tau_effAbs[i] / slip_r[i] ),_q1);
		deltaG  = deltaG0*( std::pow(inner,_q2) );
		exp_arg = deltaG / (boltz*_temp);
		_dv_dtau[i] = 0.00; 
	
		}
	else
	  {
	  _dislo_velocity_screw[_qp][i] = 0.00;
	  _dv_dtau[i] = 0.00;
	  }
  }

}