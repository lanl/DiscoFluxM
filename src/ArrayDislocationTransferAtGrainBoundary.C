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

#include "ArrayDislocationTransferAtGrainBoundary.h"


#include "StressUpdateCPBase.h"
#include "CrystalPlasticityDiscoflux.h"

registerMooseObject("TensorMechanicsApp", ArrayDislocationTransferAtGrainBoundary);

InputParameters
ArrayDislocationTransferAtGrainBoundary::validParams()
{
  InputParameters params = ArrayInterfaceKernel::validParams();
  params.addClassDescription("Transfers dislocation density as flux across the grain boundary.");

  params.addParam<Real>("density_critical", 1.0,"Critical density beyond which there will be dislocation transfer across Grain Boundary"); 
  params.addParam<Real>("tau_critical", 0.0,"Critical resolved shear stres beyond which there will be dislocation transfer across Grain Boundary");   
  params.addParam<Real>("scale_factor", 1.0,"scaling factor"); 
  MooseEnum dislocation_character("edge screw", "edge");
  params.addParam<MooseEnum>("dislocation_character", dislocation_character, "Character of dislocation");
  MooseEnum dislocation_sign("positive negative", "positive");
  params.addParam<MooseEnum>("dislocation_sign", dislocation_sign, "Sign of dislocation");
  params.addClassDescription("The kernel is utilized to establish flux equivalence on an interface for variables.");
  return params;
}

ArrayDislocationTransferAtGrainBoundary::ArrayDislocationTransferAtGrainBoundary(const InputParameters & parameters)
  : ArrayInterfaceKernel(parameters),
  
	_density_critical(getParam<Real>("density_critical")),
	_tau_critical(getParam<Real>("tau_critical")),
	_scale_factor(getParam<Real>("scale_factor")),
	_DD(_number_slip_systems),
    _dislo_velocity_CP(getMaterialProperty<std::vector<Real>>("dislo_velocity_edge")), // Velocity value (signed)
	_dislocationcharacter(getParam<MooseEnum>("dislocation_character").getEnum<DislocationCharacter>()),
	_dislocationsign(getParam<MooseEnum>("dislocation_sign").getEnum<DislocationSign>()), 
	_slip_direction_edge(getMaterialProperty<std::vector<RealVectorValue>>("slip_direction_edge")),
	_slip_plane_normalboth(getMaterialProperty<std::vector<RealVectorValue>>("slip_plane_normalboth")),
	_slip_direction_edge_neighbor(getNeighborMaterialProperty<std::vector<RealVectorValue>>("slip_direction_edge")),
	_slip_plane_normalboth_neighbor(getNeighborMaterialProperty<std::vector<RealVectorValue>>("slip_plane_normalboth")),
	
	_tau(getMaterialPropertyByName<std::vector<Real>>("applied_shear_stress")),
	_slip_resistance(getMaterialProperty<std::vector<Real>>("slip_resistance"))
	
{
}


void
ArrayDislocationTransferAtGrainBoundary::computeInterfaceAdvCoeff()
{
	 Real scalarProduct =0.00, density_initial, _crt;
	 Real S_GB, N_GB, L_GB;
	std::vector<std::vector<Real>> M_mod_GB, M_mod_GB_Norm, N_mod_GB;
	std::vector<Real> MaxValue_i, MaxValue_j;
	std::vector<int> Max_i, Max_j;
	RealVectorValue l1, l2, _slip_direction_rotated, _slip_direction_edge_rotated,
			 _slip_direction_rotated_neighbor, _slip_plane_normal_rotated, _slip_plane_normal_rotated_neighbor;

	//S_GB.resize(_number_slip_systems, std::vector<Real>(_number_slip_systems, 0.00));
	//N_GB.resize(_number_slip_systems, std::vector<Real>(_number_slip_systems, 0.00));
	//L_GB.resize(_number_slip_systems, std::vector<Real>(_number_slip_systems, 0.00));
	M_mod_GB.resize(_number_slip_systems, std::vector<Real>(_number_slip_systems, 0.00));
	M_mod_GB_Norm.resize(_number_slip_systems, std::vector<Real>(_number_slip_systems, 0.00));
	GB_Transfer_Coeff.resize(_number_slip_systems, std::vector<Real>(_number_slip_systems, 0.00));
	_index_j.resize(_number_slip_systems,0); 
	_DD.resize(_number_slip_systems,0.00);
	MaxValue_i.resize(_number_slip_systems,0.00);
	MaxValue_j.resize(_number_slip_systems,0.00);
	Max_i.resize(_number_slip_systems,0);
	Max_j.resize(_number_slip_systems,0);
	for (unsigned int i = 0; i < _number_slip_systems; i++)
	{
	_slip_plane_normal_rotated = _slip_plane_normalboth[_qp][i];
	_index_j[i] = 0;
	_DD[i] =0.00;
	_slip_direction_rotated = _slip_direction_edge[_qp][i];

	l1 = _slip_plane_normal_rotated.cross(_normals[_qp]);
	l1 /= l1.norm();
	for (unsigned int j = 0; j < _number_slip_systems; j++)
	{
		_slip_direction_rotated_neighbor = _slip_direction_edge_neighbor[_qp][j]; // Already rotated
		_slip_plane_normal_rotated_neighbor = _slip_plane_normalboth_neighbor[_qp][j];	
	
		l2 = _slip_plane_normal_rotated_neighbor.cross(-_normals[_qp]); 
		l2 /= l2.norm(); 
		L_GB = std::abs(l1 * l2);
		N_GB =  (_slip_plane_normal_rotated * _slip_plane_normal_rotated_neighbor);
		S_GB = (_slip_direction_rotated * _slip_direction_rotated_neighbor);
		M_mod_GB[i][j] = (L_GB * N_GB * S_GB);
	}
	Real max_coeff = 0.00;
	unsigned int index_max_coeff = 0;
	for (unsigned int j = 0; j < _number_slip_systems; j++)
	{
		GB_Transfer_Coeff[i][j] = M_mod_GB[i][j]; 
		if(GB_Transfer_Coeff[i][j] > max_coeff ) // && N_GB>0.00 // check for transfer between opposite sign
		{	index_max_coeff = j;
			max_coeff = GB_Transfer_Coeff[i][j];
		}
	}
	_index_j[i] = index_max_coeff;

		_slip_direction_rotated = _slip_direction_edge[_qp][i];
		switch (_dislocationsign) 
		{	case DislocationSign::negative:
			_slip_direction_rotated *= (-1);
			break;
		}
		velocity = _dislo_velocity_CP[_qp][i] * _slip_direction_rotated;
	
	//_crt = 1.00 + (_density_critical - 1.00) * (GB_Transfer_Coeff[i][i] - 1.0)/(0.5 - 1.0); 
	if((_u[_qp][i] > _density_critical) && std::abs(_tau[_qp][i])>_tau_critical) 
	         _DD[i] = _scale_factor * velocity * _normals[_qp]; 
	}

}

RealEigenVector
ArrayDislocationTransferAtGrainBoundary::computeQpResidual(Moose::DGResidualType type)
{  
  RealEigenVector res ; // = 0;
  res.resize(_count);
  res.setZero();

	if(_current_elem->subdomain_id() == _neighbor_elem->subdomain_id())  return res;
	
	computeInterfaceAdvCoeff(); 

  switch (type)
  {
    case Moose::Element:
	  for (unsigned int i = 0; i < _count; i++) 
	    res[i] += _u[_qp][i] * _DD[i] * _test[_i][_qp];
		//res[_index_j[i]] += _u[_qp][i] * _DD[i] * _test[_i][_qp];
    break;

    case Moose::Neighbor:
      for (unsigned int i = 0; i < _count; i++) 
	    res[i] += -_u[_qp][i] * _DD[i] * _test_neighbor[_i][_qp];
	    //res[_index_j[i]] += -_u[_qp][i] * _DD[i] * _test_neighbor[_i][_qp];
      break;
  }

  return res;
}

RealEigenVector
ArrayDislocationTransferAtGrainBoundary::computeQpJacobian(Moose::DGJacobianType type)
{
   RealEigenVector jac ;
  jac.resize(_count);
  jac.setZero();
  if(_current_elem->subdomain_id() == _neighbor_elem->subdomain_id())  return jac;
	
	computeInterfaceAdvCoeff();
	
  switch (type)
  {
    case Moose::ElementElement:
	 //for (unsigned int i = 0; i < _count; i++) jac[i] += _phi[_j][_qp] * _DD[i] * _test[_i][_qp];
     break;

    case Moose::NeighborElement:
	  if(true)
	  //for (unsigned int i = 0; i < _count; i++) jac[i] += -_phi_neighbor[_j][_qp] * _DD[i] * _test_neighbor[_i][_qp];
     break;

  }

  return jac;
}

