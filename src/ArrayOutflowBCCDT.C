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
#include "ArrayOutflowBCCDT.h"

registerMooseObject("TensorMechanicsApp", ArrayOutflowBCCDT);

InputParameters
ArrayOutflowBCCDT::validParams()
{
  InputParameters params = ArrayIntegratedBC::validParams();
  params.addClassDescription("Imposes free flow of flux of a quantity(array variable) across the boundary.");
  MooseEnum dislocation_character("edge screw");
  params.addRequiredParam<MooseEnum>("dislocation_character", dislocation_character, "Character of dislocation");
  MooseEnum dislocation_sign("positive negative");
  params.addRequiredParam<MooseEnum>("dislocation_sign", dislocation_sign, "Sign of dislocation");
  params.addClassDescription(
      "Enforces a outflow boundary condition for advection kernels");
  return params;
}

ArrayOutflowBCCDT::ArrayOutflowBCCDT(const InputParameters & parameters)
  : ArrayIntegratedBC(parameters),
    _Dislo_Velocity(_count),
    _dislo_velocity_CP(getMaterialProperty<std::vector<Real>>("dislo_velocity_edge")), // Velocity value (signed)
	_dislocationcharacter(getParam<MooseEnum>("dislocation_character").getEnum<DislocationCharacter>()),
	_dislocationsign(getParam<MooseEnum>("dislocation_sign").getEnum<DislocationSign>()), 
	_slip_direction_edge(getMaterialProperty<std::vector<RealVectorValue>>("slip_direction_edge")),
	_slip_plane_normalboth(getMaterialProperty<std::vector<RealVectorValue>>("slip_plane_normalboth"))
{
}

void
ArrayOutflowBCCDT::computeQpResidual(RealEigenVector & residual)
{
	computeSlipDirection();
  for (unsigned int i = 0; i < _count; i++)
  {
	  residual[i] = _test[_i][_qp] * (_u[_qp][i] - 1.0) * (_Dislo_Velocity[i] * _normals[_qp]);
  }
}

RealEigenVector
ArrayOutflowBCCDT::computeQpJacobian()
{
	computeSlipDirection();
  RealEigenVector jac ;
  jac.resize(_count);
  jac.setZero();
  for (unsigned int i = 0; i < _count; i++)
  {
	  jac[i] = _test[_i][_qp] * _phi[_j][_qp] * (_Dislo_Velocity[i] * _normals[_qp]);
  }
  return jac;
}


void
ArrayOutflowBCCDT::computeSlipDirection()
{
  RealVectorValue _slip_direction_rotated;
	
	for (unsigned int i = 0; i < _count; i++)
	{	
		_slip_direction_rotated = _slip_direction_edge[_qp][i];
		_Dislo_Velocity[i] = _dislo_velocity_CP[_qp][i] * _slip_direction_rotated;
	}
}

