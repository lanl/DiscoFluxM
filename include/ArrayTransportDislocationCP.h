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
#pragma once

#include "ArrayKernel.h"

class ArrayTransportDislocationCP : public ArrayKernel
{
public:
  static InputParameters validParams();

  ArrayTransportDislocationCP(const InputParameters & parameters);

protected:
  virtual void computeQpResidual(RealEigenVector & residual) override;
  virtual RealEigenVector computeQpJacobian() override;
  virtual RealEigenMatrix computeQpOffDiagJacobian(const MooseVariableFEBase & jvar) override;

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar_num) override;
  void computeDisloVelocity();
  virtual void rateSourceTerm() override;
  RealEigenVector negativeVelocityGradTestQp();  

  const MaterialProperty<std::vector<Real>> & _dislo_velocity_CP_edge;
  std::vector<RealVectorValue> _Dislo_Velocity_atQP ;

  enum class JacRes
  {
    CALCULATE_RESIDUAL = 0,
    CALCULATE_JACOBIAN = 1
  };
  const enum class UpwindingType { none, full } _stabilization_method;
  const ArrayVariableValue & _u_nodal;
  std::vector<std::vector<bool>> _upwind_node;
  std::vector<std::vector<Real>> _dtotal_mass_out;
  unsigned int jvar_num02;
  
  void stabilizedScheme(JacRes res_or_jac);
  RealEigenVector work_vector, work_vector02;
  RealEigenMatrix work_matrix;
  
  // Dislocation character
  const enum class DislocationCharacter { edge, screw } _dislocationcharacter;  
  // Dislocation sign
  const enum class DislocationSign { positive, negative } _dislocationsign;
  // data related to slip_systems
  const MaterialProperty<std::vector<RealVectorValue>> & _slip_direction_edge;
  const MaterialProperty<std::vector<RealVectorValue>> & _slip_plane_normalboth;


  RealEigenVector dislo_velocity;
};
