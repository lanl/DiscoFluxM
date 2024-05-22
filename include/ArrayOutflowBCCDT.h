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

#include "ArrayIntegratedBC.h"

class ArrayOutflowBCCDT : public ArrayIntegratedBC
{
public:
  static InputParameters validParams();

  ArrayOutflowBCCDT(const InputParameters & parameters);

protected:
  virtual void computeQpResidual(RealEigenVector & residual) override;
  virtual RealEigenVector computeQpJacobian() override;
  void computeSlipDirection();

private:
  //RealVectorValue _velocity;
  std::vector<RealVectorValue> _Dislo_Velocity;
  const MaterialProperty<std::vector<Real>> & _dislo_velocity_CP;
  
  // Dislocation character
  const enum class DislocationCharacter { edge, screw } _dislocationcharacter;  
  // Dislocation sign
  const enum class DislocationSign { positive, negative } _dislocationsign;
  
  const MaterialProperty<std::vector<RealVectorValue>> & _slip_direction_edge;
  const MaterialProperty<std::vector<RealVectorValue>> & _slip_plane_normalboth;
};
