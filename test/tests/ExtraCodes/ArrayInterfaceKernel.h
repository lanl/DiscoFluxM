//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// local includes
#include "InterfaceKernelBase.h"

#define TemplateVariableValue typename OutputTools<T>::VariableValue
#define TemplateVariableGradient typename OutputTools<T>::VariableGradient
#define TemplateVariablePhiValue typename OutputTools<T>::VariablePhiValue
#define TemplateVariablePhiGradient typename OutputTools<T>::VariablePhiGradient
#define TemplateVariableTestValue typename OutputTools<T>::VariableTestValue
#define TemplateVariableTestGradient typename OutputTools<T>::VariableTestGradient

// Forward Declarations
template <typename T>
class InterfaceKernelTempl02;

typedef InterfaceKernelTempl02<RealEigenVector> ArrayInterfaceKernel;
template <typename T>
class InterfaceKernelTempl02 : public InterfaceKernelBase, public NeighborMooseVariableInterface<T>
{
public:
  static InputParameters validParams();

  InterfaceKernelTempl02(const InputParameters & parameters);

  /// The primary variable that this interface kernel operates on
  virtual const MooseVariableFE<T> & variable() const override { return _var; }
  virtual const MooseVariableFE<T> & neighborVariable() const override { return _neighbor_var; }

  // For computing residuals and jacobians
  virtual void computeElemNeighResidual(Moose::DGResidualType type);
  virtual void computeElemNeighJacobian(Moose::DGJacobianType type);
  virtual void computeOffDiagElemNeighJacobian(Moose::DGJacobianType type, unsigned int jvar) ;
  virtual void computeElementOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;
  virtual void computeNeighborOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) override;
  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeResidualAndJacobian() override;

  /// Compute residuals at quadrature points
  virtual RealEigenVector computeQpResidual(Moose::DGResidualType type) = 0;
  virtual RealEigenVector computeQpJacobian(Moose::DGJacobianType type) = 0;


  // initialized required array for residual and jacobian
  virtual void initQpResidual(Moose::DGResidualType /* type */) {}
  virtual void initQpJacobian(Moose::DGJacobianType /* type */) {}
  virtual void initQpOffDiagJacobian(Moose::DGJacobianType /* type */, unsigned int /* jvar */) {}
  
protected:
  MooseVariableFE<T> & _var;
  const MooseArray<Point> & _normals;
  const TemplateVariableValue & _u;
  const TemplateVariableGradient & _grad_u;
  const TemplateVariablePhiValue & _phi;
  const TemplateVariablePhiGradient & _grad_phi;
  const TemplateVariableTestValue & _test;
  const TemplateVariableTestGradient & _grad_test;
  const MooseVariableFE<T> & _neighbor_var;
  const TemplateVariableValue & _neighbor_value;
  const TemplateVariableGradient & _grad_neighbor_value;
  const TemplateVariablePhiValue & _phi_neighbor;
  const TemplateVariablePhiGradient & _grad_phi_neighbor;
  const TemplateVariableTestValue & _test_neighbor;
  const TemplateVariableTestGradient & _grad_test_neighbor;
  DenseMatrix<Number> _local_kxx;
  const unsigned int _count;
  
  private:
  RealEigenVector _work_vector;

};
