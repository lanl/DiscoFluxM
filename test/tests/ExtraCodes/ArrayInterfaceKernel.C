//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html


#include "ArrayInterfaceKernel.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

#include "libmesh/quadrature.h"

template <typename T>
InputParameters
InterfaceKernelTempl02<T>::validParams()
{
  InputParameters params = InterfaceKernelBase::validParams();
  if (std::is_same<T, RealEigenVector>::value)
    params.registerBase("ArrayInterfaceKernel");
  else
    ::mooseError("unsupported InterfaceKernelTempl02 specialization");
  return params;
}

template <typename T>
InterfaceKernelTempl02<T>::InterfaceKernelTempl02(const InputParameters & parameters)
  : InterfaceKernelBase(parameters),
    NeighborMooseVariableInterface<T>(this,
                                      false,
                                      Moose::VarKindType::VAR_NONLINEAR,
                                      std::is_same<T, Real>::value
                                          ? Moose::VarFieldType::VAR_FIELD_STANDARD
                                          : Moose::VarFieldType::VAR_FIELD_ARRAY),
    _var(*this->mooseVariable()),
    _normals(_assembly.normals()),
    _u(_is_implicit ? _var.sln() : _var.slnOld()),
    _grad_u(_is_implicit ? _var.gradSln() : _var.gradSlnOld()),
    _phi(_assembly.phiFace(_var)),
    _grad_phi(_assembly.gradPhiFace(_var)),
    _test(_var.phiFace()),
    _grad_test(_var.gradPhiFace()),
    _neighbor_var(*getVarHelper<MooseVariableFE<T>>("neighbor_var", 0)),
    _neighbor_value(_is_implicit ? _neighbor_var.slnNeighbor() : _neighbor_var.slnOldNeighbor()),
    _grad_neighbor_value(_neighbor_var.gradSlnNeighbor()),
    _phi_neighbor(_assembly.phiFaceNeighbor(_neighbor_var)),
    _grad_phi_neighbor(_assembly.gradPhiFaceNeighbor(_neighbor_var)),
    _test_neighbor(_neighbor_var.phiFaceNeighbor()),
    _grad_test_neighbor(_neighbor_var.gradPhiFaceNeighbor()),
	_count(_var.count()),
	_work_vector(_count)

{
  addMooseVariableDependency(this->mooseVariable());

  if (!parameters.isParamValid("boundary"))
    mooseError(
        "In order to use an interface kernel, you must specify a boundary where it will live.");

  if (parameters.isParamSetByUser("save_in"))
  {
    if (_save_in_strings.size() != _save_in_var_side.size())
      mooseError("save_in and save_in_var_side must be the same length");
    else
    {
      for (unsigned i = 0; i < _save_in_strings.size(); ++i)
      {
        ArrayMooseVariable * var = &_subproblem.getArrayVariable(_tid, _save_in_strings[i]);

        if (_sys.hasVariable(_save_in_strings[i]))
          mooseError("Trying to use solution variable " + _save_in_strings[i] +
                     " as a save_in variable in " + name());

        if (_save_in_var_side[i] == "m")
        {
          if (var->feType() != _var.feType())
            mooseError(
                "Error in " + name() +
                ". There is a mismatch between the fe_type of the save-in Auxiliary variable "
                "and the fe_type of the the primary side nonlinear "
                "variable this interface kernel object is acting on.");
          _primary_save_in_residual_variables.push_back(var);
        }
        else
        {
          if (var->feType() != _neighbor_var.feType())
            mooseError(
                "Error in " + name() +
                ". There is a mismatch between the fe_type of the save-in Auxiliary variable "
                "and the fe_type of the the secondary side nonlinear "
                "variable this interface kernel object is acting on.");
          _secondary_save_in_residual_variables.push_back(var);
        }

        var->sys().addVariableToZeroOnResidual(_save_in_strings[i]);
        addMooseVariableDependency(var);
      }
    }
  }

  _has_primary_residuals_saved_in = _primary_save_in_residual_variables.size() > 0;
  _has_secondary_residuals_saved_in = _secondary_save_in_residual_variables.size() > 0;

  if (parameters.isParamSetByUser("diag_save_in"))
  {
    if (_diag_save_in_strings.size() != _diag_save_in_var_side.size())
      mooseError("diag_save_in and diag_save_in_var_side must be the same length");
    else
    {
      for (unsigned i = 0; i < _diag_save_in_strings.size(); ++i)
      {
        ArrayMooseVariable * var = &_subproblem.getArrayVariable(_tid, _diag_save_in_strings[i]);

        if (_sys.hasVariable(_diag_save_in_strings[i]))
          mooseError("Trying to use solution variable " + _diag_save_in_strings[i] +
                     " as a save_in variable in " + name());

        if (_diag_save_in_var_side[i] == "m")
        {
          if (var->feType() != _var.feType())
            mooseError(
                "Error in " + name() +
                ". There is a mismatch between the fe_type of the save-in Auxiliary variable "
                "and the fe_type of the the primary side nonlinear "
                "variable this interface kernel object is acting on.");
          _primary_save_in_jacobian_variables.push_back(var);
        }
        else
        {
          if (var->feType() != _neighbor_var.feType())
            mooseError(
                "Error in " + name() +
                ". There is a mismatch between the fe_type of the save-in Auxiliary variable "
                "and the fe_type of the the secondary side nonlinear "
                "variable this interface kernel object is acting on.");
          _secondary_save_in_jacobian_variables.push_back(var);
        }

        var->sys().addVariableToZeroOnJacobian(_diag_save_in_strings[i]);
        addMooseVariableDependency(var);
      }
    }
  }

  _has_primary_jacobians_saved_in = _primary_save_in_jacobian_variables.size() > 0;
  _has_secondary_jacobians_saved_in = _secondary_save_in_jacobian_variables.size() > 0;
}


template <typename T>
void
InterfaceKernelTempl02<T>::computeResidual()
{ 
  if (!_var.activeOnSubdomain(_current_elem->subdomain_id()) ||
      !_neighbor_var.activeOnSubdomain(_neighbor_elem->subdomain_id()))
    return;
    
  precalculateResidual();
  // Compute the residual for this element
  computeElemNeighResidual(Moose::Element);
}

template <typename T>
void
InterfaceKernelTempl02<T>::computeElemNeighResidual(Moose::DGResidualType type)
{
  bool is_elem;
  if (type == Moose::Element)
    is_elem = true;
  else
    is_elem = false;

  const TemplateVariableTestValue & test_space = is_elem ? _test : _test_neighbor;

  if (is_elem)
    prepareVectorTag(_assembly, _var.number());
  else
    prepareVectorTagNeighbor(_assembly, _neighbor_var.number());


  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    initQpResidual(type);
    for (_i = 0; _i < test_space.size(); _i++)
	{
		_work_vector.setZero();
		_work_vector = computeQpResidual(type);
		mooseAssert(_work_vector.size() == _count,
                  "Size of local residual is not equal to the number of array variable compoments");
		_work_vector *= _JxW[_qp] * _coord[_qp];
		_assembly.saveLocalArrayResidual(_local_re, _i, test_space.size(), _work_vector);
	}
  }
	
  accumulateTaggedLocalResidual();

  if (_has_primary_residuals_saved_in && is_elem)
  {
    Threads::spin_mutex::scoped_lock lock(_resid_vars_mutex);
    for (const auto & var : _primary_save_in_residual_variables)
    {
      var->sys().solution().add_vector(_local_re, var->dofIndices());
    }
  }
  else if (_has_secondary_residuals_saved_in && !is_elem)
  {
    Threads::spin_mutex::scoped_lock lock(_resid_vars_mutex);
    for (const auto & var : _secondary_save_in_residual_variables)
      var->sys().solution().add_vector(_local_re, var->dofIndicesNeighbor());
  }
}
template <typename T>
void
InterfaceKernelTempl02<T>::computeJacobian()
{
  if (!_var.activeOnSubdomain(_current_elem->subdomain_id()) ||
      !_neighbor_var.activeOnSubdomain(_neighbor_elem->subdomain_id()))
    return;

  computeElemNeighJacobian(Moose::ElementElement);
  computeElemNeighJacobian(Moose::NeighborNeighbor);
}

// place holders
template <typename T> void InterfaceKernelTempl02<T>::computeOffDiagElemNeighJacobian(Moose::DGJacobianType type, unsigned int jvar) {}
template <typename T> void InterfaceKernelTempl02<T>::computeElementOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) 
{
}
template <typename T> void InterfaceKernelTempl02<T>::computeNeighborOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar) 
{
}

template <typename T>
void
InterfaceKernelTempl02<T>::computeResidualAndJacobian()
{
  computeResidual();

  if (!isImplicit())
    return;

  for (const auto & [ivariable, jvariable] : _fe_problem.couplingEntries(_tid, _sys.number()))
  {
    if (ivariable->isFV())
      continue;

    const unsigned int ivar = ivariable->number();
    const unsigned int jvar = jvariable->number();

    prepareShapes(jvar);
    prepareNeighborShapes(jvar);

    if (_var.number() == ivar)
      computeElementOffDiagJacobian(Moose::ElementElement, jvar);

    if (_neighbor_var.number() == ivar)
      computeNeighborOffDiagJacobian(Moose::NeighborNeighbor, jvar);
  }
}

template <typename T>
void
InterfaceKernelTempl02<T>::computeElemNeighJacobian(Moose::DGJacobianType type)
{
  const TemplateVariableTestValue & test_space =
      (type == Moose::ElementElement || type == Moose::ElementNeighbor) ? _test : _test_neighbor;
  const TemplateVariableTestValue & loc_phi =
      (type == Moose::ElementElement || type == Moose::NeighborElement) ? _phi : _phi_neighbor;

  unsigned int ivar, jvar;

  switch (type)
  {
    case Moose::ElementElement:
      ivar = jvar = _var.number();
      break;
    case Moose::ElementNeighbor:
      ivar = _var.number(), jvar = _neighbor_var.number();
      break;
    case Moose::NeighborElement:
      ivar = _neighbor_var.number(), jvar = _var.number();
      break;
    case Moose::NeighborNeighbor:
      ivar = _neighbor_var.number(), jvar = _neighbor_var.number();
      break;
    default:
      mooseError("Unknown DGJacobianType ", type);
  }

if (_has_primary_jacobians_saved_in && type == Moose::ElementElement)
  {
    auto rows = _local_ke.m();
    DenseVector<Number> diag(rows);
    for (decltype(rows) i = 0; i < rows; i++)
      diag(i) = _local_ke(i, i);

    Threads::spin_mutex::scoped_lock lock(_jacoby_vars_mutex);
    for (const auto & var : _primary_save_in_jacobian_variables)
      var->sys().solution().add_vector(diag, var->dofIndices());
  }
  else if (_has_secondary_jacobians_saved_in && type == Moose::NeighborNeighbor)
  {
    auto rows = _local_ke.m();
    DenseVector<Number> diag(rows);
    for (decltype(rows) i = 0; i < rows; i++)
      diag(i) = _local_ke(i, i);

    Threads::spin_mutex::scoped_lock lock(_jacoby_vars_mutex);
    for (const auto & var : _secondary_save_in_jacobian_variables)
      var->sys().solution().add_vector(diag, var->dofIndicesNeighbor());
  }
}

template class InterfaceKernelTempl02<RealEigenVector>;
