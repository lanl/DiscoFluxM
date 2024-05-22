
#pragma once

#include "ComputeRotatedElasticityTensorBase.h"
#include "PropertyReadFile.h"
#include "RankTwoTensor.h"
#include "RotationTensor.h"

/**
 * This class computes the elasticity tensor and apply the rotation corresponding o 4'th order tensor
 */
class ElasticityTensorRotated : public ComputeRotatedElasticityTensorBaseTempl<false>
{
public:
  static InputParameters validParams();

  ElasticityTensorRotated(const InputParameters & parameters);

  // modulo operator for symmetry operation
  inline int op_modulo(int a, int b) {
  const int result = a % b;
  return result >= 0 ? result : result + b;
}
protected:
  virtual void initQpStatefulProperties() override;
  virtual void getEulerAngles();
  virtual void computeQpElasticityTensor() override;
  virtual void MaterialPropertyInterface() override;

// The material property that this class will compute
  RankFourTensor _Cijkl;

  const PropertyReadFile * const _read_prop_user_object;
  MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop;
  MaterialProperty<RealVectorValue> & _Euler_angles_mat_prop_intf; // for the both sides of interface/GB
  MaterialProperty<RankTwoTensor> & _crysrot;
  RotationTensor _R;
  bool _user_provided_rotation_matrix;

  using ComputeRotatedElasticityTensorBaseTempl<false>::isParamValid;
  using ComputeRotatedElasticityTensorBaseTempl<false>::_elasticity_tensor_name;
  using ComputeRotatedElasticityTensorBaseTempl<false>::_Euler_angles;
  using ComputeRotatedElasticityTensorBaseTempl<false>::_elasticity_tensor;
  using ComputeRotatedElasticityTensorBaseTempl<false>::_qp;
  using ComputeRotatedElasticityTensorBaseTempl<false>::issueGuarantee;
  using ComputeRotatedElasticityTensorBaseTempl<false>::_rotation_matrix;
};


