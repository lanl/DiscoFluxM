
#include "ElasticityTensorRotated.h"
#include "RotationTensor.h"

registerMooseObject("TensorMechanicsApp", ElasticityTensorRotated);

InputParameters
ElasticityTensorRotated::validParams()
{
  InputParameters params = ComputeRotatedElasticityTensorBaseTempl<false>::validParams();
  params.addClassDescription("Computes elasticity tensor in rotated configuration.");
  params.addRequiredParam<std::vector<Real>>("C_ijkl", "Stiffness tensor for material before rotation");
  params.addParam<MooseEnum>("fill_method", RankFourTensor::fillMethodEnum() = "symmetric9", "The fill method");
  params.addRequiredParam<UserObjectName>("read_prop_user_object", "Read the euler angles from a file");
  return params;
}

ElasticityTensorRotated::ElasticityTensorRotated(
    const InputParameters & parameters)
  : ComputeRotatedElasticityTensorBaseTempl<false>(parameters),
    _Cijkl(this->template getParam<std::vector<Real>>("C_ijkl"),
           (RankFourTensor::FillMethod)(int)this->template getParam<MooseEnum>("fill_method")),

   
    _read_prop_user_object(isParamValid("read_prop_user_object")
                               ? &getUserObject<PropertyReadFile>("read_prop_user_object")
                               : nullptr),
    _Euler_angles_mat_prop(declareProperty<RealVectorValue>("Euler_angles")),
    _Euler_angles_mat_prop_intf(declareProperty<RealVectorValue>("Euler_angles_intf")),
    _crysrot(declareProperty<RankTwoTensor>(_base_name + "crysrot")),
    _R(_Euler_angles)

{
      // Compute rotation matrix according to the euler angle provided
      RotationTensor R(_Euler_angles);

      // rotate elasticity tensor
      _Cijkl.rotate(R);

}

void
ElasticityTensorRotated::MaterialPropertyInterface()
{
    // get the euler angles for the current element
    _Euler_angles_mat_prop_intf[_qp](0) = _read_prop_user_object->getData(_current_elem, 0);
    _Euler_angles_mat_prop_intf[_qp](1) = _read_prop_user_object->getData(_current_elem, 1);
    _Euler_angles_mat_prop_intf[_qp](2) = _read_prop_user_object->getData(_current_elem, 2); 
}
void
ElasticityTensorRotated::getEulerAngles()
{ 

    // get the euler angles for the current element
    _Euler_angles_mat_prop[_qp](0) = _read_prop_user_object->getData(_current_elem, 0);
    _Euler_angles_mat_prop[_qp](1) = _read_prop_user_object->getData(_current_elem, 1);
    _Euler_angles_mat_prop[_qp](2) = _read_prop_user_object->getData(_current_elem, 2);

  _R.update(_Euler_angles_mat_prop[_qp]);
}

void
ElasticityTensorRotated::initQpStatefulProperties()
{
    getEulerAngles();
    _crysrot[_qp] = _R.transpose();
}

void
ElasticityTensorRotated::computeQpElasticityTensor()
{

    getEulerAngles();
    _crysrot[_qp] = _R.transpose();

  // cmpute the rotated elasticity tensor.
  _elasticity_tensor[_qp] = _Cijkl;
  _elasticity_tensor[_qp].rotate(_crysrot[_qp]);
}
