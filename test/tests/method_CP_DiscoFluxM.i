[Mesh]
 [./MeshFile]
  type = FileMeshGenerator
  file = './Mesh/Mesh_n2_4X4X2.e'
 [../]
[]

[Outputs]
    file_base = Data2_DiscofluxM_
    csv = true
    exodus = true
    console = true
[]

[UserObjects]
  [./prop_read]
    type = PropertyReadFile
    prop_file_name = './euler_ang_file.inp'
    nprop = 3
    read_type = block
    nblock= 2
  [../]
[]

[Variables]
  [./DD_EdgePositive]
    components = 12
  [../]
  [./DD_EdgeNegative]
    components = 12
  [../]
  [./DD_ScrewPositive]
    components = 12
  [../]
  [./DD_ScrewNegative]
    components = 12
  [../]
[]
  
[ICs]
  [./IC_DD_EdgePositive]
    type = ArrayFunctionIC  
    variable = DD_EdgePositive
    function = '1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0'
  [../]  
  [./IC_DD_EdgeNegative]
    type = ArrayFunctionIC  
    variable = DD_EdgeNegative
   function = '1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0'
  [../]
[]

[Kernels]
#inactive = 'Source_DD_EdgePositive Source_DD_EdgeNegative '
  [./dot_DD_EdgePositive]
    type = ArrayTimeDerivative  
    variable = DD_EdgePositive
  [../]
  [./Advection_DD_EdgePositive]
    type = ArrayTransportDislocationCP  
    variable = DD_EdgePositive
    dislocation_character = edge
    dislocation_sign = positive
  [../]
  [./Source_DD_EdgePositive]
    type = ArraySourceDislocationVolume  
    variable = DD_EdgePositive
    dislocation_character = edge
    dislocation_sign = positive
  [../]
  
  [./dot_DD_EdgeNegative]
    type = ArrayTimeDerivative  
    variable = DD_EdgeNegative
  [../]
  [./Advection_DD_EdgeNegative]
    type = ArrayTransportDislocationCP
    variable = DD_EdgeNegative
    dislocation_character = edge
    dislocation_sign = negative
  [../]
  [./Source_DD_EdgeNegative]
    type = ArraySourceDislocationVolume  
    variable = DD_EdgeNegative
    dislocation_character = edge
    dislocation_sign = negative
  [../]
  
  
  [./dot_DD_ScrewPositive]
    type = ArrayTimeDerivative  
    variable = DD_ScrewPositive
  [../]
  
  [./dot_DD_ScrewNegative]
    type = ArrayTimeDerivative  
    variable = DD_ScrewNegative
  [../]
[]

[InterfaceKernels]
  [./Intf_Avd_DD_EdgePositive]
    type = ArrayDislocationTransferAtGrainBoundary
    variable = DD_EdgePositive
    neighbor_var = DD_EdgePositive
    dislocation_character = edge
    dislocation_sign = positive
    boundary = surface_GB
    density_critical = 1.5
    tau_critical = 50
  [../]
  [./Intf_Avd_DD_EdgeNegative]
    type = ArrayDislocationTransferAtGrainBoundary
    variable = DD_EdgeNegative
    neighbor_var = DD_EdgeNegative
    dislocation_character = edge
    dislocation_sign = negative
    boundary = surface_GB
    density_critical = 1.5
    tau_critical = 50
  [../]
  
[]

[BCs]
  [./BC_Loading]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'top bottom '
    function = '-0.1*y*t'
  [../]

  [./BC_X]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'bottom '
    function = 0.0
  [../]
   [./BC_Y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'bottom'
    function = 0.0
  [../]
   [./BC_Z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'bottom '
    function = 0.0
  [../]
  
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules/TensorMechanics/Master/all]
  strain = FINITE
  add_variables = true
  new_system = true
  formulation = total
  incremental = true
[]

[Materials]
  [./elasticity_tensor]
    type = ElasticityTensorRotated
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
    read_prop_user_object = prop_read
  [../]
  [./compute_stress]
    type = StressUpdateCP
    discoflux_model_name = 'CP_DiscoFlux'
  [../]
  [./CP_DiscoFlux]
    type = CrystalPlasticityDiscoflux
    number_slip_systems = 12
    slip_sys_file_name = ./input_slip_sys.inp
    lattice_friction = 15
    Coeff_hardening = 0.85
    slip_increment_tolerance = 2.0e-2
    DD_EdgePositive = DD_EdgePositive
    DD_EdgeNegative = DD_EdgeNegative
    DD_ScrewPositive = DD_ScrewPositive
    DD_ScrewNegative = DD_ScrewNegative
  [../]
  [updated_euler_angle]
    type = ComputeUpdatedEulerAngle  #ComputeUpdatedEulerAngleDF
    radian_to_degree = true
  []
[]

[GlobalParams]
    dislo_density_initial = 1.0e+05
    C_multi = 8.5e-06
    C_trap = 5.5e-03
    C_m_ann = 0.5
    C_im_ann = 0.5
    burgers_vector_mag = 2.52e-07
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  line_search = 'none'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'
  nl_abs_tol = 1e-5 
  nl_max_its = 10
  l_abs_tol = 1e-5
  l_max_its = 20

  dtmax = 0.002
  dtmin = 0.00000001
  end_time = 1
  [./TimeStepper]
    type = ConstantDT
    dt = 0.001
    growth_factor = 1.01
  [../]
[]

[AuxVariables]
  [./pk2_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./E_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./pk2_yy] 
   type = RankTwoAux
   variable = pk2_yy
   rank_two_tensor = stress
   index_j = 1
   index_i = 1
  [../]
  [./E_yy]
   type = RankTwoAux
   variable = E_yy
   rank_two_tensor = total_strain
   index_j = 1
   index_i = 1
  [../]
[]

[Postprocessors]
  [./pk2_yy]
    type = ElementAverageValue
    variable = pk2_yy
  [../]
  [./E_yy]
    type = ElementAverageValue
    variable = E_yy
  [../]
[]

[AuxVariables]
  [./DD_EP_00]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_01]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_02]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_03]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_04]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_05]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_06]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_07]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_08]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_09]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_10]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EP_11]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]


[AuxKernels]
  [./DD_EP_00]
   type = ArrayVariableComponent
   variable = DD_EP_00
   array_variable = DD_EdgePositive  
   component = 0
  [../]
  [./DD_EP_01]
   type = ArrayVariableComponent
   variable = DD_EP_01
   array_variable = DD_EdgePositive  
   component = 1
  [../]
  [./DD_EP_02]
   type = ArrayVariableComponent
   variable = DD_EP_02
   array_variable = DD_EdgePositive  
   component = 2
  [../]
  [./DD_EP_03]
   type = ArrayVariableComponent
   variable = DD_EP_03
   array_variable = DD_EdgePositive  
   component = 3
  [../]
  [./DD_EP_04]
   type = ArrayVariableComponent
   variable = DD_EP_04
   array_variable = DD_EdgePositive  
   component = 4
  [../]
  [./DD_EP_05]
   type = ArrayVariableComponent
   variable = DD_EP_05
   array_variable = DD_EdgePositive  
   component = 5
  [../]
  [./DD_EP_06]
   type = ArrayVariableComponent
   variable = DD_EP_06
   array_variable = DD_EdgePositive  
   component = 6
  [../]
  [./DD_EP_07]
   type = ArrayVariableComponent
   variable = DD_EP_07
   array_variable = DD_EdgePositive  
   component = 7
  [../]
  [./DD_EP_08]
   type = ArrayVariableComponent
   variable = DD_EP_08
   array_variable = DD_EdgePositive  
   component = 8
  [../]
  [./DD_EP_09]
   type = ArrayVariableComponent
   variable = DD_EP_09
   array_variable = DD_EdgePositive  
   component = 9
  [../]
  [./DD_EP_10]
   type = ArrayVariableComponent
   variable = DD_EP_10
   array_variable = DD_EdgePositive  
   component = 10
  [../]
  [./DD_EP_11]
   type = ArrayVariableComponent
   variable = DD_EP_11
   array_variable = DD_EdgePositive  
   component = 11
  [../]
[]
  

[AuxVariables]
  [./DD_EN_00]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_01]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_02]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_03]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_04]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_05]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_06]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_07]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_08]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_09]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_10]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_EN_11]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]


[AuxKernels]
  [./DD_EN_00]
   type = ArrayVariableComponent
   variable = DD_EN_00
   array_variable = DD_EdgeNegative  
   component = 0
  [../]
  [./DD_EN_01]
   type = ArrayVariableComponent
   variable = DD_EN_01
   array_variable = DD_EdgeNegative  
   component = 1
  [../]
  [./DD_EN_02]
   type = ArrayVariableComponent
   variable = DD_EN_02
   array_variable = DD_EdgeNegative  
   component = 2
  [../]
  [./DD_EN_03]
   type = ArrayVariableComponent
   variable = DD_EN_03
   array_variable = DD_EdgeNegative  
   component = 3
  [../]
  [./DD_EN_04]
   type = ArrayVariableComponent
   variable = DD_EN_04
   array_variable = DD_EdgeNegative  
   component = 4
  [../]
  [./DD_EN_05]
   type = ArrayVariableComponent
   variable = DD_EN_05
   array_variable = DD_EdgeNegative  
   component = 5
  [../]
  [./DD_EN_06]
   type = ArrayVariableComponent
   variable = DD_EN_06
   array_variable = DD_EdgeNegative 
   component = 6
  [../]
  [./DD_EN_07]
   type = ArrayVariableComponent
   variable = DD_EN_07
   array_variable = DD_EdgeNegative
   component = 7
  [../]
  [./DD_EN_08]
   type = ArrayVariableComponent
   variable = DD_EN_08
   array_variable = DD_EdgeNegative
   component = 8
  [../]
  [./DD_EN_09]
   type = ArrayVariableComponent
   variable = DD_EN_09
   array_variable = DD_EdgeNegative
   component = 9
  [../]
  [./DD_EN_10]
   type = ArrayVariableComponent
   variable = DD_EN_10
   array_variable = DD_EdgeNegative
   component = 10
  [../]
  [./DD_EN_11]
   type = ArrayVariableComponent
   variable = DD_EN_11
   array_variable = DD_EdgeNegative
   component = 11
  [../]
[]
  


[AuxVariables]
  [./shear_stress_00]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_00]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_00]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]


[AuxKernels]
[./shear_stress_00]
   type = MaterialStdVectorAux
   variable = shear_stress_00
   property = applied_shear_stress  
   index = 0
  [../]
[./back_stress_00]
   type = MaterialStdVectorAux
   variable = back_stress_00
   property = back_stress
   index = 0
  [../]
[./DD_mobile_00]
   type = MaterialStdVectorAux
   variable = DD_mobile_00
   property = dislocation_immobile
   index = 0
  [../]
[]

[AuxVariables]
  [./shear_stress_01]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_01]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_01]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_01]
   type = MaterialStdVectorAux
   variable = shear_stress_01
   property = applied_shear_stress  
   index = 1
  [../]
[./back_stress_01]
   type = MaterialStdVectorAux
   variable = back_stress_01
   property = back_stress
   index = 1
  [../]
[./DD_mobile_01]
   type = MaterialStdVectorAux
   variable = DD_mobile_01
   property = dislocation_immobile
   index = 1
  [../]
[]

[AuxVariables]
  [./shear_stress_02]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_02]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_02]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_02]
   type = MaterialStdVectorAux
   variable = shear_stress_02
   property = applied_shear_stress  
   index = 2
  [../]
[./back_stress_02]
   type = MaterialStdVectorAux
   variable = back_stress_02
   property = back_stress
   index = 2
  [../]
[./DD_mobile_02]
   type = MaterialStdVectorAux
   variable = DD_mobile_02
   property = dislocation_immobile
   index = 2
  [../]
[]

[AuxVariables]
  [./shear_stress_03]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_03]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_03]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_03]
   type = MaterialStdVectorAux
   variable = shear_stress_03
   property = applied_shear_stress  
   index = 3
  [../]
[./back_stress_03]
   type = MaterialStdVectorAux
   variable = back_stress_03
   property = back_stress
   index = 3
  [../]
[./DD_mobile_03]
   type = MaterialStdVectorAux
   variable = DD_mobile_03
   property = dislocation_immobile
   index = 3
  [../]
[]

[AuxVariables]
  [./shear_stress_04]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_04]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_04]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_04]
   type = MaterialStdVectorAux
   variable = shear_stress_04
   property = applied_shear_stress  
   index = 4
  [../]
[./back_stress_04]
   type = MaterialStdVectorAux
   variable = back_stress_04
   property = back_stress
   index = 4
  [../]
[./DD_mobile_04]
   type = MaterialStdVectorAux
   variable = DD_mobile_04
   property = dislocation_immobile
   index = 4
  [../]
[]

[AuxVariables]
  [./shear_stress_05]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_05]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_05]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_05]
   type = MaterialStdVectorAux
   variable = shear_stress_05
   property = applied_shear_stress  
   index = 5
  [../]
[./back_stress_05]
   type = MaterialStdVectorAux
   variable = back_stress_05
   property = back_stress
   index = 5
  [../]
[./DD_mobile_05]
   type = MaterialStdVectorAux
   variable = DD_mobile_05
   property = dislocation_immobile
   index = 5
  [../]
[]

[AuxVariables]
  [./shear_stress_06]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_06]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_06]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_06]
   type = MaterialStdVectorAux
   variable = shear_stress_06
   property = applied_shear_stress  
   index = 6
  [../]
[./back_stress_06]
   type = MaterialStdVectorAux
   variable = back_stress_06
   property = back_stress
   index = 6
  [../]
[./DD_mobile_06]
   type = MaterialStdVectorAux
   variable = DD_mobile_06
   property = dislocation_immobile
   index = 6
  [../]
[]


[AuxVariables]
  [./shear_stress_07]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_07]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_07]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_07]
   type = MaterialStdVectorAux
   variable = shear_stress_07
   property = applied_shear_stress  
   index = 7
  [../]
[./back_stress_07]
   type = MaterialStdVectorAux
   variable = back_stress_07
   property = back_stress
   index = 7
  [../]
[./DD_mobile_07]
   type = MaterialStdVectorAux
   variable = DD_mobile_07
   property = dislocation_immobile
   index = 7
  [../]
[]



[AuxVariables]
  [./shear_stress_08]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_08]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_08]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_08]
   type = MaterialStdVectorAux
   variable = shear_stress_08
   property = applied_shear_stress  
   index = 8
  [../]
[./back_stress_08]
   type = MaterialStdVectorAux
   variable = back_stress_08
   property = back_stress
   index = 8
  [../]
[./DD_mobile_08]
   type = MaterialStdVectorAux
   variable = DD_mobile_08
   property = dislocation_immobile
   index = 8
  [../]
[]


[AuxVariables]
  [./shear_stress_09]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_09]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_09]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_09]
   type = MaterialStdVectorAux
   variable = shear_stress_09
   property = applied_shear_stress  
   index = 9
  [../]
[./back_stress_09]
   type = MaterialStdVectorAux
   variable = back_stress_09
   property = back_stress
   index = 9
  [../]
[./DD_mobile_09]
   type = MaterialStdVectorAux
   variable = DD_mobile_09
   property = dislocation_immobile
   index = 9
  [../]
[]

[AuxVariables]
  [./shear_stress_10]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_10]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_10]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_10]
   type = MaterialStdVectorAux
   variable = shear_stress_10
   property = applied_shear_stress  
   index = 10
  [../]
[./back_stress_10]
   type = MaterialStdVectorAux
   variable = back_stress_10
   property = back_stress
   index = 10
  [../]
[./DD_mobile_10]
   type = MaterialStdVectorAux
   variable = DD_mobile_10
   property = dislocation_immobile
   index = 10
  [../]
[]

[AuxVariables]
  [./shear_stress_11]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./back_stress_11]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./DD_mobile_11]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./shear_stress_11]
   type = MaterialStdVectorAux
   variable = shear_stress_11
   property = applied_shear_stress  
   index = 11
  [../]
[./back_stress_11]
   type = MaterialStdVectorAux
   variable = back_stress_11
   property = back_stress
   index = 11
  [../]
[./DD_mobile_11]
   type = MaterialStdVectorAux
   variable = DD_mobile_11
   property = dislocation_immobile
   index = 11
  [../]
[]


[Materials]
  [./compute_stress_wrapper]
    type = ComputeLagrangianWrappedStress
  [../]
[]

[AuxVariables]
  [./GND_density_00]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_01]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_02]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_03]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_04]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_05]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_06]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_07]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_08]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_09]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_10]
   order = CONSTANT
   family = MONOMIAL
  [../]
  [./GND_density_11]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[AuxKernels]
[./GND_density_00]
   type = MaterialStdVectorAux
   variable = GND_density_00
   property = kappa  
   index = 0
  [../]
[./GND_density_01]
   type = MaterialStdVectorAux
   variable = GND_density_01
   property = kappa  
   index = 1
  [../]
[./GND_density_02]
   type = MaterialStdVectorAux
   variable = GND_density_02
   property = kappa  
   index = 2
  [../]
[./GND_density_03]
   type = MaterialStdVectorAux
   variable = GND_density_03
   property = kappa  
   index = 3
  [../]
[./GND_density_04]
   type = MaterialStdVectorAux
   variable = GND_density_04
   property = kappa  
   index = 4
  [../]
[./GND_density_05]
   type = MaterialStdVectorAux
   variable = GND_density_05
   property = kappa  
   index = 5
  [../]
[./GND_density_06]
   type = MaterialStdVectorAux
   variable = GND_density_06
   property = kappa  
   index = 6
  [../]
[./GND_density_07]
   type = MaterialStdVectorAux
   variable = GND_density_07
   property = kappa  
   index = 7
  [../]
[./GND_density_08]
   type = MaterialStdVectorAux
   variable = GND_density_08
   property = kappa  
   index = 8
  [../]
[./GND_density_09]
   type = MaterialStdVectorAux
   variable = GND_density_09
   property = kappa  
   index = 9
  [../]
[./GND_density_10]
   type = MaterialStdVectorAux
   variable = GND_density_10
   property = kappa  
   index = 10
  [../]
[./GND_density_11]
   type = MaterialStdVectorAux
   variable = GND_density_11
   property = kappa  
   index = 11
  [../]
[]

