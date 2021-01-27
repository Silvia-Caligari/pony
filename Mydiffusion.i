[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
[]

[Variables]
  [./u]
  [../]
[]

[Functions]
  [./parsed_function]
    type = ParsedFunction
    value = 'sin(x)-cos(y/2)'
  [../]
[]

[ICs]
  [./u_ic]
    type = FunctionIC
    variable = 'u'
    function = parsed_function
  [../]
[]

[Materials]
[./diffusion] 
type = FlowAndTransport 
diffusion1 = 1.0 
diffusion2 = 100.0 
[../]
[]

[Kernels]
  [./diff]
    type = MyDiffusion
    variable = u
    [../]
[./timederivative]
    type = Mytimederivative
    variable = u
    first_coefficient = 1
    second_coefficient = 1
  [../]
[]



[Preconditioning]
  [./pre]
    type = SMP 
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'LINEAR'
  start_time = 0.0
  end_time = 1.0
  dt = 0.05
#  petsc_options = '-pc_svd_monitor -ksp_view_pmat'
#  petsc_options_iname = '-pc_type'
#  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]
