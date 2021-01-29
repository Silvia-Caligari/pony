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
    value = '(x<0.05)*(y<0.05)'
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
[./Myproblem] 
type = CardiacProblem 
diffusion1 = 0.001 
diffusion2 = 0.001
C_m = 1 
surf_per_volume = 1
alpha = 18.5150
[../]
#[./timederivative]
#type = FlowAndTransport
#C_m = 1
surf_per_volume = 1
#[../]
#[./nonlinear]
#type = FlowAndTransport
#_alpha = 18.5150
#[../]
[]

[Kernels]
  [./diff]
    type = MyDiffusion
    variable = u
    [../]
[./timederivative]
    type = Mytimederivative
    variable = u

  [../]
[./Nonlinear]
 type = Mynonlinear
variable = u
beta = 0.2383
delta = 1
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
  end_time = 5.0
  dt = 0.05
#  petsc_options = '-pc_svd_monitor -ksp_view_pmat'
#  petsc_options_iname = '-pc_type'
#  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]
