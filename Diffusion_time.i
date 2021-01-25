[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5 
[]

[Variables]
  [./u]
  [../]
[]

[Materials]
[./diffusion] 
type = FlowAndTransport 
diffusion1 = 1.0 
diffusion2 = 10.0 
[../]
[]

[Kernels]
  [./diff]
    type = MyDiffusion
    variable = u
    [../]
[./source]
    type = BodyForce
    variable = u
    value = 1
    function = '(0.5*cos(2*pi*t)-x)^2 + (0.5*sin(2*pi*t)-y)^2 < 0.0625'
  [../]
[]

[BCs]
  [./bottom]
    type = DirichletBC
   # preset = false
    variable = u
    boundary = bottom
    value = 0
    [../]
  [./top]
    type = DirichletBC
   # preset = false
    variable = u
    boundary = top
    value = 0
    [../]
  [./left]
    type = DirichletBC
   # preset = false
    variable = u
    boundary = left
    value = 0
  [../]
  [./right]
    type = DirichletBC
   # preset = false
    variable = u
    boundary = right
    value = 0
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

