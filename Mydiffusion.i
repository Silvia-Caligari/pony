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
[]

[BCs]
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
    value = 1
  [../]
[]

[Preconditioning]
  [./pre]
    type = SMP 
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'LINEAR'
#  petsc_options = '-pc_svd_monitor -ksp_view_pmat'
#  petsc_options_iname = '-pc_type'
#  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]
