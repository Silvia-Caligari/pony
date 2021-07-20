#[Problem]
#type = MyFEproblem
#[]

[Mesh]
  [gmg]
   type = GeneratedMeshGenerator
   dim = 2
   nx = 20
   ny = 2
   xmax = 10
   ymax = 1 
  []
  #[./subdomain]
   #type = ParsedSubdomainMeshGenerator
    #input = gmg
    #combinatorial_geometry = 'x < 1.0 & y < 0.5'
    #block_id = 1
   #[../]
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[BCs]
 [./left_x_dir]
  type = DirichletBC
  variable = disp_x
  boundary = 'left'
  value = 0 
  [../] 
 [./left_y_dir]
  type = DirichletBC
  variable = disp_y
  boundary = 'left'
  value = 0 
  [../] 
 [./y_dir]
  type = NeumannBC
  variable = disp_y
  boundary = 'bottom'
  value = 0.0009
  [../] 
[]

[Materials]
  [./elasticity_tensor]
     type = Isotropic_Elasticity_Tensor
     disp_x = disp_x
     disp_y = disp_y
     mu = 1.0
     lambda = 2.0  
  [../] 

[]

[Kernels]
  [./elasticity_x]
     type = Linear_Elasticity
     variable = disp_x
     disp_x = disp_x
     disp_y = disp_y
  [../]
  [./elasticity_y]
     type = Linear_Elasticity
     variable = disp_y
     disp_x = disp_x
     disp_y = disp_y
  [../]
[]


[Preconditioning]
  [./pre]
     type = FDP 
     full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
  petsc_options_value='  preonly   lu       NONZERO               mumps '
  #petsc_options = '-pc_svd_monitor -ksp_view_pmat'
#  petsc_options_iname = '-pc_type'
#  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]
