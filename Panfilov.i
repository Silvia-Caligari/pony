[Problem]
   type = FEproblem_Panfilov
   alpha = 0.5
   theta = 0.5 
   k = 8.0
   a = 0.15
   e0 = 0.002
   mu_1 = 0.2
   mu_2 = 0.3
   method = 0 #(method = 0 predictor-corrector, method =1 Newton-Method)
   how_solve = 1
   material_name = 'materials_new'
   kernel_coverage_check = false
[]

[Mesh]
  [gmg]
   type = GeneratedMeshGenerator
   dim = 2
   nx = 100
   ny = 100
   xmax = 1
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
  [./u]
  [../]
[]

#[AuxVariables]
#[./activation_time]
#[../]
#[]

[Functions]
  [./ic_function1]
    type = ParsedFunction
    value = '(x<0.05)*(y<0.05)'
  [../]
[]

[ICs]
  [./u_ic]
    type = FunctionIC
    variable = 'u'
    function = ic_function1
  [../]
[]

[Materials]
  [./materials_new]
     type = EPmaterialsNew
     a_l = '1 0 0'
     a_t = '0 1 0'
     sigma_i = '0.001 0.001 0'
     C_m = 1.0 #membrane conductance
     Chi = 1.0 #surface per volume 
  [../]
[]

#[Kernels]
  #[./diff]
     #type = EPdiffusion
     #variable = u
  #[../]
  #[./timederivative]
     #type = EPtimederivative
     #variable = u
     #order = 2 #order of time discretization (=1 or =2)
  #[../]
  #[./Nonlinear]
     #type = Ionicfunction
     #variable = u
     #explicit = true # =true if explicit or =false if implicit
  #[../]
#[]

#[AuxKernels]
#[./activation]
#type = TimeActivation
#variable = activation_time
#coupled_variable = u
#uthresh = 0.25
#[../]     
#[]

[Preconditioning]
  [./pre]
     type = SMP 
     full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'LINEAR'
  #scheme = 'crank-nicolson' 
  start_time = 0.0
  end_time = 100.0
  dt = 0.05
  petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
  petsc_options_value='  preonly   lu       NONZERO               mumps '
  #petsc_options = '-pc_svd_monitor -ksp_view_pmat'
#  petsc_options_iname = '-pc_type'
#  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]
