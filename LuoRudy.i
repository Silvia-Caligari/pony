[Mesh]
  [gmg]
   type = GeneratedMeshGenerator
   dim = 2
   nx = 64
   ny = 64
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
    value = '-84.00'
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
  [./fibers_base]
     type = EPfibersdirections
     a_l = '1 0 0'
     a_t = '0 1 0'
     block = 0
  [../] 
  #[./fibers_base_subdomain]
     #type = EPfibersdirections
     #a_l = '1 0 0'
     #a_t = '0 1 0'
     #block = 1
  #[../] 
  [./materials_electrophysiology] 
     type = EPmaterials
     sigma_i = '0.001 0.001 0'
     C_m = 1.0 #membrane conductance
     Chi = 1.0 #surface per volume 
     block = 0 
  [../]
  #[./materials_electrophysiology_subdomain] 
     #type = EPmaterials
     #sigma_i = '0.0001 0.001 0'
     #C_m = 1.0 #membrane conductance
     #Chi = 1.0 #surface per volume 
     #block = 1
  #[../]
  [./LRmaterials]
     type = LuoRudy
     potential = u
  [../]
[]

[Kernels]
  [./diff]
     type = EPdiffusion
     variable = u
  [../]
  [./timederivative]
     type = EPtimederivative
     variable = u
     order = 1 #order of time discretization (=1 or =2)
  [../]
  [./Nonlinear]
     type = Ionicfunction
     variable = u
     #explicit = true # =true if explicit or =false if implicit
  [../]
[]

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
    #petsc_options = '-snes_converged_reason'
    #petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_max_it -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol '
    #petsc_options_value = 'gmres asm lu 100 NONZERO 2 1E-14 1E-12'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  start_time = 0.0
  end_time = 200
  dt = 0.01
  petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
  petsc_options_value='  gmres  lu       NONZERO               mumps '
  #petsc_options = '-pc_svd_monitor -ksp_view_pmat'
#  petsc_options_iname = '-pc_type'
#  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]
