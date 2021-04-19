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
  [./w]
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
  [./ic_function2]
    type = ParsedFunction
    value = '0.0*x*y'
  [../]
[]

[ICs]
  [./u_ic]
    type = FunctionIC
    variable = 'u'
    function = ic_function1
  [../]
  [./w_ic]
     type = FunctionIC
     variable = 'w'
     function = ic_function2
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
     type =  EPmaterials
     sigma_i = '0.001 0.001 0'
     C_m = 1.0 #membrane conductance
     Chi = 1.0 #surface per volume 
     block = 0 
  [../]
  #[./materials_electrophysiology_subdomain] 
     #type =  EPmaterials
     #sigma_i = '0.0001 0.001 0'
     #C_m = 1.0 #membrane conductance
     #Chi = 1.0 #surface per volume 
     #block = 1
  #[../]
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
     type = FHNionicfunction
     variable = u
     uthresh = 0.1
     udepol = 1.0
     urest = 0.0
     alpha =  5.0 #ionic function coefficient f(u)=alpha*(u-u_rest)(u-u_thresh)(u-u_depol)
     explicit = true # =true if explicit or =false if implicit
  [../]
  [./coupling1]
     type = CoupledGating
     variable = u
     coupled_variable = w
     coef = 1.0
  [../]
  [./timegating]
     type = EPtimederivative
     variable = w
     order = 1
  [../]
  [./gating]
     type = Gating
     variable = w
  [../] 
  [./coupling2]
     type = Coupledpotential
     variable = w
     coupled_variable = u
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
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  start_time = 0.0
  end_time = 200.0
  dt = 0.25
  #petsc_options_iname=' -ksp_type -pc_type -pc_factor_shift_type -pc_factor_mat_solver_package '
  #petsc_options_value='  preonly   lu       NONZERO               mumps '
  #petsc_options = '-pc_svd_monitor -ksp_view_pmat'
#  petsc_options_iname = '-pc_type'
#  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]
