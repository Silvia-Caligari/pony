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

#[AuxVariables]
#[./activation_time]
#[../]
#[]

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
  [./fibers_base]
     type = EPfibersdirections
     a_l = '1 0 0'
     a_t = '0 1 0'
  [../] 
  [./materials_electrophysiology] 
     type =  EPmaterials
     sigma_i = '0.001 0.001 0'
     C_m = 1.0 #membrane conductance
     Chi = 1.0 #surface per volume
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
     type = FHNionicfunction
     variable = u
     uthresh = 0.25
     udepol = 1.0
     urest = 0.0
     alpha =  5.0 #ionic function coefficient f(u)=alpha*(u-u_rest)(u-u_thresh)(u-u_depol)
     explicit = true # =true if explicit or =false if implicit
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
  end_time = 100.0
  dt = 0.002
#  petsc_options = '-pc_svd_monitor -ksp_view_pmat'
#  petsc_options_iname = '-pc_type'
#  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]
