[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1000
  ny = 1000
[]

[Variables]
  [./u]
  [../]
[]

[Functions]
  [./parsed_function]
    type = ParsedFunction
    value = '-85.0 + 95.0*(x<0.05)*(y<0.05)'
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
  [./FHN_EP] 
     type =  EPmaterials
     diffusion_i = 0.17 
     diffusion_e = 0.62
     C_m = 0.01 #membrane conductance
     Chi = 140.0 #surface per volume
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
     uthresh = -57.6
     udepol = 30.0
     urest = -85.0
     alpha =  0.000014 #ionic function coefficient f(u)=alpha*(u-u_rest)(u-u_thresh)(u-u_depol)
     explicit = true # =true if explicit or =false if implicit
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
  solve_type = 'NEWTON'
  start_time = 0.0
  end_time = 1.0
  dt = 0.000125
#  petsc_options = '-pc_svd_monitor -ksp_view_pmat'
#  petsc_options_iname = '-pc_type'
#  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]
