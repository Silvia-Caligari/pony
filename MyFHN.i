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
  [./FHN_EP] 
     type =  EPmaterials
     diffusion_i = 0.001 
     diffusion_e = 0.001
     C_m = 1 #membrane conductance
     Chi = 1 #surface per volume
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
     order = 2 #order of time discretization (=1 or =2)
  [../]
  [./Nonlinear]
     type = FHNionicfunction
     variable = u
     uthresh = 0.2383
     udepol = 1
     urest = 0
     alpha = 18.32 #ionic function coefficient f(u)=alpha*(u-u_rest)(u-u_thresh)(u-u_depol)
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
  end_time = 0.05
  dt = 0.05
#  petsc_options = '-pc_svd_monitor -ksp_view_pmat'
#  petsc_options_iname = '-pc_type'
#  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]
