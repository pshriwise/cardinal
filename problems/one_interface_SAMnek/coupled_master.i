[GlobalParams]
  global_init_P = 0.0
  global_init_V = 2.5527051E-01
  global_init_T = 1.0
  scaling_factor_var = '1 1 1'
[]

[EOS]
  [./eos]                                       
    type   = PTConstantEOS                      
    beta   = 1.0                                
    cp     = 1.0                                
    h_0    = 1.0                                
    T_0    = 1.0                                
    rho_0  = 9.9756100e02
    mu     = 8.8871000e-04
    k      = 1.0                                
    p_0    = 0.0
  [../]
[]

[Components]
  v_inlet = 2.5527051E-01

  [./pipe1]
    type = PBOneDFluidComponent
    eos = eos
    position = '0 0 0'
    orientation = '1 0 0'
    Dh = 0.1                           
    length = 0.9                       
    n_elems = 100                     
    A = 7.85398e-3
  [../]

  [./Branch1]
    type = CoupledPPSTDV
    input = 'pipe1(out)'               
    eos = eos                          
    postprocessor_pbc = nekRS_pres_drop
    postprocessor_Tbc = temp_inlet_interface
  [../]

  [./inlet]
    type = PBTDJ
    input = 'pipe1(in)'
    eos = eos
    v_bc = ${v_inlet}
    T_bc = 1.0
  [../]
[]


[Postprocessors]
  [vel_interface]
    type = ComponentBoundaryVariableValue
    variable = velocity
    input = pipe1(out)
    execute_on = 'INITIAL NONLINEAR'
  [../]

  [SAM_mFlow]
    type = ComponentBoundaryFlow
    input = pipe1(out)
    execute_on = 'INITIAL NONLINEAR'
  [../]

  [SAM_area_pp]
    type = ComponentBoundaryArea
    input = pipe1(out)
    execute_on = 'INITIAL'
  [../]

  [nekRS_pres_drop]
    type = Receiver
  [../]

  [temp_inlet_interface]
    type = Receiver
  [../]

  [./p_inlet]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe1(in)
  [../]
[]



[Preconditioning]
  active = 'SMP_PJFNK'
  [./SMP_PJFNK]
    type = SMP
    full = true
    solve_type = 'PJFNK'
    petsc_options_iname = '-pc_type -ksp_gmres_restart'
    petsc_options_value = 'lu 101'
  [../]
[]

[MultiApps]
  [nek]
    type = TransientMultiApp
    app_type = NekApp
    input_files = 'nek.i'
    sub_cycling = true
    execute_on = timestep_end
  []
[]

[Transfers]
  [vel_inlet_interface]
    type = MultiAppPostprocessorTransfer
    direction = to_multiapp
    multi_app = nek
    from_postprocessor = vel_interface
    to_postprocessor = vel_interface
  []

  [SAM_mflow_transfer]
    type = MultiAppPostprocessorTransfer
    direction = to_multiapp
    multi_app = nek
    from_postprocessor = SAM_mFlow
    to_postprocessor = SAM_mflow_inlet_interface
  []

  [nekRS_pres_drop]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    multi_app = nek
    reduction_type = average
    from_postprocessor = nekRS_pres_drop
    to_postprocessor = nekRS_pres_drop
  [../]

  [temp_inlet_interface]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    multi_app = nek
    reduction_type = average
    from_postprocessor = temp_inlet_interface
    to_postprocessor = temp_inlet_interface
  [../]
[]

[Executioner]
  type = Transient

  dt    = 1e-1                        # Targeted time step size
  dtmin = 1e-1                        # The allowed minimum time step size
  num_steps = 150 

  petsc_options_iname = '-ksp_gmres_restart'  # Additional PETSc settings, name list
  petsc_options_value = '100'                 # Additional PETSc settings, value list

  nl_rel_tol = 1e-7                   # Relative nonlinear tolerance for each Newton solve
  nl_abs_tol = 1e-6                   # Relative nonlinear tolerance for each Newton solve
  nl_max_its = 20                     # Number of nonlinear iterations for each Newton solve

  l_tol = 1e-4                        # Relative linear tolerance for each Krylov solve
  l_max_its = 100                     # Number of linear iterations for each Krylov solve


  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]
[]


[Outputs]
  print_linear_residuals = false
  perf_graph = true
  [./csv]
    type = CSV
  [../]
  [./console]
    type = Console
  [../]
[]
