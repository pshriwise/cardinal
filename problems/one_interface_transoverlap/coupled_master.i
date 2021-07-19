[GlobalParams]
  global_init_P = 0.0
  global_init_V = 0.001
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
    rho_0  = 1.0  #9.9756100e02
    mu     = 1e-4 #8.8871000e-04
    k      = 1.0
    p_0    = 0.0
  [../]
[]

[Components]
  [./inlet]
    type  = PBTDV
    input = 'pipe1(in)'
    eos  = eos
    p_bc = 0.32
    T_bc = 1.0
  [../]

  [./pipe1]
    type = PBOneDFluidComponent
    eos  = eos
    position    = '0 0 0'
    orientation = '1 0 0'
    Dh      = 0.1
    length  = 0.666666
    n_elems = 50
    A       = 7.85398e-3

    WF_user_option = User
    User_defined_WF_parameters = '0 64 -1'
  [../]

  [./junction1]
    type = PBSingleJunction
    eos  = eos
    inputs  = 'pipe1(out)'
    outputs = 'pipe2(in)'
  [../]

  [./pipe2]
    type = PBOneDFluidComponent
    eos  = eos
    position    = '0.3333 0 0'
    orientation = '1 0 0'
    Dh      = 0.1
    length  = 0.333333
    n_elems = 1
    A       = 7.85398e-3

    overlap_coupled = True
  [../]

  [./outlet]
    type  = PBTDV
    input = 'pipe2(out)'
    eos  = eos
    p_bc = 0.0
    T_bc = 1.0
  [../]
[]

[Postprocessors]
  [vel_interface]
    type = ComponentBoundaryVariableValue
    variable = velocity
    input = pipe1(out)
    execute_on = 'INITIAL NONLINEAR TIMESTEP_BEGIN TIMESTEP_END'
  [../]

  [SAM_mFlow]
    type = ComponentBoundaryFlow
    input = pipe1(out)
    execute_on = 'INITIAL NONLINEAR'
  [../]

  [nekRS_pres_drop]
    type = Receiver
    execute_on = 'NONLINEAR TIMESTEP_BEGIN TIMESTEP_END'
  [../]

  [nekRS_posGradPres] # required for overlap coupling in SAM
    type = ScalePostprocessor
    value = nekRS_pres_drop
    scaling_factor = 3.0 # one divided by length of coupled section
    execute_on = 'NONLINEAR TIMESTEP_BEGIN TIMESTEP_END'
  [../]

  [nekRS_avgUvolNew]
    type = Receiver
    default = 0.001 # initial condition velocity
    execute_on = 'TIMESTEP_BEGIN'
  [../]

  [nekRS_avgUvolOld]
    type = ScalePostprocessor
    value = nekRS_avgUvolNew
    scaling_factor = 1.0
    execute_on = 'TIMESTEP_END' # want to update this after calculation f_cfd
  [../]

  [nekRS_UNewminusUold]
    type = DifferencePostprocessor
    value1 = nekRS_avgUvolNew
    value2 = nekRS_avgUvolOld
    execute_on = 'NONLINEAR TIMESTEP_BEGIN TIMESTEP_END'
  [../]

  [nekRS_accelGradPres]
    type = ScalePostprocessor
    value = nekRS_UNewminusUold
    scaling_factor = 20 # rho/dt (used in SAM) = 1/5e-2
    execute_on = 'NONLINEAR TIMESTEP_BEGIN TIMESTEP_END'
  [../]

  [coupled_nekdP]
    type = DifferencePostprocessor
    value1 = nekRS_posGradPres
    value2 = nekRS_accelGradPres
    execute_on = 'NONLINEAR TIMESTEP_BEGIN TIMESTEP_END'
  [../]

#  [coupled_nekdP] # required for overlap coupling in SAM
#    type = ScalePostprocessor
#    value = nekRS_pres_drop
#    scaling_factor = 3.0 # one divided by length of coupled section
#    execute_on = 'NONLINEAR TIMESTEP_BEGIN TIMESTEP_END'
#  [../]

  [temp_inlet_interface]
    type = Receiver
  [../]

  [temp_outlet_interface]
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
    type = FDP
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
    execute_on = timestep_begin
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
    to_postprocessor   = nekRS_pres_drop
  [../]

  [nekRS_avgUvol]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    multi_app = nek
    reduction_type = average
    from_postprocessor = nekRS_avgUvol
    to_postprocessor   = nekRS_avgUvolNew
  [../]

  [temp_inlet_interface]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    multi_app = nek
    reduction_type = average
    from_postprocessor = temp_inlet_interface
    to_postprocessor = temp_inlet_interface
  [../]

  [temp_outlet_interface]
    type = MultiAppPostprocessorTransfer
    direction = from_multiapp
    multi_app = nek
    reduction_type = average
    from_postprocessor = temp_outlet_interface
    to_postprocessor = temp_outlet_interface
  [../]
[]

[Executioner]
  type = Transient

  dt    = 5e-2                        # Targeted time step size
  dtmin = 5e-2                        # The allowed minimum time step size
  num_steps = 2

  petsc_options_iname = '-ksp_gmres_restart'  # Additional PETSc settings, name list
  petsc_options_value = '200'                 # Additional PETSc settings, value list

  nl_rel_tol = 1e-6                   # Relative nonlinear tolerance for each Newton solve
  nl_abs_tol = 1e-5                   # Relative nonlinear tolerance for each Newton solve
  nl_max_its = 200                     # Number of nonlinear iterations for each Newton solve

  l_tol = 1e-4                        # Relative linear tolerance for each Krylov solve
  l_max_its = 200                     # Number of linear iterations for each Krylov solve
[]

[Outputs]
  print_linear_residuals = true
  perf_graph = true
  [./csv]
    type = CSV
  [../]
  [./console]
    type = Console
  [../]
[]
