[GlobalParams]
  global_init_P = 0.0
  global_init_V = 1.0
  global_init_T = 1.0
  scaling_factor_var = '1 1 1'
[]

[EOS]
  [./eos]                                       # EOS name
    type   = PTConstantEOS                      # Using the sodium equation-of-state
    rho_0  = 1.0                                #Dimensionless
    beta   = 1.0                                #Dimensionless
    cp     = 1.0                                #Dimensionless
    h_0    = 1.0                                #Dimensionless
    T_0    = 1.0                                #Dimensionless
    mu     = 5.0e-4                             #Re = 100
    k      = 1.0                                #Dimensionless
#    p_0    = 0.0
  [../]
[]

[Components]

  [./pipe2]
    type = PBOneDFluidComponent
    eos = eos
    position = '0 1 0'
    orientation = '1 0 0'

    WF_user_option = User
    User_defined_WF_parameters = '0 64 -1'
 
    Dh = 0.05                          # Equivalent hydraulic diameter
    length = 1.0                       # Length of the component
    n_elems = 100                      # Number of elements used in discretization
    A = 1.96e-3                        # Area of the One-D fluid component
  [../]

  [./pipe3]
    type = PBOneDFluidComponent
    eos = eos
    position = '1 1 0'
    orientation = '0 -1 0'

    WF_user_option = User
    User_defined_WF_parameters = '0 64 -1'
    Dh = 0.05                          # Equivalent hydraulic diameter
    length = 1.0                       # Length of the component
    n_elems = 100                      # Number of elements used in discretization
    A = 1.96e-3                        # Area of the One-D fluid component
  [../]

  [./pipe4]
    type = PBOneDFluidComponent
    eos = eos
    position = '1 0 0'
    orientation = '-1 0 0'

    WF_user_option = User
    User_defined_WF_parameters = '0 64 -1'
    Dh = 0.05                          # Equivalent hydraulic diameter
    length = 1.0                       # Length of the component
    n_elems = 100                      # Number of elements used in discretization
    A = 1.96e-3                        # Area of the One-D fluid component
  [../]

  [./Branch1]
    type = CoupledPPSTDV
    input = 'pipe2(in)'                # Name of the connected components and the end type
    eos = eos                          # The equation-of-state
    postprocessor_pbc = pres_outlet_interface
    postprocessor_Tbc = temp_outlet_interface
  [../]

  [./Branch2]
    type = PBBranch
    inputs = 'pipe2(out)'
    outputs = 'pipe3(in) pipe5(in)'
    K = '0.0 0.0 0.0'
    Area = 1.96e-3
    initial_P = 0.0
    eos = eos
  [../]

  [./Pump_p]
    type = PBPump                               # This is a PBPump component
    eos = eos
    inputs = 'pipe3(out)'
    outputs = 'pipe4(in)'
    K = '0. 0.'                                 # Form loss coefficient at pump inlet and outlet
    Area = 1.96e-3                              # Reference pump flow area
    initial_P = 30.0                             # Initial pressure
#    Head_fn = pump_startup
    Desired_mass_flow_rate = 1.96e-3
    Head = 25.6
    Response_interval=0.1
    Mass_flow_rate_tolerance = 0.01
  [../]

  [./Branch3]
    type = CoupledPPSTDV
    input = 'pipe4(out)'                # Name of the connected components and the end type
    eos = eos                          # The equation-of-state
    postprocessor_pbc = pres_inlet_interface
    postprocessor_Tbc = temp_inlet_interface
  [../]

  [./pipe5]
    type = PBOneDFluidComponent
    eos = eos
    position = '1 1 0'
    orientation = '1 0 0'
    A = 1.96e-3
    Dh = 0.05
    length = 0.1
    n_elems = 10
  [../]

#  [./TDV_p]
#    type = PBTDV
#    input = 'pipe5(out)'
#    eos = eos
#    p_bc = 100.0
#    T_bc = 1.0
#  [../]
  [./TDV_v]
    type = PBTDJ
    input = 'pipe5(out)'
    eos = eos
    v_bc = 0.0
    T_bc = 1.0
  [../]
[]

[Postprocessors]
  [vel_inlet_interface]
    type = ComponentBoundaryVariableValue
    variable = velocity
    input = pipe4(out)
  [../]
  [nekRS_pres_drop]
    type = Receiver
  [../]
  [negativenekRS_pres_drop]
    type = ScalePostprocessor
    value = nekRS_pres_drop
    scaling_factor = -1.0
  []

  [pres_outlet_interface]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe2(in)
  [../]

  [pres_inlet_interface]
    type = DifferencePostprocessor
    value1 = pres_outlet_interface
    value2 = negativenekRS_pres_drop
  []

  [temp_inlet_interface]
    type = Receiver
  [../]

  [temp_outlet_interface]
    type = Receiver
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
    execute_on = timestep_begin
  []
[]

[Transfers]
  [vel_inlet_interface]
    type = MultiAppPostprocessorTransfer
    direction = to_multiapp
    multi_app = nek
    from_postprocessor = vel_inlet_interface
    to_postprocessor = vel_interface
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


  dt    = 1e-3                        # Targeted time step size
  dtmin = 1e-4                        # The allowed minimum time step size

  petsc_options_iname = '-ksp_gmres_restart'  # Additional PETSc settings, name list
  petsc_options_value = '100'                 # Additional PETSc settings, value list

  nl_rel_tol = 1e-7                   # Relative nonlinear tolerance for each Newton solve
  nl_abs_tol = 1e-6                   # Relative nonlinear tolerance for each Newton solve
  nl_max_its = 20                     # Number of nonlinear iterations for each Newton solve

  l_tol = 1e-4                        # Relative linear tolerance for each Krylov solve
  l_max_its = 100                     # Number of linear iterations for each Krylov solve

  start_time = 0.0                    # Physical time at the beginning of the simulation
#  num_steps = 1000                    # Max. simulation time steps
  end_time = 5.                     # Max. physical time at the end of the simulation


  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]
[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  [./out_displaced]
    type = Exodus
    use_displaced = true
    execute_on = 'initial timestep_end final'
    sequence = false
  [../]
  [./csv]
    type = CSV
  [../]
  [./console]
    type = Console
  [../]
[]
