[GlobalParams]
    global_init_P = 0.0                         # Global initial fluid pressure
    global_init_V = 1.0                         # Global initial fluid velocity
    global_init_T = 1.0                         # Global initial temperature for fluid and solid 
    scaling_factor_var = '1 1 1'          # Scaling factors for fluid variables (p, v, T)
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
  [../]
[]

[Components]
  [./pipe1]
    type = PBOneDFluidComponent
    eos = eos                          # The equation-of-state name
    position = '0 0 0'                 # The origin position of this component
    orientation = '0 0 1'              # The orientation of the component
    heat_source = 0                    # Volumetric heat source
    f = 0.64                           # Specified friction coefficient
    Dh = 0.05                          # Equivalent hydraulic diameter
    length = 1.0                       # Length of the component
    n_elems = 1000                      # Number of elements used in discretization
    A = 1.96e-3                        # Area of the One-D fluid component
  [../]

  [./inlet]
    type = PBTDJ
    input = 'pipe1(in)'                # Name of the connected components and the end type
    eos = eos                          # The equation-of-state
    v_bc = 1.0                         # Velocity boundary condition
    T_bc = 1.0
  [../]

  [./outlet]
    type = PressureOutlet
    input = 'pipe1(out) '              # Name of the connected components and the end type
    eos = eos                          # The equation-of-state
    p_bc = '0.0'                     # Pressure boundary condition
  [../]
[]

[Preconditioning]
  active = 'SMP_PJFNK'
  [./SMP_PJFNK]
    type = SMP                         # Single-Matrix Preconditioner
    full = true                        # Using the full set of couplings among all variables
    solve_type = 'PJFNK'               # Using Preconditioned JFNK solution method
    petsc_options_iname = '-pc_type'   # PETSc option, using preconditiong
    petsc_options_value = 'lu'         # PETSc option, using ‘LU’ precondition type in Krylov solve
  [../]
[]

[Executioner]
  type = Transient                    # This is a transient simulation

  dt = 0.02                           # Targeted time step size
  dtmin = 1e-5                        # The allowed minimum time step size

  petsc_options_iname = '-ksp_gmres_restart'  # Additional PETSc settings, name list
  petsc_options_value = '100'                 # Additional PETSc settings, value list

  nl_rel_tol = 1e-7                   # Relative nonlinear tolerance for each Newton solve
  nl_abs_tol = 1e-6                   # Relative nonlinear tolerance for each Newton solve
  nl_max_its = 20                     # Number of nonlinear iterations for each Newton solve

  l_tol = 1e-4                        # Relative linear tolerance for each Krylov solve
  l_max_its = 100                     # Number of linear iterations for each Krylov solve

  start_time = 0.0                    # Physical time at the beginning of the simulation
  num_steps = 200                     # Max. simulation time steps
  end_time = 100.                     # Max. physical time at the end of the simulation

  [./Quadrature]
    type = TRAP                       # Using trapezoid integration rule
    order = FIRST                     # Order of the quadrature
  [../]
[] # close Executioner section

[Postprocessors]
  [./Pin]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe1(in)
  [../]
  [./Pout]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe1(out)
  [../]
[]

[Outputs]
  perf_graph = true                      # Output the performance log
  csv = true
  [./CSV]
    type = CSV
    execute_on = 'timestep_end'
  []
[]
