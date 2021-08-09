[GlobalParams]
  global_init_P      = 0.0
  global_init_T      = 1.0
  scaling_factor_var = '1 1 1'
  gravity            = '0 0 0' # Gravity vector

  global_init_V = 0.23 # 3% above steady solution
#  global_init_V = 0.2234 # ~ steady solution
[]

[EOS]
  [./eos]
    type   = PTConstantEOS
    beta   = 0.0
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
  [./pipe1]
    type = PBOneDFluidComponent
    eos  = eos
    position    = '0 0 0'
    orientation = '1 0 0'
    Dh      = 0.1
    length  = 1.0
    n_elems = 100
    A       = 7.85398e-3
  [../]

  [./pump]
    type = PBPump                               # This is a PBPump component
    eos  = eos
    inputs  = 'pipe1(out)'
    outputs = 'pipe2(in)'
    K       = '0.0 0.0'                        # Form loss coefficient at pump inlet and outlet
    Area      = 7.85398e-3                     # Reference pump flow area
    initial_P = -100                              # Initial pressure
    Head      = 25                             # Specified constant pump head
  [../]

  [./pipe2]
    type = PBOneDFluidComponent
    eos  = eos
    position    = '1 0 0'
    orientation = '1 0 0'
    Dh      = 0.1
    length  = 1.0
    n_elems = 100
    A       = 7.85398e-3
  [../]

  [./junction1]
    type = PBSingleJunction
    eos  = eos
    inputs = 'pipe2(out)'
    outputs = 'pipe3(in)'
  [../]

  [./pipe3] # replaced in coupling test methods
    type = PBOneDFluidComponent
    eos  = eos
    position    = '2 0 0'
    orientation = '1 0 0'
    Dh      = 0.1
    length  = 1.0
    n_elems = 100
    A       = 7.85398e-3
  [../]

  [./junction2]
    type = PBSingleJunction
    eos  = eos
    inputs = 'pipe3(out)'
    outputs = 'pipe4(in)'
  [../]

  [./pipe4]
    type = PBOneDFluidComponent
    eos  = eos
    position    = '3 0 0'
    orientation = '1 0 0'
    Dh      = 0.1
    length  = 1.0
    n_elems = 100
    A       = 7.85398e-3
  [../]
 
  [./branch1]
    type  = PBBranch
    eos   = eos
    inputs  = 'pipe4(out)'
    outputs = 'pipe1(in) pipe5(in)'
    K       = '0.0 0.0 0.0'
    Area    = 7.85398e-3
  [../]

  [./pipe5]
    type = PBOneDFluidComponent
    eos  = eos
    position    = '4 0 0'
    orientation = '1 0 0'
    Dh      = 0.1
    length  = 0.1
    n_elems = 10
    A       = 7.85398e-3
  [../]

  [./pbtdv]
    type  = PBTDV
    input = 'pipe5(out)'
    eos  = eos
    p_bc = 0.0
    T_bc = 1.0
  [../]
[]

[Postprocessors]
  [mFlow_pipe1out]
    type = ComponentBoundaryFlow
    input = pipe1(out)
    execute_on = 'INITIAL NONLINEAR'
  [../]
  [mFlow_pipe2out]
    type = ComponentBoundaryFlow
    input = pipe2(out)
    execute_on = 'INITIAL NONLINEAR'
  [../]
  [mFlow_pipe3out]
    type = ComponentBoundaryFlow
    input = pipe3(out)
    execute_on = 'INITIAL NONLINEAR'
  [../]
  [mFlow_pipe4out]
    type = ComponentBoundaryFlow
    input = pipe4(out)
    execute_on = 'INITIAL NONLINEAR'
  [../]
  [mFlow_pipe5out]
    type = ComponentBoundaryFlow
    input = pipe5(out)
    execute_on = 'INITIAL NONLINEAR'
  [../]

  [vel_pipe5out]
    type = ComponentBoundaryVariableValue
    variable = velocity
    input = pipe5(out)
  [../]

  [./p_pipe1in]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe1(in)
  [../]
  [./p_pipe1out]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe1(out)
  [../]
  [./p_pipe2in]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe2(in)
  [../]
  [./p_pipe2out]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe2(out)
  [../]
  [./p_pipe3in]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe3(in)
  [../]
  [./p_pipe3out]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe3(out)
  [../]
  [./p_pipe4in]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe4(in)
  [../]
  [./p_pipe4out]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe4(out)
  [../]
  [./p_pipe5in]
    type = ComponentBoundaryVariableValue
    variable = pressure
    input = pipe5(in)
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

[Executioner]
  type = Transient
  scheme = explicit-euler

  dt    = 0.01                        # Targeted time step size
  dtmin = 0.01                        # The allowed minimum time step size
  num_steps = 10000                    # Max. simulation time steps

  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-4
  nl_max_its = 100
  l_tol      = 1e-5
  l_max_its  = 200

#  type = Steady
#  petsc_options_iname = '-pc_type -ksp_gmres_restart'  # Additional PETSc settings, name list
#  petsc_options_value = 'lu 101'                 # Additional PETSc settings, value list
#
#  nl_rel_tol = 1e-8
#  nl_abs_tol = 1e-7
#  nl_max_its = 50
#  l_tol      = 1e-6
#  l_max_its  = 100

#  nl_rel_tol = 1e-6
#  nl_abs_tol = 1e-3
#  nl_max_its = 50
#  l_tol      = 1e-5
#  l_max_its  = 100

[]

[Outputs]
  print_linear_residuals = false
  perf_graph = true
  [./csv]
    type = CSV
  [../]
  [./exodus]
    type = Exodus
  [../]
  [./console]
    type = Console
  [../]
[]
