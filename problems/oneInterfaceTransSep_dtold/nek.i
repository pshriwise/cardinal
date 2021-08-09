[Mesh]
  type = NekRSMesh
  boundary = '1 2'
[]

[Problem]
  type = NekRSProblem1d
  SAMtoNek_interface = true # specify that there is coupling interface from SAM to nekRS
[]

[Executioner]
  type = Transient
  timestep_tolerance = 1e-10

  [./TimeStepper]
    type = NekTimeStepper
  [../]
[]

[Outputs]
  [out]
    type = CSV
    execute_on = 'final'
  []
[]

[Postprocessors]
  [SAM_mflow_inlet_interface] # required for SAM-nekRS coupling interface
    type = Receiver
  []

  [nekRS_mflow_inlet] # check inlet nekRS mass flow rate = mass flwo rate from SAM
    type = NekSideMassFlowRate
    boundary = '1'
  []

  [nekRS_mflow_outlet] # check nekRS outlet mass flow rate = inlet mass flow rate
    type = NekSideMassFlowRate
    boundary = '2'
  []

  [nekRS_pres_drop] # the oulet pressure in the nekRS simulation is always zero so the pressure drop is simply the inlet pressure 
    type = NekSideAverage
    field = pressure
    boundary = '1'
  []

  [temp_inlet_interface]
    type = NekSideAverage
    field = temperature
    boundary = '1'
  []

  [temp_outlet_interface]
    type = NekSideAverage
    field = temperature
    boundary = '2'
  []

[]
