[Mesh]
  type = NekRSMesh
  boundary = '1 2'
[]

[Problem]
  type = NekRSProblem
#  transfer_in = synchronization_in
[]

[Executioner]
  type = Transient

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

  [vel_interface]
    type = Receiver
  []

  [nekRS_pres_drop]           #The pressure drop in the nekRS simulation is always zero so the pressure drop is simply the inlet pressure 
    type = NekSideAverage
    field = pressure
    boundary = '1'
  []

#  [pres_outlet_interface]
#    type = Receiver
#  []

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
