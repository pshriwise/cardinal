[Mesh]
  type = NekRSMesh
  boundary = '1 2'
[]

[Problem]
  type = NekRSProblem
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
  [p_interface]
    type = NekSideAverage
    field = pressure
    boundary = '1'
  []

  [p_out]
    type = NekSideAverage
    field = pressure
    boundary = '2'
  []

  [temp_interface]
    type = NekSideAverage
    field = temperature
    boundary = '1'
  []
  [vel_interface]
    type = Receiver
  []
[]
