[Mesh]
  type = NekRSMesh
  boundary = '1'
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
    hide = 'flux_integral'
    execute_on = 'final'
  []
[]

[Postprocessors]
  [p_average]
    type = NekSideAverage
    field = pressure
    boundary = '1'
  []

  [temp_average]
    type = NekSideAverage
    field = temperature
    boundary = '1'
  []
  [flux_integral]
    type = Receiver
  []
  [vel_interface]
    type = Receiver
  []
[]
