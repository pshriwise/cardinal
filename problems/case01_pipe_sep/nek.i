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
  exodus = true
[]

[Postprocessors]
  [p_average]
    type = NekSideIntegral
    field = pressure
    boundary = '1'
  []

  [temp_average]
    type = NekSideIntegral
    field = temperature
    boundary = '1'
  []
[]
