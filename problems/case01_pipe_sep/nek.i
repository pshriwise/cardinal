[Mesh]
  type = NekRSMesh
  boundary = '1'

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

