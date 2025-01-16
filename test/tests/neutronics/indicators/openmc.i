[Mesh]
  [sphere]
    type = FileMeshGenerator
    file = ../meshes/sphere.e
  []
  [solid1]
    type = SubdomainIDGenerator
    input = sphere
    subdomain_id = '100'
  []
  [sphereb]
    type = FileMeshGenerator
    file = ../meshes/sphere.e
  []
  [solid2]
    type = SubdomainIDGenerator
    input = sphereb
    subdomain_id = '200'
  []
  [spherec]
    type = FileMeshGenerator
    file = ../meshes/sphere.e
  []
  [solid3]
    type = SubdomainIDGenerator
    input = spherec
    subdomain_id = '300'
  []
  [combine]
    type = CombinerGenerator
    inputs = 'solid1 solid2 solid3'
    positions = '0 0 0
                 0 0 4
                 0 0 8'
  []

  allow_renumbering = false
  parallel_type = replicated
[]

[Adaptivity]
  [Indicators]
    [optical_depth_hmin]
      type = ElementOpticalDepthIndicator
      rxn_rate = 'total'
      h_type = 'min'
    []
    [optical_depth_hmax]
      type = ElementOpticalDepthIndicator
      rxn_rate = 'total'
      h_type = 'max'
    []
    [optical_depth_cuberoot]
      type = ElementOpticalDepthIndicator
      rxn_rate = 'total'
      h_type = 'cube_root'
    []
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  temperature_blocks = '100 200 300'
  initial_properties = xml
  verbose = true
  cell_level = 0
  normalize_by_global_tally = false

  power = 100.0
  source_rate_normalization = 'kappa_fission'

  [Tallies]
    [Mesh]
      type = MeshTally
      score = 'kappa_fission flux total'
    []
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  execute_on = final
  exodus = true
  hide = 'temp cell_instance cell_id kappa_fission flux total'
[]
