#include "BulkEnergyConservationICAction.h"
#include "IntegralPreservingFunctionIC.h"
#include "FEProblem.h"

registerMooseAction("CardinalApp", BulkEnergyConservationICAction, "add_bulk_fluid_temperature_ic");
registerMooseAction("CardinalApp", BulkEnergyConservationICAction, "add_bulk_fluid_temperature_user_object");

InputParameters
BulkEnergyConservationICAction::validParams()
{
  InputParameters params = CardinalAction::validParams();
  params.addParam<VariableName>("variable", "Name of the fluid temperature variable");
  params.addRequiredParam<unsigned int>("num_layers", "Number of layers to use for integrating the heat source");

  MooseEnum directions("x y z");
  params.addRequiredParam<MooseEnum>("flow_direction", directions,
    "The flow direction along which to integrate the heat source");
  params.addParam<bool>("positive_flow_direction", true,
    "Whether the flow is along the positive 'direction' or negative 'direction'");
  params.addParam<Real>("direction_min",
    "Minimum coordinate along 'direction' that bounds the layers");
  params.addParam<Real>("direction_max",
    "Maximum coordinate along 'direction' that bounds the layers");

  params.addRequiredParam<Real>("mass_flowrate", "Mass flowrate of the fluid");
  params.addRequiredParam<Real>("cp", "Fluid isobaric specific heat capacity");
  params.addRequiredRangeCheckedParam<Real>("inlet_T", "inlet_T >= 0.0", "Inlet temperature");
  return params;
}

BulkEnergyConservationICAction::BulkEnergyConservationICAction(const InputParameters & parameters)
  : CardinalAction(parameters),
    _variable(getParam<VariableName>("variable")),
    _mdot(getParam<Real>("mass_flowrate")),
    _cp(getParam<Real>("cp")),
    _inlet_T(getParam<Real>("inlet_T")),
    _num_layers(getParam<unsigned int>("num_layers")),
    _direction(parameters.get<MooseEnum>("flow_direction")),
    _positive_flow_direction(parameters.get<bool>("positive_flow_direction")),
    _has_direction_min(isParamValid("direction_min")),
    _has_direction_max(isParamValid("direction_max")),
    _direction_min(_has_direction_min ? &getParam<Real>("direction_min") : nullptr),
    _direction_max(_has_direction_max ? &getParam<Real>("direction_max") : nullptr)
{
}

void
BulkEnergyConservationICAction::act()
{
  // by the nature of the task dependencies we set, we can get the function name
  // from the prerequisite task IC
  const auto & ic_warehouse = _problem->getInitialConditionWarehouse();

  if (!ic_warehouse.hasActiveObject("cardinal_heat_source_ic"))
    mooseError("To use the 'BulkEnergyConservation' action syntax, you must also have "
      "a 'VolumetricHeatSource' action!\nYour input file should contain a section like:\n\n"
      "  [Cardinal]\n"
      "    [ICs]\n"
      "      type = VolumetricHeatSource\n"
      "      ...\n"
      "    []\n"
      "  []");

  std::shared_ptr<InitialConditionBase> prereq = ic_warehouse.getObject("cardinal_heat_source_ic");
  std::shared_ptr<IntegralPreservingFunctionIC> ic = std::dynamic_pointer_cast<IntegralPreservingFunctionIC>(prereq);

  if (_current_task == "add_bulk_fluid_temperature_ic")
  {
    const std::string ic_type = "BulkEnergyConservationIC";
    InputParameters params = _factory.getValidParams(ic_type);
    params.set<VariableName>("variable") = _variable;
    params.set<UserObjectName>("layered_integral") = "cardinal_heat_source_layered_integral";
    params.set<Real>("mass_flowrate") = _mdot;
    params.set<Real>("cp") = _cp;
    params.set<Real>("inlet_T") = _inlet_T;
    params.set<PostprocessorName>("integral") = "cardinal_heat_source_integral";
    params.set<Real>("magnitude") = ic->magnitude();

    setObjectBlocks(params, _blocks);
    setObjectBoundaries(params, _boundary);

    _problem->addInitialCondition(ic_type, "cardinal_fluid_temp_ic", params);
  }

  if (_current_task == "add_bulk_fluid_temperature_user_object")
  {
    const std::string uo_type = "FunctionLayeredIntegral";
    InputParameters params = _factory.getValidParams(uo_type);
    params.set<FunctionName>("function") = ic->functionName();
    params.set<MooseEnum>("direction") = _direction;
    params.set<bool>("cumulative") = true;
    params.set<bool>("positive_cumulative_direction") = _positive_flow_direction;
    params.set<unsigned int>("num_layers") = _num_layers;

    if (_has_direction_min)
      params.set<Real>("direction_min") = *_direction_min;

    if (_has_direction_max)
       params.set<Real>("direction_max") = *_direction_max;

    setObjectBlocks(params, _blocks);
    setObjectBoundaries(params, _boundary);

    params.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;
    _problem->addUserObject(uo_type, "cardinal_heat_source_layered_integral", params);
  }
}
