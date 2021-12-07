#include "CardinalAction.h"

InputParameters
CardinalAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addParam<std::vector<SubdomainName>>("block",
    "The list of block ids (SubdomainID) to which this action will be applied");
  params.addParam<std::vector<BoundaryName>>("boundary",
    "The list of boundary IDs (BoundaryName) to which this action will be applied");
  return params;
}

CardinalAction::CardinalAction(const InputParameters & parameters)
  : Action(parameters),
    _blocks(getParam<std::vector<SubdomainName>>("block")),
    _boundary(getParam<std::vector<BoundaryName>>("boundary"))
{
}

void
CardinalAction::setObjectBlocks(InputParameters & params, const std::vector<SubdomainName> & blocks)
{
  for (const auto & id: blocks)
    params.set<std::vector<SubdomainName>>("block").push_back(Moose::stringify(id));
}

void
CardinalAction::setObjectBoundaries(InputParameters & params, const std::vector<BoundaryName> & boundaries)
{
  for (const auto & id: boundaries)
    params.set<std::vector<BoundaryName>>("boundary").push_back(Moose::stringify(id));
}
