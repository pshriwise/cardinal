//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NekSideMassFlowRate.h"

registerMooseObject("CardinalApp", NekSideMassFlowRate);

defineLegacyParams(NekSideMassFlowRate);

InputParameters
NekSideMassFlowRate::validParams()
{
  InputParameters params = NekSidePostprocessor::validParams();
  return params;
}

NekSideMassFlowRate::NekSideMassFlowRate(const InputParameters & parameters) :
  NekSidePostprocessor(parameters)
{
}

Real
NekSideMassFlowRate::getValue()
{
  return nekrs::massFlowrate(_boundary);
}
