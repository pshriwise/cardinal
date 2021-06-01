//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NekMassFluxWeightedSideIntegral.h"

class NekSideMassFlowRate;

template <>
InputParameters validParams<NekSideMassFlowRate>();

/**
 * mass flow rate from nekRS at a boundary
 */
class NekSideMassFlowRate : public NekSidePostprocessor
{
public:
  static InputParameters validParams();

  NekSideMassFlowRate(const InputParameters & parameters);

  virtual Real getValue() override;
};

