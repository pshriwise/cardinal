//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NekSidePostprocessor.h"
#include "CardinalEnums.h"

class NekSideAverageNormVelocity;

template <>
InputParameters validParams<NekSideAverageNormVelocity>();

/**
 * computes normal velocity going into a boundary,
 * assuming density = 1
 **/
class NekSideAverageNormVelocity : public NekSidePostprocessor
{
public:
  static InputParameters validParams();

  NekSideAverageNormVelocity(const InputParameters & parameters);

  virtual Real getValue() override;

protected:
  // Area by which to normalize
  Real _area;

  // rho
//  Real rho;
  };
