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

class NekSideArea;

template <>
InputParameters validParams<NekSideArea>();

/**
 * computes area of a boundary
 **/
class NekSideArea : public NekSidePostprocessor
{
public:
  static InputParameters validParams();

  NekSideArea(const InputParameters & parameters);

  virtual Real getValue() override;

protected:
  Real _area;
  };
