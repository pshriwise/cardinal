//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NekSideArea.h"

registerMooseObject("NekApp", NekSideArea);

defineLegacyParams(NekSideArea);

InputParameters
NekSideArea::validParams()
{
  InputParameters params = NekSidePostprocessor::validParams();
  return params;
}

NekSideArea::NekSideArea(const InputParameters & parameters) :
  NekSidePostprocessor(parameters)
{
  if (_fixed_mesh)
    _area = nekrs::area(_boundary);
}

Real
NekSideArea::getValue()
{
  Real area = _fixed_mesh ? _area : nekrs::area(_boundary);

  return area;
}
