//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NekSideAverageNormVelocity.h"

registerMooseObject("NekApp", NekSideAverageNormVelocity);

defineLegacyParams(NekSideAverageNormVelocity);

InputParameters
NekSideAverageNormVelocity::validParams()
{
  InputParameters params = NekSidePostprocessor::validParams();
  return params;
}

NekSideAverageNormVelocity::NekSideAverageNormVelocity(const InputParameters & parameters) :
  NekSidePostprocessor(parameters)
{
  if (_fixed_mesh)
    _area = nekrs::area(_boundary);
}

Real
NekSideAverageNormVelocity::getValue()
{
  Real mflowrate = nekrs::massFlowrate(_boundary);
  Real area = _fixed_mesh ? _area : nekrs::area(_boundary);
//  Real rho  = nrs->options.getArgs("DENSITY", rho);

  return mflowrate / area;// / rho;
}
