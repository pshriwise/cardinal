//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SAMCorrectedVelocity.h"

registerMooseObject("NekApp", SAMCorrectedVelocity);

defineLegacyParams(SAMCorrectedVelocity);

InputParameters
SAMCorrectedVelocity::validParams()
{
  InputParameters params = NekSidePostprocessor::validParams();
  params.addRequiredParam<PostprocessorName>("SAM_mflow_inlet_pp", "Inlet SAM mass flowrate");
  return params;
}

SAMCorrectedVelocity::SAMCorrectedVelocity(const InputParameters & parameters) :
  NekSidePostprocessor(parameters),
  _SAM_mFlow_inlet(getPostprocessorValue("SAM_mflow_inlet_pp"))
{
}

Real
SAMCorrectedVelocity::getValue()
{
  Real _rhoArea = nekrs::rhoArea(_boundary);

  return _SAM_mFlow_inlet/_rhoArea;
}
