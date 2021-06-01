//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CorrectedVelocityToSAM.h"

registerMooseObject("NekApp", CorrectedVelocityToSAM);

defineLegacyParams(CorrectedVelocityToSAM);

InputParameters
CorrectedVelocityToSAM::validParams()
{
  InputParameters params = NekSidePostprocessor::validParams();
  params.addRequiredParam<PostprocessorName>("SAM_rhoArea_outlet_pp", "Outlet to SAM rho*Area");
  return params;
}

CorrectedVelocityToSAM::CorrectedVelocityToSAM(const InputParameters & parameters) :
  NekSidePostprocessor(parameters),
  _SAM_rhoArea_outlet(getPostprocessorValue("SAM_rhoArea_outlet_pp"))
{
}

Real
CorrectedVelocityToSAM::getValue()
{
  Real _nek_mflow = nekrs::massFlowrate(_boundary);

  return _nek_mflow/_SAM_rhoArea_outlet;
}
