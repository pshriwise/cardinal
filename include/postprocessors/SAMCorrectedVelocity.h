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

class SAMCorrectedVelocity;

template <>
InputParameters validParams<SAMCorrectedVelocity>();

/**
 * computes area-corrected velocity from SAM for a boundary
 **/
class SAMCorrectedVelocity : public NekSidePostprocessor
{
public:
  static InputParameters validParams();

  SAMCorrectedVelocity(const InputParameters & parameters);

  virtual Real getValue() override;

protected:
  Real _rhoArea;
  const PostprocessorValue & _SAM_mFlow_inlet;
  };
