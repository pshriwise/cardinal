#pragma once

#include "CardinalAction.h"

class VolumetricHeatSourceICAction;

/**
 * Action that automatically sets up a volumetric heat source initial
 * condition as a wrapping of a FunctionElementIntegral and an
 * IntegralPreservingFunctionIC.
 */
class VolumetricHeatSourceICAction : public CardinalAction
{
public:
  VolumetricHeatSourceICAction(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void act();

protected:
  /// Variable name to apply the initial condition to
  const VariableName & _variable;

  /// Functional form for the heat source
  const FunctionName & _function;

  /// Total magnitude of the heat source upon integration
  const Real & _magnitude;
};
