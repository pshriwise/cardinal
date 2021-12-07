#pragma once

#include "InitialCondition.h"
#include "FunctionLayeredIntegral.h"

class BulkEnergyConservationIC;
class InputParameters;

template <typename T>
InputParameters validParams();

template <>
InputParameters validParams<BulkEnergyConservationIC>();

/**
 * Applies a temperature initial condition based on bulk energy conservation
 * in a fluid without any losses
 */
class BulkEnergyConservationIC : public InitialCondition
{
public:
  static InputParameters validParams();

  BulkEnergyConservationIC(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual Real value(const Point & p) override;

  /// Cumulative integral of the heat source in the direction of flow
  const FunctionLayeredIntegral & _layered_integral;

  /// Fluid mass flowrate
  const Real & _mdot;

  /// Fluid isobaric specific heat capacity
  const Real & _cp;

  /// Fluid inlet temperature
  const Real & _inlet_T;

  /// Name of postprocessor providing the integral of the heat source
  const PostprocessorName & _pp_name;

  /// Value of the postprocessor providing the integral of the heat source
  const PostprocessorValue & _integral;

  /// Total magnitude of the heat source
  const Real & _magnitude;
};
