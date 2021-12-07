#pragma once

#include "Action.h"

class CardinalAction;

/**
 * This is a base Action class for derived classes that set up syntax
 * in Cardinal.
 */
class CardinalAction : public Action
{
public:
  CardinalAction(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  /**
   * Set the blocks to which an object created by this action applies
   * @param[in] params input parameters
   * @param[in] blocks block that the object applies to
   */
  virtual void setObjectBlocks(InputParameters & params, const std::vector<SubdomainName> & blocks);
  /**
   * Set the boundaries to which an object created by this action applies
   * @param[in] params input parameters
   * @param[in] boundaries boundaries that the object applies to
   */
  virtual void setObjectBoundaries(InputParameters & params, const std::vector<BoundaryName> & boundaries);

  /// subdomains to which this action applies
  std::vector<SubdomainName> _blocks;

  /// boundaries to which this action applies
  std::vector<BoundaryName> _boundary;
};
