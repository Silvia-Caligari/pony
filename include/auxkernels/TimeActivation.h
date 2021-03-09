//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// forward declarations
class TimeActivation;

template <>
InputParameters validParams<TimeActivation>();

/**
 * An AuxKernel that can be used to dteerminate activation time (the potential is greater than the threshold potential) on each qp point.
 */
class TimeActivation : public AuxKernel
{
public:
  static InputParameters validParams();

  TimeActivation(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _coupled_var;

  const VariableValue & _coupled_var_old;
    
  const Real & _u_thresh;
};

