//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "MooseVariableScalar.h"

class Mytimederivative;

template <>
InputParameters validParams<Mytimederivative>();


/**
 * This kernel implements the weak form of time derivative:
 * $\partial u \cdot  \phi_i$
 */
class Mytimederivative : public Kernel
{
public:
  static InputParameters validParams();

  Mytimederivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
    
  MaterialProperty<Real> const & _time_coefficient;
    
  const VariableValue & _u_old; 
    
};
