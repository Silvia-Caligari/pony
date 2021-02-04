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

class EPtimederivative;

template <>
InputParameters validParams<EPtimederivative>();


/**
 * This kernel implements the weak form of time derivative:
 * $\partial u \cdot  \phi_i$
 */
class EPtimederivative : public Kernel
{
public:
  static InputParameters validParams();

  EPtimederivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
    
  MaterialProperty<Real> const &time_coefficient;

    
  const int &order; //value = 1 Euler, value = 2 BDF2
    
    
  const VariableValue &_u_old;
    
  const VariableValue &_u_older;
    
    
};
