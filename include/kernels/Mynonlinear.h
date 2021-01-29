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

class Mynonlinear;

template <>
InputParameters validParams<Mynonlinear>();


/**
 * This kernel implements the non linear term alpha*u*(u-\beta)(\delta-u)
 */
class Mynonlinear : public Kernel
{
public:
  static InputParameters validParams();

  Mynonlinear(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  //virtual Real computeQpJacobian() override;
    
  MaterialProperty<Real> const & _nonlinear_coefficient;
    
  const VariableValue & _u_old;
    
  const Real & _beta;
    
  const Real & _delta;
    
};

