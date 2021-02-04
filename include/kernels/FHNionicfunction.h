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

class FHNionicfunction;

template <>
InputParameters validParams<FHNionicfunction>();


/**
 * This kernel implements the non linear term alpha*u*(u-\beta)(\delta-u)
 */
class FHNionicfunction : public Kernel
{
public:
  static InputParameters validParams();

  FHNionicfunction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
    
  const VariableValue & _u_old;
    
  const Real & _u_thresh;
    
  const Real & _u_depol;

  const Real & _u_rest;
    
  const Real & alpha; 
    
  const bool _explicit;
    
    
};

