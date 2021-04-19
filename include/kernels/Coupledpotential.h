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

// Forward Declaration
class Coupledpotential;

template <>
InputParameters validParams<Coupledpotential>();


class Coupledpotential : public Kernel
{
public:
  static InputParameters validParams();

  Coupledpotential(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  //virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
    
  //unsigned int _w_ind;
  const VariableValue & _v;
  Real _coef;
};
