//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Gating.h"

registerMooseObject("ponyApp",Gating);

defineLegacyParams(Gating);

InputParameters
Gating::validParams()
{
  InputParameters params = Kernel::validParams();
    params.addClassDescription(
        "Implements a simple consuming reaction term with weak form $(\\psi_i, \\coef w_h)$.");
    params.addParam<Real>(
        "coef", 1.0, "The $(\\lambda)$ multiplier, the relative amount consumed per unit time.");
    return params;
  }

  Gating::Gating(const InputParameters & parameters)
    : Kernel(parameters),
      _coef(getParam<Real>("coef"))

  {
  }

  Real
  Gating::computeQpResidual()
  {
      _coef = 0.025;
    return _test[_i][_qp] * _coef * _u[_qp];
  }

  Real
  Gating::computeQpJacobian()
  {
      _coef = 0.025;
    return _test[_i][_qp] * _coef * _phi[_j][_qp];
  }

