//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Mynonlinear.h"

registerMooseObject("ponyApp",Mynonlinear);

defineLegacyParams(Mynonlinear);

InputParameters
Mynonlinear::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The nonlinear term f(u) = alpha*u*(u-beta)(gamma-u), with the weak "
                             "form of $(\\phi_i, f(u))$.");
  params.addParam<Real>("beta", 1.0, "uthresh");
  params.addParam<Real>("delta", 1.0, "udepol");
  return params;
}

Mynonlinear::Mynonlinear(const InputParameters & parameters) :
Kernel(parameters),
_nonlinear_coefficient(getMaterialProperty<Real>("nonlinear_coefficient")),
_u_old(_var.slnOld()),
_beta(getParam<Real>("beta")),
_delta(getParam<Real>("delta"))
{}



Real
Mynonlinear::computeQpResidual()
{
  return _nonlinear_coefficient[_qp] * (_u_old[_qp] * (_u_old[_qp] - _beta) * (_delta - _u_old[_qp])) * _test[_i][_qp];
}

