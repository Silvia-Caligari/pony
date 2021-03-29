//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoupledGating.h"

#include "MooseVariable.h"

registerMooseObject("ponyApp",CoupledGating);

defineLegacyParams(CoupledGating);

InputParameters
CoupledGating::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addClassDescription("Implements gating coupling c*w "
                             " Weak form: $(\\psi_i, \\c w)$.");
  params.addRequiredCoupledVar("coupled_variable", "The coupled gating variable");
  params.addRequiredParam<Real>("coef", "Coefficent ($\\c$).");

  return params;
}

CoupledGating::CoupledGating(const InputParameters & parameters) :
Kernel(parameters),
_w_ind(coupled("coupled_variable")),
_w(coupledValue("coupled_variable")),
_coef(getParam<Real>("coef"))
{
  if (_var.number() == _w_ind)
    mooseError("Coupled variable 'w' needs to be different from 'variable'"
               "consider using Reaction or somethig similar");
}

Real
CoupledGating::computeQpResidual()
{
  return _coef * _w[_qp] * _test[_i][_qp];
}

Real
CoupledGating::computeQpJacobian()
{
  return 0;
}

Real
CoupledGating::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _w_ind)
    return _coef * _phi[_j][_qp] * _test[_i][_qp];
  //else
    return 0.0;
}
