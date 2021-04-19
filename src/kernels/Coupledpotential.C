//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Coupledpotential.h"

#include "MooseVariable.h"

registerMooseObject("ponyApp",Coupledpotential);

defineLegacyParams(Coupledpotential);

InputParameters
Coupledpotential::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addClassDescription("Implements -coef*u "
                             " Weak form: $(\\psi_i, -\\c u)$.");
  params.addRequiredCoupledVar("coupled_variable", "The coupled variable");
  params.addParam<Real>("coef",1.0, "Coefficent ($\\c$).");

  return params;
}

Coupledpotential::Coupledpotential(const InputParameters & parameters) :
Kernel(parameters),
_v(coupledValueOld("coupled_variable")),
_coef(getParam<Real>("coef"))
{}

Real
Coupledpotential::computeQpResidual()
{
    _coef = 0.1;
  return -_coef * _v[_qp] * _test[_i][_qp];
}

Real
Coupledpotential::computeQpJacobian()
{
  return 0;
}

/*Real
CoupledGating::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _w_ind)
    return -_coef * _phi[_j][_qp] * _test[_i][_qp];
  //else
    return 0.0;
}*/
