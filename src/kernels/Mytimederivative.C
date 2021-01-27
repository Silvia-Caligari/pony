//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Mytimederivative.h"
#include "Transient.h"
#include "MooseVariableScalar.h"


registerMooseObject("ponyApp",Mytimederivative);

defineLegacyParams(Mytimederivative);

InputParameters
Mytimederivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The time derivative operator with the weak form of $(\\psi_i,\\frac{\\partial u_h}{\\partial t})$, where partial derivative is evaluated ad \\frac{_u - _u_old}{dt}");
  params.addParam<Real>("first_coefficient", 1.0, " First Time Coefficient (C_m)");
  params.addParam<Real>("second_coefficient", 1.0, " Second Time Coefficient (\\Chi)");
  return params;
}

Mytimederivative::Mytimederivative(const InputParameters & parameters) :
Kernel(parameters),
//Monodomain coefficient C_m and Chi
_first_coefficient(getParam<Real>("first_coefficient")),
_second_coefficient(getParam<Real>("second_coefficient")),
_u_old(_var.slnOld())
{}



Real
Mytimederivative::computeQpResidual()
{
    Real _time_coefficient = _first_coefficient / _second_coefficient;
    
    _time_coefficient = _time_coefficient / _dt;
    
  return _time_coefficient * (_u[_qp] - _u_old[_qp]) * _test[_i][_qp];
}

Real
Mytimederivative::computeQpJacobian()
{
    Real _time_coefficient = _first_coefficient / _second_coefficient;
    
    _time_coefficient = _time_coefficient / _dt;
    
    return _time_coefficient * _phi[_j][_qp] * _test[_i][_qp];
}
