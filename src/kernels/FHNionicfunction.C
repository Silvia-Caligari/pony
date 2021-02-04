//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FHNionicfunction.h"

registerMooseObject("ponyApp",FHNionicfunction);

defineLegacyParams(FHNionicfunction);

InputParameters
FHNionicfunction::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The nonlinear term f(u) = alpha*u*(u-beta)(gamma-u), with the weak "
                             "form of $(\\phi_i, f(u))$.");
  params.addParam<Real>("uthresh", 1.0, "uthresh");
  params.addParam<Real>("udepol", 1.0, "udepol");
  params.addParam<Real>("urest", 0.0, "urest");
  params.addParam<Real>("alpha", 1.0, "alpha");
  params.addParam<bool>("explicit", true , "define nonlinear discretization explicit/implicit.");
  return params;
}

FHNionicfunction::FHNionicfunction(const InputParameters & parameters) :
Kernel(parameters),
_u_old(_var.slnOld()),
_u_thresh(getParam<Real>("uthresh")),
_u_depol(getParam<Real>("udepol")),
_u_rest(getParam<Real>("urest")),
alpha(getParam<Real>("alpha")),
_explicit(getParam<bool>("explicit"))
{}



Real
FHNionicfunction::computeQpResidual()
{
    if(_explicit)
    {
        return alpha * ((_u_old[_qp] - _u_rest) * (_u_old[_qp] - _u_thresh) * (_u_old[_qp] - _u_depol)) * _test[_i][_qp];
    }
    
    else
    {
        return alpha * ((_u[_qp] - _u_rest) * (_u[_qp] - _u_thresh) * (_u[_qp] - _u_depol)) * _test[_i][_qp];
    }

}

Real
FHNionicfunction::computeQpJacobian()
{
    if(_explicit)
        
        return 0;
    
    else
    {
        Real _u_derivative = alpha * ((_u[_qp] - _u_thresh)*(_u[_qp] - _u_depol) + (_u[_qp] - _u_rest)*(_u[_qp] - _u_depol) + (_u[_qp] - _u_rest)*(_u[_qp] - _u_thresh));
        return  _u_derivative * _phi[_j][_qp] * _test[_i][_qp];
    }
        
    
}

