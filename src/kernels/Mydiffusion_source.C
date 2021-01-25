//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Mydiffusion_source.h"

registerMooseObject("ponyApp",Mydiffusion_source);

defineLegacyParams(Mydiffusion_source);

InputParameters
Mydiffusion_source::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$ with source $(\\psi_i, -f)$");
  params.addParam<Real>("constant", 1.0, "Constant source");
    
    return params;
}

Mydiffusion_source::Mydiffusion_source(const InputParameters & parameters) :
Kernel(parameters),
_diffusion(getMaterialProperty<Real>("diffusionProperty")),
_scale(getParam<Real>("constant"))
{}



Real
Mydiffusion_source::computeQpResidual()
{
  return _diffusion[_qp] * _grad_u[_qp] * _grad_test[_i][_qp] + _test[_i][_qp] * -_scale;
}


//Jacobian source equal to zero.

Real
Mydiffusion_source::computeQpJacobian()
{
  return _diffusion[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
