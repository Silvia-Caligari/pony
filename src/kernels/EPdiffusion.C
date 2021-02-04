//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EPdiffusion.h"

registerMooseObject("ponyApp",EPdiffusion);

defineLegacyParams(EPdiffusion);

InputParameters
EPdiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  return params;
}

EPdiffusion::EPdiffusion(const InputParameters & parameters) :
Kernel(parameters),
diffusion(getMaterialProperty<Real>("diffusionProperty"))
{}



Real
EPdiffusion::computeQpResidual()
{
  return diffusion[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
EPdiffusion::computeQpJacobian()
{
  return diffusion[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
