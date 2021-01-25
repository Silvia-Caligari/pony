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

class Mydiffusion_source;

template <>
InputParameters validParams<Mydiffusion_source>();


/**
 * This kernel implements the Laplacian operator with source:
 * $\nabla u \cdot \nabla \phi_i + (\\psi_i, -f)$
 */
class Mydiffusion_source : public Kernel
{
public:
  static InputParameters validParams();

  Mydiffusion_source(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
    
  MaterialProperty<Real> const & _diffusion;
  
  const Real & _scale;
    
};

