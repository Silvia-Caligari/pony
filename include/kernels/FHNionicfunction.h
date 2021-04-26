#pragma once

#include "Kernel.h"

class FHNionicfunction;

template <>
InputParameters validParams<FHNionicfunction>();


/**
 * This kernel implements the non linear term alpha*u*(u-\beta)(\delta-u)
 */
class FHNionicfunction : public Kernel
{
public:
  static InputParameters validParams();

  FHNionicfunction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
    
  MaterialProperty<Real> const &_I_ion;
    
    
};
