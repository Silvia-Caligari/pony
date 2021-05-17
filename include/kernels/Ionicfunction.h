#pragma once

#include "Kernel.h"

class Ionicfunction;

template <>
InputParameters validParams<Ionicfunction>();


/**
 * This kernel implements the non linear term alpha*u*(u-\beta)(\delta-u)
 */
class Ionicfunction : public Kernel
{
public:
  static InputParameters validParams();

  Ionicfunction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
    
  MaterialProperty<Real> const &_I_ion;
    
    
};
