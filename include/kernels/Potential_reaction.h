#pragma once

#include "Kernel.h"

class Potential_reaction;

template <>
InputParameters validParams<Potential_reaction>();


/**
 * This kernel implements the non linear term alpha*u*(u-\beta)(\delta-u)
 */
class Potential_reaction : public Kernel
{
public:
  static InputParameters validParams();

  Potential_reaction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
    
  const VariableValue & _u_old;
    
  const Real & _u_thresh;
    
  const Real & _u_depol;

  const Real & _u_rest;
    
  const Real & alpha; 
    
  const bool _explicit;
    
    
};

