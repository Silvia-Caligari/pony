#pragma once

#include "Kernel.h"

class Linear_Elasticity;

template <>
InputParameters validParams<Linear_Elasticity>();

class Linear_Elasticity : public Kernel
{
public:
  static InputParameters validParams();

  Linear_Elasticity(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  //virtual Real computeQpJacobian() override;
    
    int _ind_disp_x;
    
    int _ind_disp_y;
    
    int _ind_disp_z;
    
    unsigned int _ind;
    
    MaterialProperty<RealTensorValue> const &_S;
    
    void computebasisTensor(RealTensorValue &T, const RealVectorValue &phi, const int _ind);
    
};
