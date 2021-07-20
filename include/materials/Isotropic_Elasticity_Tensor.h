#pragma once

#include "Material.h"

class Isotropic_Elasticity_Tensor;

template <>
InputParameters validParams<Isotropic_Elasticity_Tensor>();

class Isotropic_Elasticity_Tensor : public Material
{
public:
  static InputParameters validParams();

  Isotropic_Elasticity_Tensor(const InputParameters &parameters);

protected:

    const VariableGradient &_grad_disp_x;
    
    const VariableGradient &_grad_disp_y;
    
    const VariableGradient &_grad_disp_z;
    
    const Real mu;
    
    const Real lambda;
    
    RealTensorValue _U;
    
    RealTensorValue _id;
    
    RealTensorValue _E;
    
    MaterialProperty<RealTensorValue> &_S;
    
    virtual void computeQpProperties();
    
    //virtual void initQpStatefulProperties();
        


};
