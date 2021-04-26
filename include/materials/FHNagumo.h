#pragma once

#include "Material.h"

class FHNagumo;

template <>
InputParameters validParams<FHNagumo>();

class FHNagumo : public Material
{
public:
  static InputParameters validParams();

  FHNagumo(const InputParameters &parameters);

protected:

    const VariableValue &_V_old;
     
    virtual void computeQpProperties();
    
    virtual void initQpStatefulProperties();
    
    const Real &V_rest;
    
    const Real &V_thresh;
    
    const Real &V_depol;
    
    const Real &alpha;
    
    const Real &beta;
    
    const Real &delta;
    
    const Real &gamma;
    
    Real V;
    Real w;
    
    MaterialProperty<Real> &_I_ion;
    
    MaterialProperty<Real> &_w;
    
    const MaterialProperty<Real> &_w_old;

};

