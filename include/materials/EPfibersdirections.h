#pragma once

#include "Material.h"

class EPfibersdirections;

template <>
InputParameters validParams<EPfibersdirections>();

class EPfibersdirections : public Material
{
public:
  static InputParameters validParams();

  EPfibersdirections(const InputParameters &parameters);

protected:

    RealVectorValue a_l;
    RealVectorValue a_t;
    //RealVectorValue a_n;
    
    MaterialProperty<RealVectorValue> &_a_l;
    MaterialProperty<RealVectorValue> &_a_t;
    MaterialProperty<RealVectorValue> &_a_n;
    
    
protected:
  virtual void computeQpProperties();

};
