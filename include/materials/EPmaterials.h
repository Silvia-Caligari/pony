#pragma once

#include "Material.h"

class EPmaterials;

template <>
InputParameters validParams<EPmaterials>();

class EPmaterials : public Material
{
public:
  static InputParameters validParams();

  EPmaterials(const InputParameters &parameters);

protected:

   // MeshGeneratorName       _meshGeneratorName;
   // bool              const _hasMeshGenerator;

    RealVectorValue const sigma_i;
    //RealVectorValue const sigma_e;
    
    Real const C_m;
    Real const Chi;
    
    MaterialProperty<Real> &sigma_l;
    MaterialProperty<Real> &sigma_t;
    MaterialProperty<Real> &sigma_n;
    
    
    
    MaterialProperty<RealVectorValue> const &_a_l;
    MaterialProperty<RealVectorValue> const &_a_t;
    MaterialProperty<RealVectorValue> const &_a_n;
    
    MaterialProperty<RealTensorValue> &_K;
    
    MaterialProperty<Real> &time_coefficient;
        
protected:
  virtual void computeQpProperties();

};
