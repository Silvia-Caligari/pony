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

    Real const diffusion_i;
    Real const diffusion_e;
    Real const C_m;
    Real const Chi;
    

    //Real const _poroInput;
    //Real _condFracture;
    
    //Real _dim;
    //RealTensorValue _id;
    //RealVectorValue _u_elem;

    //bool const _isPressureValid;
    //bool const _conservativeScheme;

    //const VariableGradient &_gradP;
  
    //MaterialProperty<Real> &_poro;
    //MaterialProperty<RealTensorValue> &_K;
    MaterialProperty<Real> &diffusion;
    MaterialProperty<Real> &time_coefficient;
    //MaterialProperty<RealVectorValue> &_U;
    
protected:
  virtual void computeQpProperties();

};
