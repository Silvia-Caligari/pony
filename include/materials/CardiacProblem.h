#pragma once

#include "Material.h"

class CardiacProblem;

template <>
InputParameters validParams<CardiacProblem>();

class CardiacProblem : public Material
{
public:
  static InputParameters validParams();

  CardiacProblem(const InputParameters &parameters);

/*protected:
  virtual void computeQpProperties();*/
  

protected:

   // MeshGeneratorName       _meshGeneratorName;
   // bool              const _hasMeshGenerator;

    Real const _diffusion_1;
    Real const _diffusion_2;
    Real const _C_m;
    Real const _Chi;
    Real const _alpha;

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
    MaterialProperty<Real> &_diffusion;
    MaterialProperty<Real> &_time_coefficient;
    MaterialProperty<Real> &_nonlinear_coefficient;
    //MaterialProperty<RealVectorValue> &_U;
    
protected:
  virtual void computeQpProperties();

};
