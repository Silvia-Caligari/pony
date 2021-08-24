#pragma once

#include "Material.h"

class EPmaterialsNew;

template <>
InputParameters validParams<EPmaterialsNew>();

class EPmaterialsNew : public Material
{
public:
  static InputParameters validParams();

  EPmaterialsNew(const InputParameters &parameters);
    
    void getDiffusion(std::vector<Point> const & p , std::vector<RealTensorValue> & diffusion);
    
    void EvaluateDiffusion(Point const & p, RealTensorValue & diffusion);

protected:

   // MeshGeneratorName       _meshGeneratorName;
   // bool              const _hasMeshGenerator;
    RealVectorValue a_l;
    RealVectorValue a_t;
    RealVectorValue const sigma_i;
    //RealVectorValue const sigma_e;
    
    Real const C_m;
    Real const Chi;
    
    /*MaterialProperty<Real> &sigma_l;
    MaterialProperty<Real> &sigma_t;
    MaterialProperty<Real> &sigma_n;
    
    
    
    MaterialProperty<RealVectorValue> const &_a_l;
    MaterialProperty<RealVectorValue> const &_a_t;
    MaterialProperty<RealVectorValue> const &_a_n;*/
    
    MaterialProperty<RealTensorValue> &_K;
    
    MaterialProperty<Real> &time_coefficient;
        
protected:
  virtual void computeQpProperties();

};
