#pragma once

#include "Material.h"

class BeelerReuter;

template <>
InputParameters validParams<BeelerReuter>();

class BeelerReuter : public Material
{
public:
  static InputParameters validParams();

  BeelerReuter(const InputParameters &parameters);

protected:

    const VariableValue &_V_old;
     
    virtual void computeQpProperties();
    
    virtual void initQpStatefulProperties();
    
    MaterialProperty<Real> &_I_ion;
    
    MaterialProperty<Real> &_m;
    
    MaterialProperty<Real> &_h;
    
    MaterialProperty<Real> &_j;
    
    MaterialProperty<Real> &_d;
    
    MaterialProperty<Real> &_f;
    
    MaterialProperty<Real> &_x1;
    
    MaterialProperty<Real> &_Ca_i;
    
    const MaterialProperty<Real> &_m_old;
    
    const MaterialProperty<Real> &_h_old;
    
    const MaterialProperty<Real> &_j_old;
    
    const MaterialProperty<Real> &_d_old;
    
    const MaterialProperty<Real> &_f_old;
    
    const MaterialProperty<Real> &_x1_old;
    
    const MaterialProperty<Real> &_Ca_i_old;
    
    Real V;
    
    Real i_K1;
    
    Real i_x1;
    
    Real I_X1;
    
    Real i_Na;
    
    Real i_s;
    
    Real E_s;
    
    Real alpha_m, alpha_j, alpha_h, alpha_d, alpha_f, alpha_x1;
    
    Real beta_m, beta_h, beta_j, beta_d, beta_f, beta_x1;
    
    const Real G_Na;
    
    const Real G_NaC;
    
    const Real E_Na;
    
    const Real G_s;
    
    bool explicit;
    
    inline Real update_m(Real V, Real m_old, Real dt);
    
    inline Real update_h(Real V, Real h_old, Real dt);
    
    inline Real update_j(Real V, Real j_old, Real dt);
    
    inline Real update_d(Real V, Real d_old, Real dt);
    
    inline Real update_f(Real V, Real f_old, Real dt);
    
    inline Real update_x1(Real V, Real x1_old, Real dt);
    
    inline Real update_Ca(Real Ca_old, Real dt);

};


