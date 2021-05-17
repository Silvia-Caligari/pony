#pragma once

#include "Material.h"

class LuoRudy;

template <>
InputParameters validParams<LuoRudy>();

class LuoRudy : public Material
{
public:
  static InputParameters validParams();

  LuoRudy(const InputParameters &parameters);

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
    
    MaterialProperty<Real> &_x;
    
    MaterialProperty<Real> &_Ca_i;
    
    const MaterialProperty<Real> &_m_old;
    
    const MaterialProperty<Real> &_h_old;
    
    const MaterialProperty<Real> &_j_old;
    
    const MaterialProperty<Real> &_d_old;
    
    const MaterialProperty<Real> &_f_old;
    
    const MaterialProperty<Real> &_x_old;
    
    const MaterialProperty<Real> &_Ca_i_old;
    
    Real V;
    
    Real i_K1;
    
    Real i_K;
    
    Real i_KP;
    
    Real i_Na;
    
    Real i_si;
    
    Real i_b;
    
    Real alpha_m, alpha_j, alpha_h, alpha_d, alpha_f, alpha_x, alpha_k1;
    
    Real beta_m, beta_h, beta_j, beta_d, beta_f, beta_x, beta_k1;
    
    Real k1_inf, Xi, KP;
    
    Real K_o, K_i, Na_o, Na_i, PRNaK, E_Na, E_si, E_K, Ek1, EKP, GK, GK1;
    
    bool _explicit;
    
    inline Real update_m(Real V, Real m_old, Real dt);
    
    inline Real update_h(Real V, Real h_old, Real dt);
    
    inline Real update_j(Real V, Real j_old, Real dt);
    
    inline Real update_d(Real V, Real d_old, Real dt);
    
    inline Real update_f(Real V, Real f_old, Real dt);
    
    inline Real update_x(Real V, Real x1_old, Real dt);
    
    inline Real update_Ca(Real Ca_old, Real dt, Real I_s);

};


