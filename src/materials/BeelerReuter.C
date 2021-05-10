#include "BeelerReuter.h"

registerMooseObject("ponyApp",BeelerReuter);

defineLegacyParams(BeelerReuter);

InputParameters BeelerReuter::validParams()
{
	InputParameters params = Material::validParams();

	
    params.addRequiredCoupledVar("potential", "forgot membrane potential FHN.");
    return params;
}

BeelerReuter::BeelerReuter(const InputParameters &parameters) :
Material(parameters),
_V_old(coupledValueOld("potential")),
_I_ion(declareProperty<Real>("I_ion")),
_m(declareProperty<Real>("m")),
_h(declareProperty<Real>("h")),
_j(declareProperty<Real>("j")),
_d(declareProperty<Real>("d")),
_f(declareProperty<Real>("f")),
_x1(declareProperty<Real>("x1")),
_Ca_i(declareProperty<Real>("Ca_i")),
_m_old(getMaterialPropertyOld<Real>("m")),
_h_old(getMaterialPropertyOld<Real>("h")),
_j_old(getMaterialPropertyOld<Real>("j")),
_d_old(getMaterialPropertyOld<Real>("d")),
_f_old(getMaterialPropertyOld<Real>("f")),
_x1_old(getMaterialPropertyOld<Real>("x1")),
_Ca_i_old(getMaterialPropertyOld<Real>("Ca_i"))


                                                       
{
    G_Na = 4.0; //activated sodium conductance
    G_NaC = 0.003; // background sodium conductance
    E_Na = 50.0; // (mV) Sodium Nernst equilibrium potential
    G_s = 0.009; // slow inward calcium current conductance
    
}


void FHNagumo::initQpStatefulProperties()
{
    _m[_qp] = 0.0;
    _h[_qp] = 0.0;
    _j[_qp] = 0.0;
    _d[_qp] = 0.0;
    _f[_qp] = 0.0;
    _x1[_qp] = 0.0;
    _Ca_i[_qp] = 0.0004;
    
}

void FHNagumo::computeQpProperties()
{
    //initialization
    
    V = _V_old[_qp];
    
    explicit = true; //= true Explicit Euler ... =false Implicit Euler
    
    
    //update Calcium dynamics

    _Ca_i[_qp] = update_Ca(Real _Ca_i_old[_qp], Real _dt);
    
    //update Calcium Equilibrium Potential
    
    E_s = -82.3 - 13.0287 * std::log(_Ca_i[_qp]);

    //update Gating
    
    _m[_qp] = update_m(V, _m_old[_qp], _dt);
    _h[_qp] = update_h(Real V, _h_old[_qp], _dt);
    _j[_qp] = update_j(Real V, Real _j_old[_qp], Real _dt);
    _d[_qp] = update_d(Real V, Real _d_old[_qp], Real _dt);
    _f[_qp] = update_f(Real V, Real _f_old[_qp], Real _dt);
    _x1[_qp] = update_x1(Real V, Real _x1_old[_qp], Real _dt);
    
    // Compute Ionic Functions
    
    i_K1 = 0.35*((4*(std::exp(0.04*(V+85.0))-1.0))/(std::exp(0.08*(V+53.0)) + std::exp(0.04*(V+53.0)))+(0.2*(V+23.0))/(1.0 - std::exp(-0.04*(V+23.0))));
    I_X1 = 0.8*(std::exp(0.04*(V+77.0))-1.0)/(std::exp(0.04*(V+35.0)));
    i_x1 = I_X1 * _x1[_qp];
    i_Na = (G_Na * pow(_m[_qp], 3) * _h[_qp] * _j[_qp] + G_NaC) * (V - E_Na);
    i_s = G_s * _d[_qp] * _f[_qp] * (V - E_s);
    
    _I_ion[_qp] = i_K1 + i_x1 + i_Na + i_s;
    

}

inline Real BeelerReuter::update_m(Real V, Real m_old, Real dt){
    
    alpha_m = ((0.0*std::exp(0.0*(V+47.0))-1.0*(V+47.0))/(std::exp(-0.1*(V+47.0))-1.0));
    beta_m = ((40.0*std::exp(-0.056*(V+72.0))+0.0*(V+0.0))/(std::exp(0.0*(V+72.0))+0.0));
    
    if (explicit){
    //Explicit Euler
           m = (1.0 - dt * (alpha_m+beta_m)) * m_old + dt * alpha_m;
    }
    
    else{
        //Implicit Euler
        m = (m_old + dt * alpha_m)/(1.0 + dt * (alpha_m+beta_m));
    }
    return m;
    
}

inline Real BeelerReuter::update_h(Real V, Real h_old, Real dt){
    
    alpha_h = ((0.126*std::exp(-0.25*(V+77.0))+0.0*(V+0.0))/(std::exp(0.0*(V+77.0))+0.0));
    beta_h = ((1.7*std::exp(0.0*(V+22.5))+0.0*(V+0.0))/(std::exp(-0.082*(V+22.5))+1.0));
    
    if (explicit){
    //Explicit Euler
           h = (1.0 - dt * (alpha_h+beta_h)) * h_old + dt * alpha_h;
    }
    
    else{
        //Implicit Euler
        h = (h_old + dt * alpha_h)/(1.0 + dt * (alpha_h+beta_h));
    }
    
    return h;
}

inline Real BeelerReuter::update_j(Real V, Real h_old, Real dt){
    
    alpha_j = ((0.055*std::exp(-0.25*(V+78.0))+0.0*(V+0.0))/(std::exp(-0.2*(V+78.0))+1.0));
    beta_j = ((0.3*std::exp(0.0*(V+32.0))+0.0*(V+0.0))/(std::exp(-0.1*(V+32.0))+1.0));
    
    if (explicit){
    //Explicit Euler
           j = (1.0 - dt * (alpha_j+beta_j)) * j_old + dt * alpha_j;
    }
    
    else{
        //Implicit Euler
        j = (j_old + dt * alpha_j)/(1.0 + dt * (alpha_j+beta_j));
    }
    
    return j;
    
}

inline Real BeelerReuter::update_d(Real V, Real d_old, Real dt){
    
    alpha_d = ((0.095*std::exp(-0.01*(V-5.0))+0.0*(V+0.0))/(std::exp(-0.072*(V-5.0))+1.0));
    beta_d = ((0.07*std::exp(-0.017*(V+44.0))+0.0*(V+0.0))/(std::exp(0.05*(V+44.0))+1.0));
    
    if (explicit){
    //Explicit Euler
           d = (1.0 - dt * (alpha_d+beta_d)) * d_old + dt * alpha_d;
    }
    
    else{
        //Implicit Euler
        d = (d_old + dt * alpha_d)/(1.0 + dt * (alpha_d+beta_d));
    }
    
    return d;
    
}

inline Real BeelerReuter::update_f(Real V, Real f_old, Real dt){
    
    alpha_f = ((0.012*std::exp(-0.008*(V+28.0))+0.0*(V+0.0))/(std::exp(0.15*(V+28.0))+1.0));
    beta_f = ((0.0065*std::exp(-0.02*(V+30.0))+0.0*(V+0.0))/(std::exp(-0.2*(V+30.0))+1.0));
    
    if (explicit){
    //Explicit Euler
           f = (1.0 - dt * (alpha_f+beta_f)) * f_old + dt * alpha_f;
    }
    
    else{
        //Implicit Euler
        f = (f_old + dt * alpha_f)/(1.0 + dt * (alpha_f+beta_f));
    }
    
    return f;
    
}

inline Real BeelerReuter::update_x1(Real V, Real x1_old, Real dt){
    
    alpha_x1 = ((0.0005*std::exp(0.083*(V+50.0))+0.0*(V+0.0))/(std::exp(0.057*(V+50.0))+1.0));
    beta_x1 = ((0.0013*std::exp(-0.06*(V+20.0))+0.0*(V+0.0))/(std::exp(-0.04*(V+20.0))+1.0));
    
    if (explicit){
    //Explicit Euler
           x1 = (1.0 - dt * (alpha_x1+beta_x1)) * x1_old + dt * alpha_x1;
    }
    
    else{
        //Implicit Euler
        x1 = (x1_old + dt * alpha_x1)/(1.0 + dt * (alpha_x1+beta_x1));
    }
    
    return x1;
    
}

inline Real BeelerReuter::update_Ca(Real V, Real Ca_old, Real dt){
    
    if (explicit){
    //Explicit Euler
        Ca = - pow(10, -7) * (i_s + 0.07) * _dt + (1 - 0.07) * _Ca_i_old[_qp];
    }
    
    else{
        //Implicit Euler
        x1 = (x1_old + dt * alpha_x1)/(1.0 + dt * (alpha_x1+beta_x1));
    }
    
    return Ca;
    
}
