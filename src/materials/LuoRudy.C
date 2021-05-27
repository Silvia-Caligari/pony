#include "LuoRudy.h"

registerMooseObject("ponyApp",LuoRudy);

defineLegacyParams(LuoRudy);

InputParameters LuoRudy::validParams()
{
	InputParameters params = Material::validParams();

	
    params.addRequiredCoupledVar("potential", "forgot membrane potential Luo-Rudy.");
    return params;
}

LuoRudy::LuoRudy(const InputParameters &parameters) :
Material(parameters),
_V_old(coupledValueOld("potential")),
_I_ion(declareProperty<Real>("I_ion")),
_I_stim(declareProperty<Real>("I_stim")),
_m(declareProperty<Real>("m")),
_h(declareProperty<Real>("h")),
_j(declareProperty<Real>("j")),
_d(declareProperty<Real>("d")),
_f(declareProperty<Real>("f")),
_x(declareProperty<Real>("x")),
_Ca_i(declareProperty<Real>("Ca_i")),
_m_old(getMaterialPropertyOld<Real>("m")),
_h_old(getMaterialPropertyOld<Real>("h")),
_j_old(getMaterialPropertyOld<Real>("j")),
_d_old(getMaterialPropertyOld<Real>("d")),
_f_old(getMaterialPropertyOld<Real>("f")),
_x_old(getMaterialPropertyOld<Real>("x")),
_Ca_i_old(getMaterialPropertyOld<Real>("Ca_i"))


                                                       
{
    K_o = 5.4;
    K_i = 145.0;
    Na_o = 140.0;
    Na_i = 18.0;

    PRNaK = 0.01833;                    // Na/K Permeability Ratio

    E_Na = 54.4;
    E_K = -77.0;                           //(RT/F)*log((Ko + PRNaK*Nao)/(Ko + PRNaK*Nai)); (B-R)
    Ek1 = E_Na * (std::log(K_o/K_i))/(std::log(Na_o/Na_i));
    EKP = Ek1;
    
    GK = 0.282*(std::sqrt(K_o/5.4));
    GK1 = 0.6047*(std::sqrt(K_o/5.4));
    
    _explicit = false; //= true Explicit Euler ... =false Implicit Euler
    
   
    _m_start = 0.0;
    _h_start = 1.0;
    _j_start = 1.0;
    _d_start = 0.0;
    _f_start = 1.0;
    _x_start = 0.0;
    _Ca_start = 0.0002;
    
    
}


void LuoRudy::initQpStatefulProperties()
{
 
   
    _m[_qp] = _m_start;
    _h[_qp] = _h_start;
    _j[_qp] = _j_start;
    _d[_qp] = _d_start;
    _f[_qp] = _f_start;
    _x[_qp] = _x_start;
    _Ca_i[_qp] = _Ca_start;
    
    
}

void LuoRudy::computeQpProperties()
{
    
    V = _V_old[_qp];
    
    //update Calcium Equilibrium Potential
    
    E_si = 7.7 - 13.0287 * std::log(_Ca_i[_qp]);

    //update Gating
    if (_t_step == 1){
        
        _m[_qp] = update_m(V, _m_start, _dt);
        _h[_qp] = update_h(V, _h_start, _dt);
        _j[_qp] = update_j(V, _j_start, _dt);
        _d[_qp] = update_d(V, _d_start, _dt);
        _f[_qp] = update_f(V, _f_start, _dt);
        _x[_qp] = update_x(V, _x_start, _dt);
    }
    
    else {
    _m[_qp] = update_m(V, _m_old[_qp], _dt);
    _h[_qp] = update_h(V, _h_old[_qp], _dt);
    _j[_qp] = update_j(V, _j_old[_qp], _dt);
    _d[_qp] = update_d(V, _d_old[_qp], _dt);
    _f[_qp] = update_f(V, _f_old[_qp], _dt);
    _x[_qp] = update_x(V, _x_old[_qp], _dt);
    }
    
        
    //Computer other Functions
    
    alpha_k1 = 1.02/(1.0 + std::exp(0.2385 * (V - Ek1 - 59.215)));
    
    beta_k1 = (0.49124 * std::exp(0.08032 * (V - Ek1 + 5.476)) + std::exp(0.06175 * (V - Ek1 - 594.31)))/(1.0 + std::exp(-0.5143 * (V - Ek1 + 4.753)));
    
    k1_inf = alpha_k1/(alpha_k1 + beta_k1);
    
    Xi = (2.837 * (std::exp(0.04 * (V + 77.0)) - 1.0)/((V + 77.0) * std::exp(0.04 * (V + 35.0)))) * (V > -100.0) + 1.0 *(V <= -100.0);

    
    KP = 1.0/(1.0 + std::exp((7.488 - V)/5.98));
    
    
    // Compute Ionic Functions
    
    i_Na = 23.0 * pow(_m[_qp],3) * _h[_qp] * _j[_qp] * (V - E_Na);

    
    i_si = 0.09 * _d[_qp] * _f[_qp] * (V - E_si);
        
    i_K = GK * _x[_qp] * Xi * (V - E_K);

    
    i_K1 = GK1 * k1_inf * (V - Ek1);

    
    i_KP = 0.0183 * KP * (V - EKP);

    
    i_b = 0.03921 * (V + 59.87);
    
    
    
    //Compute Stimulus current
    Real x_point = _q_point[_qp](0);
    Real y_point = _q_point[_qp](1);
    
    if(((x_point<0.04)*(y_point<0.04)) && (_t >= 0.0) && (_t <= 1.0)){
    //if((_t >= 0.0) && (_t <= 1.0)){
        _I_stim[_qp] = -100.0;
    }
    
    else
        
        _I_stim[_qp] = 0.0;
     
    
    
    //Compute Ionic current

    _I_ion[_qp] = i_Na + i_si + i_K + i_K1 + i_KP + i_b + _I_stim[_qp];
    
    
    //update Calcium dynamics
    
    if (_t_step ==1){
        
        _Ca_i[_qp] = update_Ca(_Ca_start, _dt, i_si);
        
    }
    
    else{

    _Ca_i[_qp] = update_Ca(_Ca_i_old[_qp], _dt, i_si);
    
    }


}


inline Real LuoRudy::update_m(Real V, Real m_old, Real dt){
    
    Real m;
    
    alpha_m = 0.32 * (V + 47.13)/(1.0 - std::exp(-0.1 * (V + 47.13)));
    beta_m = 0.08 * std::exp(-V/11.0);
    
    
    Real m_inf = alpha_m/(alpha_m + beta_m);
    Real tao_m = 1/(alpha_m + beta_m);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////RUSH LARSEN
    m = m_inf + (m_old - m_inf) * std::exp(-dt/tao_m);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
    /*
    if (_explicit){
    //Explicit Euler
      Real pm = 1.0 - dt*(alpha_m + beta_m);
      m = pm * m_old + dt * alpha_m;
    }
    
    else{
        //Implicit Euler
     Real pm = 1.0 + dt*(alpha_m + beta_m);
     Real qm = m_old + dt * alpha_m;
     m = qm/pm;
    }
    */
    
    return m;
    
}

inline Real LuoRudy::update_h(Real V, Real h_old, Real dt){
    
    Real h;
    
    alpha_h = 0.0 * (V >= -40.0) + 0.135 * (std::exp((80.0 + V)/(-6.8))) * (V < -40.0);
    beta_h = (1.0/(0.13 * (1.0 + std::exp((V + 10.66)/(-11.1))))) * (V >= -40.0) + (3.56 * std::exp(0.079 * V) + 3.1*pow(10.0,5) * std::exp(0.35 * V)) * (V < -40.0);
    
    
    Real h_inf = alpha_h/(alpha_h + beta_h);
    Real tao_h = 1/(alpha_h + beta_h);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //RUSH LARSEN
    h = h_inf + (h_old - h_inf) * std::exp(-dt/tao_h);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /*
    if (_explicit){
    //Explicit Euler
      Real ph = 1.0 - dt*(alpha_h + beta_h);
      h = ph * h_old + dt * alpha_h;
    }
    
    else{
        //Implicit Euler
     Real ph = 1.0 + dt*(alpha_h + beta_h);
     Real qh = h_old + dt * alpha_h;
     h = qh/ph;
    }
    */
    
    return h;
    

}

inline Real LuoRudy::update_j(Real V, Real j_old, Real dt){
    
    
    Real j;
    
    alpha_j = 0.0 * (V >= -40.0) + ((-1.2714e5 * std::exp(0.2444 * V) - 3.474e-5 * std::exp(-0.04391 * V)) * (V + 37.78)/(1 + std::exp(0.311 * (V + 79.23)))) * (V < -40);
    beta_j = (0.3 * std::exp(-2.535e-7 * V)/(1.0 + std::exp(-0.1 * (V + 32.0)))) * (V >= -40.0) + (0.1212 * std::exp(-0.01052 * V)/(1.0 + std::exp(-0.1378 * (V + 40.14)))) * (V < -40.0);
    
    
    Real j_inf = alpha_j/(alpha_j + beta_j);
    Real tao_j = 1/(alpha_j + beta_j);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //RUSH LARSEN
    j = j_inf + (j_old - j_inf) * std::exp(-dt/tao_j);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /*
    if (_explicit){
    //Explicit Euler
      Real pj = 1.0 - dt*(alpha_j + beta_j);
      j = pj * j_old + dt * alpha_j;
    }
    
    else{
        //Implicit Euler
     Real pj = 1.0 + dt*(alpha_j + beta_j);
     Real qj = j_old + dt * alpha_j;
     j = qj/pj;
    }
     
    */
    return j;
    
}

inline Real LuoRudy::update_d(Real V, Real d_old, Real dt){
    
    
    Real d;
    
    alpha_d = 0.095 * std::exp(-0.01 * (V - 5.0))/(1.0 + std::exp(-0.072 * (V - 5)));
    beta_d = 0.07 * std::exp(-0.017 * (V + 44.0))/(1.0 + std::exp(0.05 * (V + 44)));

    
    Real d_inf = alpha_d/(alpha_d + beta_d);
    Real tao_d = 1/(alpha_d + beta_d);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //RUSH LARSEN
    d = d_inf + (d_old - d_inf) * std::exp(-dt/tao_d);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /*
    if (_explicit){
    //Explicit Euler
      Real pd = 1.0 - dt*(alpha_d + beta_d);
      d = pd * d_old + dt * alpha_d;
    }
    
    else{
        //Implicit Euler
     Real pd = 1.0 + dt*(alpha_d + beta_d);
     Real qd = d_old + dt * alpha_d;
     d = qd/pd;
    }
    */
    
    return d;
    
}

inline Real LuoRudy::update_f(Real V, Real f_old, Real dt){
    
    Real f;
    
    alpha_f = 0.012 * std::exp(-0.008 * (V + 28.0))/(1.0 + std::exp(0.15 * (V + 28.0)));

    beta_f = 0.0065 * std::exp(-0.02 * (V + 30.0))/(1.0 + std::exp(-0.2 * (V + 30.0)));

    
    Real f_inf = alpha_f/(alpha_f + beta_f);
    Real tao_f = 1/(alpha_f + beta_f);
    

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //RUSH LARSEN
    f = f_inf + (f_old - f_inf) * std::exp(-dt/tao_f);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /*
    if (_explicit){
    //Explicit Euler
        Real pf = 1.0 - dt*(alpha_f + beta_f);
        
      f = pf * f_old + dt * alpha_f;
  
    }
    
    else{
        //Implicit Euler
     Real pf = 1.0 + dt*(alpha_f + beta_f);
     Real qf = f_old + dt * alpha_f;
     f = qf/pf;
    }
    */
    
    return f;
    
}

inline Real LuoRudy::update_x(Real V, Real x_old, Real dt){
    
    Real x;
    
    alpha_x = 0.0005 * std::exp(0.083 * (V + 50.0))/(1.0 + std::exp(0.057 * (V + 50)));
    beta_x = 0.0013 * std::exp(-0.06 * (V + 20.0))/(1.0 + std::exp(-0.04 * (V + 20.0)));

    
    Real x_inf = alpha_x/(alpha_x + beta_x);
    Real tao_x = 1/(alpha_x + beta_x);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //RUSH LARSEN
    x = x_inf + (x_old - x_inf) * std::exp(-dt/tao_x);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /*
    if (_explicit){
    //Explicit Euler
      Real px = 1.0 - dt*(alpha_x + beta_x);
      x = px * x_old + dt * alpha_x;
    }
    
    else{
        //Implicit Euler
     Real px = 1.0 + dt*(alpha_x + beta_x);
     Real qx = x_old + dt * alpha_x;
     x = qx/px;
    }
     */
    
    return x;
    
}

inline Real LuoRudy::update_Ca(Real Ca_old, Real dt, Real I_s){
    
    Real Ca;
    
    Real pCai = (1.0/dt) + 0.07;
    
    Real qCai = (1.0/dt) * Ca_old + 0.07e-4 - I_s * 1e-4;
    
    //std::cout<<"qCa_i:"<<std::endl;
    //std::cout<<qCai<<std::endl;
    
    Ca = qCai/pCai;

    return Ca;
    
}
