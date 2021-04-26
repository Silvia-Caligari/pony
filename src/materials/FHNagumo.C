#include "FHNagumo.h"

registerMooseObject("ponyApp",FHNagumo);

defineLegacyParams(FHNagumo);

InputParameters FHNagumo::validParams()
{
	InputParameters params = Material::validParams();

	
    params.addRequiredCoupledVar("potential", "forgot membrane potential FHN.");
    params.addRequiredParam<Real>("V_rest","forgot resting potential");
    params.addRequiredParam<Real>("V_thresh","forgot threshold potential");
    params.addRequiredParam<Real>("V_depol","forgot maximum potential");
    params.addRequiredParam<Real>("alpha","alpha*(V-V_rest)*(V-Vdepol)*(V-Vthresh) + beta*w");
    params.addRequiredParam<Real>("beta","alpha*(V-V_rest)*(V-Vdepol)*(V-Vthresh) + beta*w");
    params.addRequiredParam<Real>("delta","dw/dt = delta*V - gamma*w");
    params.addRequiredParam<Real>("gamma","dw/dt = delta*V - gamma*w");
    return params;
}

FHNagumo::FHNagumo(const InputParameters &parameters) :
Material(parameters),
_V_old(coupledValueOld("potential")),
V_rest(getParam<Real>("V_rest")),
V_thresh(getParam<Real>("V_thresh")),
V_depol(getParam<Real>("V_depol")),
alpha(getParam<Real>("alpha")),
beta(getParam<Real>("beta")),
delta(getParam<Real>("delta")),
gamma(getParam<Real>("gamma")),
_I_ion(declareProperty<Real>("I_ion")),
_w(declareProperty<Real>("w")),
_w_old(getMaterialPropertyOld<Real>("w"))

                                                       
{}


void FHNagumo::initQpStatefulProperties()
{
    _w[_qp] = 0.0;
}

void FHNagumo::computeQpProperties()
{
    //initialization
    V = _V_old[_qp];
    
    //update w
    
    _w[_qp] = (1.0 - _dt * gamma)* _w_old[_qp] + _dt* delta * V;
    
    w = _w[_qp];
    
    // Compute Ionic Function
    
    _I_ion[_qp] = alpha * (V - V_depol) * (V - V_thresh) * (V - V_rest) + beta * w;
    
   

}



