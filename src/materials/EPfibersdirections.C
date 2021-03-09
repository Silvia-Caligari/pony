#include "EPfibersdirections.h"


registerMooseObject("ponyApp", EPfibersdirections);

defineLegacyParams(EPfibersdirections);

InputParameters EPfibersdirections::validParams()
{
    InputParameters params = Material::validParams();

    
    
    params.addRequiredParam<RealVectorValue>("a_l","forgot first vector for diffusion_tensor ");
    params.addRequiredParam<RealVectorValue>("a_t","forgot second vector for diffusion_tensor "); //se secondo param obbligatorio, usa params.addRequiredParam altrimenti solo params.addParam
/*params.addParam<RealVectorValue>("a_n");
    params.addParam<RealVectorValue>("_a_l");
    params.addParam<RealVectorValue>("_a_t");
    params.addParam<RealVectorValue>("_a_n");*/
    return params;
}

EPfibersdirections::EPfibersdirections(const InputParameters &parameters) :
Material(parameters),
a_l(getParam<RealVectorValue>("a_l")),
a_t(getParam<RealVectorValue>("a_t")),
//a_n(getParam<RealVectorValue>("a_n")),
_a_l(declareProperty<RealVectorValue>("_a_l")),
_a_t(declareProperty<RealVectorValue>("_a_t")),
_a_n(declareProperty<RealVectorValue>("_a_n"))
                                                       
{}

void
EPfibersdirections::computeQpProperties()
{
    _a_l[_qp] = a_l / a_l.norm();
    
    a_t -=((a_l * a_t)/ a_l.norm()) * a_l;
    _a_t[_qp] = a_t / a_t.norm();
    
    _a_n[_qp] = {0., 0.};
    
}

