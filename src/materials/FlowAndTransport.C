#include "FlowAndTransport.h"


registerMooseObject("ponyApp", FlowAndTransport);

defineLegacyParams(FlowAndTransport);

InputParameters FlowAndTransport::validParams()
{
	InputParameters params = Material::validParams();

	
	
	params.addRequiredParam<Real>("diffusion1","forgot diffusion1");
	params.addRequiredParam<Real>("diffusion2","forgot diffusion2"); //se secondo param obbligatorio, usa params.addRequiredParam altrimenti solo params.addParam

	return params;
}

FlowAndTransport::FlowAndTransport(const InputParameters &parameters) :
Material(parameters),
_diffusion_1(getParam<Real>("diffusion1")),
_diffusion_2(getParam<Real>("diffusion2")),
_diffusion(declareProperty<Real>("diffusionProperty"))
                                                       
{}

void
FlowAndTransport::computeQpProperties()
{
	_diffusion[_qp]=_diffusion_1;
    
    Real x_point = _q_point[_qp](0);

    Real y_point = _q_point[_qp](1);

    if (x_point < 0.5 && y_point < 0.5)

            _diffusion[_qp] = _diffusion_2;

}
