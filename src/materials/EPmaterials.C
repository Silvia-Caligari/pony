#include "EPmaterials.h"


registerMooseObject("ponyApp", EPmaterials);

defineLegacyParams(EPmaterials);

InputParameters EPmaterials::validParams()
{
	InputParameters params = Material::validParams();

	
	
	params.addRequiredParam<Real>("diffusion_i","forgot intracellular diffusion coefficient");
	params.addParam<Real>("diffusion_e","forgot extracellular diffusion coefficient"); //se secondo param obbligatorio, usa params.addRequiredParam altrimenti solo params.addParam
    params.addRequiredParam<Real>("C_m","forgot membrane conductance");
    params.addRequiredParam<Real>("Chi","forgot surface per volume");
	return params;
}

EPmaterials::EPmaterials(const InputParameters &parameters) :
Material(parameters),
diffusion_i(getParam<Real>("diffusion_i")),
diffusion_e(getParam<Real>("diffusion_e")),
C_m(getParam<Real>("C_m")),
Chi(getParam<Real>("Chi")),
diffusion(declareProperty<Real>("diffusionProperty")),
time_coefficient(declareProperty<Real>("time_coefficient"))
                                                       
{}

void
EPmaterials::computeQpProperties()
{
	
    time_coefficient[_qp] = C_m * Chi;
    
    diffusion[_qp]= (diffusion_i + diffusion_e)/(diffusion_i * diffusion_e);
    
    /*Real x_point = _q_point[_qp](0);

    Real y_point = _q_point[_qp](1);

    if (((x_point*x_point)+(y_point*y_point)) < 0.25)

            _diffusion[_qp] = _diffusion_2;*/

}


