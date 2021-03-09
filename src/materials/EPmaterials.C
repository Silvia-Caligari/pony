#include "EPmaterials.h"

registerMooseObject("ponyApp", EPmaterials);

defineLegacyParams(EPmaterials);

InputParameters EPmaterials::validParams()
{
	InputParameters params = Material::validParams();

	
	
	params.addRequiredParam<RealVectorValue>("sigma_i","forgot intracellular diffusion coefficients");
    //params.addRequiredParam<Real>("sigma_e_l","forgot extracellular diffusion coefficient _l");
    //params.addParam<RealVectorValue>("sigma_e","forgot extracellular diffusion coefficients");
    params.addRequiredParam<Real>("C_m","forgot membrane conductance");
    params.addRequiredParam<Real>("Chi","forgot surface per volume");
	return params;
}

EPmaterials::EPmaterials(const InputParameters &parameters) :
Material(parameters),
sigma_i(getParam<RealVectorValue>("sigma_i")),
//sigma_e(getParam<RealVectorValue>("sigma_e")),
C_m(getParam<Real>("C_m")),
Chi(getParam<Real>("Chi")),
sigma_l(declareProperty<Real>("sigma_l")),
sigma_t(declareProperty<Real>("sigma_t")),
sigma_n(declareProperty<Real>("sigma_n")),
_a_l(getMaterialProperty<RealVectorValue>("_a_l")),
_a_t(getMaterialProperty<RealVectorValue>("_a_t")),
_a_n(getMaterialProperty<RealVectorValue>("_a_n")),
_K(declareProperty<RealTensorValue>("_K")),
time_coefficient(declareProperty<Real>("time_coefficient"))
                                                       
{}


RealTensorValue transpose_product(const RealVectorValue &v) {
    
    RealTensorValue _T;
    
    for(int i=0; i < LIBMESH_DIM ; i++ ){
        for (int j=0; j < LIBMESH_DIM ; j++){
            
            _T(i,j) = v(i) * v(j);
            
        }
    }
    
    return _T;
        
}


void
EPmaterials::computeQpProperties()
{
   
    sigma_l[_qp] = sigma_i(0);
    
    sigma_t[_qp] = sigma_i(1);
    
    sigma_n[_qp] = sigma_i(2);


    //Real lamba = sigma_i(1) / sigma_e(1);
    //sigma_l[_qp]= (sigma_i(1) * lambda)/(1+lambda);
    //sigma_t[_qp]= (sigma_i(2) * lambda)/(1+lambda);
    //sigma_n[_qp]= (sigma_i(3) * lambda)/(1+lambda);
    
    
    _K[_qp] = (sigma_l[_qp] * transpose_product(_a_l[_qp])) + (sigma_t[_qp] * transpose_product(_a_t[_qp])) + (sigma_n[_qp] * transpose_product(_a_n[_qp]));
    
    
    /*Real x_point = _q_point[_qp](0);

    Real y_point = _q_point[_qp](1);

    if (((x_point*x_point)+(y_point*y_point)) < 0.25)

            _diffusion[_qp] = _diffusion_2;*/
    
    time_coefficient[_qp] = C_m * Chi;
    
    

}


