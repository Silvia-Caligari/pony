#include "EPmaterialsNew.h"

registerMooseObject("ponyApp", EPmaterialsNew);

defineLegacyParams(EPmaterialsNew);

InputParameters EPmaterialsNew::validParams()
{
	InputParameters params = Material::validParams();

	
    params.addRequiredParam<RealVectorValue>("a_l","forgot first vector for diffusion_tensor ");
    params.addRequiredParam<RealVectorValue>("a_t","forgot second vector for diffusion_tensor ");
	params.addRequiredParam<RealVectorValue>("sigma_i","forgot intracellular diffusion coefficients");
    //params.addRequiredParam<Real>("sigma_e_l","forgot extracellular diffusion coefficient _l");
    //params.addParam<RealVectorValue>("sigma_e","forgot extracellular diffusion coefficients");
    params.addRequiredParam<Real>("C_m","forgot membrane conductance");
    params.addRequiredParam<Real>("Chi","forgot surface per volume");
	return params;
}

EPmaterialsNew::EPmaterialsNew(const InputParameters &parameters) :
Material(parameters),
a_l(getParam<RealVectorValue>("a_l")),
a_t(getParam<RealVectorValue>("a_t")),
sigma_i(getParam<RealVectorValue>("sigma_i")),
//sigma_e(getParam<RealVectorValue>("sigma_e")),
C_m(getParam<Real>("C_m")),
Chi(getParam<Real>("Chi")),
/*sigma_l(declareProperty<Real>("sigma_l")),
sigma_t(declareProperty<Real>("sigma_t")),
sigma_n(declareProperty<Real>("sigma_n")),
_a_l(getMaterialProperty<RealVectorValue>("_a_l")),
_a_t(getMaterialProperty<RealVectorValue>("_a_t")),
_a_n(getMaterialProperty<RealVectorValue>("_a_n")),*/
_K(declareProperty<RealTensorValue>("_K")),
time_coefficient(declareProperty<Real>("time_coefficient"))
                                                       
{}


RealTensorValue product_transpose(const RealVectorValue &v) {
    
    RealTensorValue _T;
    
    for(int i=0; i < LIBMESH_DIM ; i++ ){
        for (int j=0; j < LIBMESH_DIM ; j++){
            
            _T(i,j) = v(i) * v(j);
            
        }
    }
    
    return _T;
        
}


void
EPmaterialsNew::computeQpProperties()
{
   
    
    EvaluateDiffusion(_q_point[_qp], _K[_qp]);
    
    
    time_coefficient[_qp] = C_m * Chi;
    
    

}

void EPmaterialsNew::getDiffusion(std::vector<Point> const & p , std::vector<RealTensorValue> & diffusion)
{
    diffusion.resize( p.size() );
    for (int i=0; i<p.size(); ++i)
    {
        EvaluateDiffusion(p.at(i),diffusion.at(i)) ;
    }
}


void EPmaterialsNew::EvaluateDiffusion(Point const & p, RealTensorValue & diffusion){
    
    RealVectorValue _a_l, _a_t, _a_n;
    
    Real sigma_l, sigma_t, sigma_n;
    
    _a_l = a_l / a_l.norm();
    
    a_t -=((a_l * a_t)/ a_l.norm()) * a_l;
    
    _a_t = a_t / a_t.norm();
    
    _a_n = {0., 0.};
    
    sigma_l = sigma_i(0);
    
    sigma_t = sigma_i(1);
    
    sigma_n = sigma_i(2);

    diffusion = (sigma_l * product_transpose(_a_l)) + (sigma_t * product_transpose(_a_t)) + (sigma_n * product_transpose(_a_n));
    
}

