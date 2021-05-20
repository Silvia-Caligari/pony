#include "Ionicfunction.h"

registerMooseObject("ponyApp",Ionicfunction);

defineLegacyParams(Ionicfunction);

InputParameters
Ionicfunction::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The nonlinear term f(u) = ionic function "
                             "form of $(\\phi_i, f(u))$.");

  return params;
}

Ionicfunction::Ionicfunction(const InputParameters & parameters) :
Kernel(parameters),
_I_ion(getMaterialProperty<Real>("I_ion"))
{}



Real
Ionicfunction::computeQpResidual()
{
    
        return _I_ion[_qp] * _test[_i][_qp];
    
   
}

Real
Ionicfunction::computeQpJacobian()
{
        
        return 0;

    
}
