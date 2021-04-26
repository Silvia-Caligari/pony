#include "FHNionicfunction.h"

registerMooseObject("ponyApp",FHNionicfunction);

defineLegacyParams(FHNionicfunction);

InputParameters
FHNionicfunction::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The nonlinear term f(u) = alpha*u*(u-beta)(gamma-u)-c*w, with the weak "
                             "form of $(\\phi_i, f(u))$.");

  return params;
}

FHNionicfunction::FHNionicfunction(const InputParameters & parameters) :
Kernel(parameters),
_I_ion(getMaterialProperty<Real>("I_ion"))
{}



Real
FHNionicfunction::computeQpResidual()
{
    
        return _I_ion[_qp] * _test[_i][_qp];
   
}

Real
FHNionicfunction::computeQpJacobian()
{
        
        return 0;

    
}
