//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TimeActivation.h"

registerMooseObject("ponyApp", TimeActivation);

defineLegacyParams(TimeActivation);

InputParameters
TimeActivation::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Time Activation");
  params.addRequiredCoupledVar("coupled_variable", "coupled variable");
  params.addParam<Real>("uthresh", 1.0, "uthresh");
  return params;
}

TimeActivation::TimeActivation(const InputParameters & parameters)
  : AuxKernel(parameters),
_coupled_var(coupledValue("coupled_variable")),
_coupled_var_old(coupledValueOld("coupled_variable")),
_u_thresh(getParam<Real>("uthresh"))
{
}

Real
TimeActivation::computeValue()
{
    
    if(_coupled_var_old[_qp] > 0.0){
        if(_coupled_var[_qp] > _u_thresh && _coupled_var_old[_qp] < _u_thresh){
        
        return _coupled_var[_qp]*0 + _t;
            
        }
        
        else{
            return _u_old[_qp];
        }
    }
    
    else return 0;
        
}
