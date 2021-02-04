//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EPtimederivative.h"
#include "Transient.h"
#include "MooseVariableScalar.h"


registerMooseObject("ponyApp",EPtimederivative);

defineLegacyParams(EPtimederivative);

InputParameters
EPtimederivative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The time derivative operator with the weak form of $(\\psi_i,\\frac{\\partial u_h}{\\partial t})$, where partial derivative is evaluated ad \\frac{_u - _u_old}{dt}");
  params.addRequiredParam<int>("order", "need specify time-discretization");

  return params;
}

EPtimederivative::EPtimederivative(const InputParameters & parameters) :
Kernel(parameters),
time_coefficient(getMaterialProperty<Real>("time_coefficient")),
order(getParam<int>("order")),
_u_old(_var.slnOld()),
_u_older(_var.slnOlder())
{}



Real
EPtimederivative::computeQpResidual()
{ 
    int d_order = order;
    
    if(d_order == 1){
    
          return (time_coefficient[_qp]/_dt) * (_u[_qp] - _u_old[_qp]) * _test[_i][_qp];
        }

    else if(d_order == 2){
        
        if(_t_step == 1){
            
            return (time_coefficient[_qp]/_dt) * (_u[_qp] - _u_old[_qp]) * _test[_i][_qp];
        }
        
        else{
        
        Real _u_prime = ( 3.0 *_u[_qp] - 4.0*_u_old[_qp] + _u_older[_qp]) / (2.0 * _dt);
        return time_coefficient[_qp] * _u_prime * _test[_i][_qp];
            
        }
        
    }
    else
    {
        mooseError("time discretization order must be 1 or 2");
        
    }
    
}


Real
EPtimederivative::computeQpJacobian()
{
    int d_order = order;
    if(d_order == 1){
        
     return (time_coefficient[_qp]/_dt) * _phi[_j][_qp] * _test[_i][_qp];
    }
    
    else if(d_order == 2){
        
        if(_t_step == 1){
            
            return (time_coefficient[_qp]/_dt) * _phi[_j][_qp] * _test[_i][_qp];
            
                   }
        else{
            
            Real _u_prime = (3.0 * _phi[_j][_qp])/(2.0 * _dt);
            return  time_coefficient[_qp] * _u_prime * _test[_i][_qp];
            
            }
    }
    
    else
    {
        mooseError("time discretization order must be 1 or 2");
        
    }
  
}

/*Real
EPtimederivative::computeQpResidual()
{
    int d_order = order;
    
    if(d_order == 1){
    
          return (time_coefficient[_qp]/_dt) * (_u[_qp] - _u_old[_qp]) * _test[_i][_qp];
        }

    else
    {
        Real _u_prime = ( 3.0 *_u[_qp] - 4.0*_u_old[_qp] + _u_older[_qp]) / (2.0 * _dt);
        return time_coefficient[_qp] * _u_prime * _test[_i][_qp];
    }
    
}


Real
EPtimederivative::computeQpJacobian()
{
    int d_order = order;
    if(d_order == 1){
        
     return (time_coefficient[_qp]/_dt) * _phi[_j][_qp] * _test[_i][_qp];
    }
    
    else
    {
        Real _u_prime = (3.0 * _phi[_j][_qp])/(2.0 * _dt);
        return  time_coefficient[_qp] * _u_prime * _test[_i][_qp];
        
    }
  
}*/

