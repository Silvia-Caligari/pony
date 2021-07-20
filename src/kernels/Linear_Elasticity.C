//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Linear_Elasticity.h"

registerMooseObject("ponyApp",Linear_Elasticity);

defineLegacyParams(Linear_Elasticity);

InputParameters Linear_Elasticity::validParams()
{
  InputParameters params = Kernel::validParams();
  
  params.addRequiredCoupledVar("disp_x", "first displacement component");
  params.addCoupledVar("disp_y", "second displacement component");
  params.addCoupledVar("disp_z", "third displacement component");
    
  return params;
}

Linear_Elasticity::Linear_Elasticity(const InputParameters & parameters) :
Kernel(parameters),
_ind_disp_x(coupled("disp_x")),
_ind_disp_y(coupled("disp_y")),
_ind_disp_z(coupled("disp_z")),
_S(getMaterialProperty<RealTensorValue>("_S"))
{
    
    
    
    if(_var.number() == _ind_disp_x)
        
        _ind = 0;
    
    if(_var.number() == _ind_disp_y)
        
        _ind = 1;
    
    if(_var.number() == _ind_disp_z)
        
        _ind = 2;
        
    
}

    
Real Linear_Elasticity::computeQpResidual()
{
    RealTensorValue T;
    computebasisTensor(T, _grad_test[_i][_qp], _ind);
    
    return _S[_qp].contract(T);
	
}


void Linear_Elasticity::computebasisTensor(RealTensorValue &T, const RealVectorValue &phi, const int ind)
{
    
    T.zero();
    
    for (int j=0; j<LIBMESH_DIM; ++j){
        
        T(ind,j) = phi(j);
        
    }
    
    T = (T + T.transpose())/2.0;
    
}
