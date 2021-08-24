//* This file is part of the MOOSE framework
  //* https://www.mooseframework.org
  //*
  //* All rights reserved, see COPYRIGHT for full restrictions
  //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
  //*
  //* Licensed under LGPL 2.1, please see LICENSE for details
  //* https://www.gnu.org/licenses/lgpl-2.1.html

  #include "Isotropic_Elasticity_Tensor.h"

  registerMooseObject("ponyApp", Isotropic_Elasticity_Tensor);

defineLegacyParams(Isotropic_Elasticity_Tensor);
InputParameters Isotropic_Elasticity_Tensor::validParams()
  {
      InputParameters params = Material::validParams();
 

    params.addRequiredCoupledVar("disp_x", "x-component displacement");
    params.addCoupledVar("disp_y", "y-component displacement");
    params.addCoupledVar("disp_z", "z-component dispalcement");
      
    params.addRequiredParam<Real>("mu","mu");
    params.addRequiredParam<Real>("lambda","lambda");

    return params;
  }

  Isotropic_Elasticity_Tensor::Isotropic_Elasticity_Tensor(const InputParameters & parameters) :
  Material(parameters),
  //_dim(_mesh.dimension()),
  _grad_disp_x(coupledGradient("disp_x") ),
  _grad_disp_y(coupledGradient("disp_y")),
  _grad_disp_z(coupledGradient("disp_z")),
  mu(getParam<Real>("mu")),
  lambda(getParam<Real>("lambda")),
 _S(declareProperty<RealTensorValue>("_S"))

  {
      //displacement gradient
      for(int i=0; i < LIBMESH_DIM ; i++ ){
          for (int j=0; j < LIBMESH_DIM ; j++){
              
              _U(i,j) = 0.0;
              
          }
      }

      //identity
      for(int i=0; i < LIBMESH_DIM ; i++ ){
          for (int j=0; j < LIBMESH_DIM ; j++){
              
              _id(i,j) = 0.0;
              
              if(i==j){
                  _id(i,i) = 1.0;
              }
              
          }
      }
  }

  void Isotropic_Elasticity_Tensor::computeQpProperties()
  {

     for (int j=0; j<LIBMESH_DIM; ++j)
     {
       _U(0,j)=_grad_disp_x[_qp](j);
       _U(1,j)=_grad_disp_y[_qp](j);
       _U(2,j)=_grad_disp_z[_qp](j);
     }

       // small deformations tensor
      _E = (_U+_U.transpose() - _id)/2.0;
      
     _S[_qp] = 2.0 * mu * _E + lambda * _E.tr() * _id;

  }

