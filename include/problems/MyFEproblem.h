
#pragma once

#include "FEProblem.h"
#include "libmesh/transient_system.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/mesh.h"
#include "libmesh/distributed_mesh.h"

#include "EPmaterialsNew.h"

class MyFEproblem;

template <>
InputParameters validParams<MyFEproblem>();

/**
 * FEProblemBase derived class to enable convergence checking relative to a user-specified
 * postprocessor
 */
class MyFEproblem : public FEProblem
{
public:
  static InputParameters validParams();
    
  MyFEproblem(const InputParameters & parameters);
    
    virtual void initialSetup() override;
    
    virtual void solve() override;
    
    virtual void computeResidual(const NumericVector< Number > & soln, NumericVector< Number > &     residual) override;
    
    virtual void computeJacobian(const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian) override;
    


protected:
    
    TransientNonlinearImplicitSystem & _system;
    std::unique_ptr<NumericVector<Number>> & _solutionPointer;
    MeshBase & _meshSystem;
    DofMap & _dof_map;
    
    Assembly & _myAssembly;
    QBase const * const & _qrule;
    
   const std::string &_material_name;
    
    //MaterialBase &_materialBase;
    
    //EPmaterialsNew &_epmaterials;
    
    Real _alpha;
    Real _theta;
    Real _c;
    Real u_rest;
    Real u_thresh;
    Real u_dep;
    Real beta;
    Real delta;
    Real gamma;
    
    bool _method;
    
    bool _howsolve;
    
    
    
    std::vector<Number> _w;
    
    unsigned int _dim;
    
    
        
    virtual void solveL1(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep, std::vector<Number> &_w, Real beta, Real delta, Real gamma);
    
    inline Real Function(Real _alpha, Real _theta, Real _c, Real v0, Real u_rest, Real u_thresh, Real u_dep, Real v, Real w);
    /*
    virtual void solveL1(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep);
    
    inline Real Function(Real _alpha, Real _theta, Real _c, Real v0, Real u_rest, Real u_thresh, Real u_dep, Real v);*/
    
    inline Real DFunction(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep, Real v);
    
    inline Real DDFunction(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep, Real v);

};

