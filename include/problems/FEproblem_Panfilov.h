
#pragma once

#include "FEProblem.h"
#include "libmesh/transient_system.h"

#include "libmesh/sparse_matrix.h"
#include "libmesh/mesh.h"
#include "libmesh/distributed_mesh.h"

#include "EPmaterialsNew.h"

class FEproblem_Panfilov;

template <>
InputParameters validParams<FEproblem_Panfilov>();

/**
 * FEProblemBase derived class to enable convergence checking relative to a user-specified
 * postprocessor
 */
class FEproblem_Panfilov : public FEProblem
{
public:
  static InputParameters validParams();
    
  FEproblem_Panfilov(const InputParameters & parameters);
    
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
    Real _alpha;
    Real _theta;
    Real k;
    Real a;
    Real e0;
    Real mu_1;
    Real mu_2;
    
    Real v0;
    Real w0;
    
    Real v_new;
    Real w_new;


    
    bool _method;
    
    bool _howsolve;
    
    std::vector<Number> _w;
    
    unsigned int _dim;
        
    virtual void solveL1(Real _alpha, Real _theta, Real k, Real a, Real e0, Real mu_1, Real mu_2, std::vector<Number> &_w);
    
    inline Real Function_v(Real _alpha, Real _theta, Real k, Real a, Real v, Real w, Real v0);
    
    inline Real DFunction_v(Real _alpha, Real _theta, Real k, Real a, Real v, Real w);
    
    inline Real Function_w(Real _alpha, Real _theta, Real k, Real a, Real e0, Real mu_1, Real mu_2, Real v, Real w, Real w0);
    
    inline Real DFunction_w(Real _alpha, Real _theta, Real k, Real a, Real e0, Real mu_1, Real mu_2, Real v, Real w);

    inline Real G1_Function(Real k, Real a, Real v, Real w);
    
    inline Real G2_Function(Real k, Real a, Real e0, Real mu_1, Real mu_2, Real v, Real w);
    
    inline Real dG1_Function(Real k, Real a, Real v, Real w);
    
    inline Real dG2_Function(Real k, Real a, Real e0, Real mu_1, Real mu_2, Real v, Real w);
    
    inline Real Newton_Gating(Real _alpha, Real _theta, Real k, Real a, Real e0, Real mu_1, Real mu_2, Real v, Real w, Real w0, Real w_new);
    
    inline Real Newton_Potential(Real _alpha, Real _theta, Real k, Real a, Real v, Real w, Real v0, Real v_new);
    
    

};

