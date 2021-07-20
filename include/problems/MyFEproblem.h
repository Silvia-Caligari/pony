
#pragma once

#include "FEProblem.h"
#include "libmesh/transient_system.h"

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
    


protected:
    
    TransientNonlinearImplicitSystem & _system;
    std::unique_ptr<NumericVector<Number>> & _solutionPointer;
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
    
    std::vector<Number> _w;
        
    virtual void solveL1(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep, std::vector<Number> &_w, Real beta, Real delta, Real gamma);
    
    inline Real Function(Real _alpha, Real _theta, Real _c, Real v0, Real u_rest, Real u_thresh, Real u_dep, Real v, Real w);
    /*
    virtual void solveL1(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep);
    
    inline Real Function(Real _alpha, Real _theta, Real _c, Real v0, Real u_rest, Real u_thresh, Real u_dep, Real v);*/
    
    inline Real DFunction(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep, Real v);
    
    inline Real DDFunction(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep, Real v);

};

