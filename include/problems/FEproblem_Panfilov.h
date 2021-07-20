
#pragma once

#include "FEProblem.h"
#include "libmesh/transient_system.h"

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
    


protected:
    
    TransientNonlinearImplicitSystem & _system;
    std::unique_ptr<NumericVector<Number>> & _solutionPointer;
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
    
    std::vector<Number> _w;
        
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

