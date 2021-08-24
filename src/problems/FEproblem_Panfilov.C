#include "FEproblem_Panfilov.h"
#include "NonlinearSystemBase.h"
#include "DisplacedProblem.h"
#include "Assembly.h"


registerMooseObject("ponyApp",FEproblem_Panfilov);

defineLegacyParams(FEproblem_Panfilov);

InputParameters
FEproblem_Panfilov::validParams()
{
  InputParameters params = FEProblem::validParams();
//InputParameters params = SubProblem::validParams();
    params.addRequiredParam<Real>("alpha","splitting coefficient");
    params.addRequiredParam<Real>("theta","theta-method coefficient");
    params.addRequiredParam<Real>("k","first coefficient");
    params.addRequiredParam<Real>("a","second coefficient");
    params.addRequiredParam<Real>("e0","parameter 1");
    params.addRequiredParam<Real>("mu_1","parameter 2");
    params.addRequiredParam<Real>("mu_2","parameter 3");
    params.addRequiredParam<bool>("method","Newton or Halley method");
    params.addRequiredParam<bool>("how_solve","how_solve_system");
    params.addRequiredParam<std::string>("material_name", "material_name");


  return params;
}

FEproblem_Panfilov::FEproblem_Panfilov(const InputParameters & parameters):
FEProblem(parameters),
_system(_eq.get_system<TransientNonlinearImplicitSystem>("nl0")),
_solutionPointer(_system.solution),
_meshSystem(_system.get_mesh()),
_dof_map(_system.get_dof_map()),
_myAssembly(assembly(0)),
_qrule(_myAssembly.qRule()),
_material_name(getParam<std::string>("material_name")),
_alpha(getParam<Real>("alpha")),
_theta(getParam<Real>("theta")),
k(getParam<Real>("k")),
a(getParam<Real>("a")),
e0(getParam<Real>("e0")),
mu_1(getParam<Real>("mu_1")),
mu_2(getParam<Real>("mu_2")),
_method(getParam<bool>("method")),
_howsolve(getParam<bool>("how_solve"))
{}


void FEproblem_Panfilov::initialSetup(){

    FEProblem::initialSetup();
 
    NumericVector<Number> * sol=_solutionPointer.get();
    
    _w.assign(sol[0].last_local_index() - sol[0].first_local_index(), 0.0);
    
    //std::cout<<"length:"<<sol[0].last_local_index() - sol[0].first_local_index()<<std::endl;
    
    if(_howsolve){
    MaterialBase & _materialBase(*getMaterial(_material_name.c_str(),Moose::BLOCK_MATERIAL_DATA,0,true));
    EPmaterialsNew & _epmaterials(dynamic_cast<EPmaterialsNew &>(_materialBase));
 
    
    
    //std::cout<<"length:"<<sol[0].last_local_index() - sol[0].first_local_index()<<std::endl;
    
    _dim = _meshSystem.mesh_dimension();
    
    SparseMatrix<Number> & M = _system.add_matrix("mass");  //Mass Matrix
    SparseMatrix<Number> & A = _system.add_matrix("stiff");//Stiffness Matrix * Tensor Diffusion
    
    SparseMatrix<Number> & A1 = _system.add_matrix("rhs_linear_system"); //A1 = (M/dt - (1-theta)A)
    SparseMatrix<Number> & A2 = _system.get_system_matrix(); //A2 = (M/dt + theta*A)
    
    //_eq.init();
    
    
    M.zero();
    A.zero();
    A1.zero();
    A2.zero();

    
    FEType fe_type = _dof_map.variable_type(0);
    
    ////FEM information///
     
    std::unique_ptr<FEBase> fe (FEBase::build(_dim, fe_type));
    std::unique_ptr<QBase> qrule( QBase::build (_qrule->type(),_dim,_qrule->get_order()));
    fe->attach_quadrature_rule (qrule.get());
    
    std::vector<Real>                       const & JxW = fe->get_JxW();
    std::vector<std::vector<Real> >         const & phi = fe->get_phi();
    std::vector<std::vector<RealGradient> > const & dphi = fe->get_dphi();
    std::vector<Point>                      const & q_points = fe->get_xyz();
    
    std::vector<dof_id_type> dof_indices;
    
    //local matrix
    
    DenseMatrix<Number> ke;
    DenseMatrix<Number> kem;
    
    //active elements
    
    MeshBase::const_element_iterator           el = _meshSystem.active_local_elements_begin();
    MeshBase::const_element_iterator const end_el = _meshSystem.active_local_elements_end();
    
    
    std::vector<RealTensorValue> diffusion;
    
    //start element loop
    
    for ( ; el != end_el; ++el)
      {
        Elem const * elem = *el;
        
        fe->reinit (elem);
        _dof_map.dof_indices(elem, dof_indices);

        unsigned int const n_dofs = cast_int<unsigned int>(dof_indices.size());

        ke.resize (n_dofs , n_dofs);
        ke.zero();
          
        kem.resize (n_dofs , n_dofs);
        kem.zero();
         
        _epmaterials.getDiffusion(q_points, diffusion);
          
        for (unsigned int i=0; i<phi.size(); i++)
          for (unsigned int j=0; j<phi.size(); j++)
            for (unsigned int qp=0; qp<qrule->n_points(); qp++)
            {
              ke(i,j) +=  JxW[qp] * ( dphi[j][qp] * (diffusion.at(qp) *  dphi[i][qp]) );
              
              kem(i,j) +=  JxW[qp] * ( phi[j][qp] * phi[i][qp] );
            }
          
          A.add_matrix(ke, dof_indices );
          M.add_matrix(kem, dof_indices );
          
          A1.add_matrix(ke, dof_indices );
          A2.add_matrix(ke, dof_indices );
      } //end element loop
    
    A.close();
    M.close();
    A1.close();
    A2.close();
        
    
    }
  
    
}



void FEproblem_Panfilov::solve()
{
    //std::cout<<"SolveMYFe"<<std::endl;

    
  //TIME_SECTION(_solve_timer);

  // This prevents stale dof indices from lingering around and possibly leading to invalid reads and
  // writes. Dof indices may be made stale through operations like mesh adaptivity
  clearAllDofIndices();
  if (_displaced_problem)
    _displaced_problem->clearAllDofIndices();

#ifdef LIBMESH_HAVE_PETSC
#if PETSC_RELEASE_LESS_THAN(3, 12, 0)
  Moose::PetscSupport::petscSetOptions(*this); // Make sure the PETSc options are setup for this app
#else
  // Now this database will be the default
  // Each app should have only one database
  if (!_app.isUltimateMaster())
    PetscOptionsPush(_petsc_option_data_base);
  // We did not add petsc options to database yet
  if (!_is_petsc_options_inserted)
  {
    Moose::PetscSupport::petscSetOptions(*this);
    _is_petsc_options_inserted = true;
  }
#endif
  // set up DM which is required if use a field split preconditioner
  // We need to setup DM every "solve()" because libMesh destroy SNES after solve()
  // Do not worry, DM setup is very cheap
  if (_nl->haveFieldSplitPreconditioner())
    Moose::PetscSupport::petscSetupDM(*_nl);
#endif

  Moose::setSolverDefaults(*this);

  // Setup the output system for printing linear/nonlinear iteration information
  initPetscOutput();

  possiblyRebuildGeomSearchPatches();

  // reset flag so that linear solver does not use
  // the old converged reason "DIVERGED_NANORINF", when
  // we throw  an exception and stop solve
  //_fail_next_linear_convergence_check = false;
  
    //std::cout<<sol[0].size()<<std::endl
    
    

    //////////////////////////////////////////////////////////////////////////
    solveL1(_alpha, _theta, k, a, e0, mu_1, mu_2, _w);
    
    
    
    if (_solve){
    _nl->solve();
    
    Real alpha_bis = 1.0 - _alpha;
    
        solveL1(alpha_bis, _theta, k, a, e0, mu_1, mu_2, _w);
  
        
    }

  if (_solve)
    _nl->update();

  // sync solutions in displaced problem
  if (_displaced_problem)
    _displaced_problem->syncSolutions();

#if !PETSC_RELEASE_LESS_THAN(3, 12, 0)
  if (!_app.isUltimateMaster())
    PetscOptionsPop();
#endif
    
    
}

void FEproblem_Panfilov::solveL1(Real _alpha, Real _theta, Real k, Real a, Real e0, Real mu_1, Real mu_2, std::vector<Number> &_w){
//void MyFEproblem::solveL1(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep){
    
    //EquationSystems & _equationSystems=es();
    
    /*TransientNonlinearImplicitSystem & _system = _equationSystems.get_system<TransientNonlinearImplicitSystem>("nl0");
    //TransientNonlinearImplicitSystem & _system = _eq.get_system<TransientNonlinearImplicitSystem>("nl0");
    
    std::unique_ptr<NumericVector<Number>> & solutionPointer(_system.solution);*/
    
    //NumericVector<Number> * sol=solutionPointer.get();
    
    NumericVector<Number> * sol=_solutionPointer.get();
   
    
    for(int i=sol[0].first_local_index(); i<sol[0].last_local_index(); ++i){
        
    Real v = sol[0](i);
        
    Real w = _w.at(i - sol[0].first_local_index());
        
    v0 = v;
    w0 = w;
        
        if(_method){
    
        //NEWTWON METHOD (Solve Non linear system L_1)
            
            w_new = Newton_Gating(_alpha, _theta, k, a, e0, mu_1, mu_2, v, w, w0, w_new);
            
            v_new = Newton_Potential(_alpha, _theta, k, a, v, w, v0, v_new);

            
        }
        
        else{
            
            w_new = w0 + _alpha * _dt * G2_Function(k, a, e0, mu_1, mu_2, v0, w0 );
            
            v_new = v0 +_alpha * _dt * (G1_Function(k, a, v0, w0));
            
            Real w_s = 0.5*(w0+w_new);
            
            Real v_s = 0.5*(v0+v_new);
            
            w_new = w0 + _alpha * _dt * G2_Function(k, a, e0, mu_1, mu_2, v_s, w_s);
            
            v_new = v0 +_alpha * _dt * (G1_Function(k, a, v_s, w_s));
            
            
        }
    
    _w[i - sol[0].first_local_index()] = w_new;
    
    v = v_new;
        
    sol[0].set(i,v);
    
}
    
    sol[0].close();
    
}

inline Real FEproblem_Panfilov::Function_v(Real _alpha, Real _theta, Real k, Real a, Real v, Real w, Real v0){
    
    Real Fv;
    
    Fv = (1/_dt) * v - _alpha * _theta * (G1_Function(k, a, v, w)) - (1/_dt) * v0 - _alpha * (1 - _theta) * (G1_Function(k, a, v0, w));
    
    return Fv;
    
}

inline Real FEproblem_Panfilov::DFunction_v(Real _alpha, Real _theta, Real k, Real a, Real v, Real w){
    
    Real dF_v;
    
    dF_v = 1/_dt + _alpha * _theta * dG1_Function(k, a, v, w);
    
    return dF_v;
    
}


inline Real FEproblem_Panfilov::Function_w(Real _alpha, Real _theta, Real k, Real a, Real e0, Real mu_1, Real mu_2, Real v, Real w, Real w0){
    
    Real Fw;
    
    Fw = (1/_dt) * w - _alpha * _theta * G2_Function(k, a, e0, mu_1, mu_2, v, w) - (1/_dt) * w0 + _alpha * (1 - _theta) * G2_Function(k, a, e0, mu_1, mu_2, v, w0);
    
    return Fw;
    
}


inline Real FEproblem_Panfilov::DFunction_w(Real _alpha, Real _theta, Real k, Real a, Real e0, Real mu_1, Real mu_2, Real v, Real w){
    
    Real dF_w;
    
    dF_w = 1/_dt - _alpha * _theta * dG2_Function( k, a, e0, mu_1, mu_2, v, w);
    
    return dF_w;
    
}


inline Real FEproblem_Panfilov::G1_Function(Real k, Real a, Real v, Real w){
    
    Real G1;
    
    G1 = -k*v*(v-a)*(v-1.0)-v*w;
    
    return G1;
}

inline Real FEproblem_Panfilov::G2_Function(Real k, Real a, Real e0, Real mu_1, Real mu_2, Real v, Real w){
    
    Real G2;
    
    G2 = (e0+mu_1*w/(v+mu_2))*(-w -k*v*(v-a-1.0));
    
    return G2;
}

inline Real FEproblem_Panfilov::dG1_Function(Real k, Real a, Real v, Real w){
    
    Real dG1;
    
    dG1 = -k*(3*v*v - 2*(a+1)*v +a)-w;
    
    return dG1;
}

inline Real FEproblem_Panfilov::dG2_Function(Real k, Real a, Real e0, Real mu_1, Real mu_2, Real v, Real w){
    
    Real dG2;
    
    dG2 = (e0+mu_1/(v+mu_2))*(-w -k*v*(v-a-1))-(e0+mu_1*w/(v+mu_2));
    
    return dG2;
}

inline Real FEproblem_Panfilov::Newton_Gating(Real _alpha, Real _theta, Real k, Real a, Real e0, Real mu_1, Real mu_2, Real v, Real w, Real w0, Real w_new){
    
    int N = 100;
    int j = 0;
    Real tol = 1e-8;
    Real error = 1.0;
    
    while(j<=N & error >= tol){
        
       Real F = Function_w(_alpha, _theta, k, a, e0, mu_1, mu_2, v, w, w0);
   
       Real DF = DFunction_w(_alpha, _theta, k, a, e0, mu_1, mu_2, v, w);
       
        
        if(std::abs(DF) == 0)
            break;
        
        else{
                
                w_new = w - (F/DF);
            
                error = std::abs(w_new - w);
                w = w_new;
                j +=1;
            
        }
        
    }
    
    return w_new;
    
}

inline Real FEproblem_Panfilov::Newton_Potential(Real _alpha, Real _theta, Real k, Real a, Real v, Real w, Real v0, Real v_new){
    
    int N = 100;
    int j = 0;
    Real tol = 1e-8;
    Real error = 1.0;
    
    while(j<=N & error >= tol){
        
       Real F = Function_v(_alpha, _theta, k, a, v, w, v0);
   
       Real DF = DFunction_v(_alpha, _theta, k, a, v, w);
       
        
        if(std::abs(DF) == 0)
            break;
        
        else{
                
                v_new = v - (F/DF);
            
                error = std::abs(v_new - v);
                v = v_new;
                j +=1;
            
        }
        
    }
    
    return v_new;
    
}

void FEproblem_Panfilov::computeResidual(const NumericVector< Number > & soln, NumericVector<Number> &     residual){
    
    NumericVector<Number> * sol=_solutionPointer.get();
    
    if(_howsolve){
        
        
        
        SparseMatrix<Number> & M = _system.get_matrix("mass");
        
        SparseMatrix<Number> & A = _system.get_matrix("stiff");
        
        
        //SparseMatrix<Number> & RHS = _system.get_matrix("rhs_linear_system"); //A1 = (M/dt - (1-theta)A)
        
        SparseMatrix<Number> & S = _system.get_system_matrix();
        
        //RHS.zero();
        S.zero();
        
        //RHS.add(1.0/_dt,M);
        //RHS.add(-(1.0 - _theta), A);

        S.add(1.0/_dt,M);
        S.add(_theta, A);
        
        //RHS.close();
        S.close();
        
        A.vector_mult(residual,sol[0]);
        
        /*sol[0] *= -1.0;
        
        sol[0].close();
        
        S.vector_mult_add(residual,sol[0]);
        
        sol[0] *= -1.0;
         */
        
        //sol[0].close();
        
        residual.close();
        
        
        //residual.print_matlab("residual1.m");
        //sol[0].print_matlab("sol1.m");
        //soln.print_matlab("soln1.m");
    
        
    }
     
    else{
        
    FEProblemBase::computeResidual(soln, residual);
    
   
    
    //residual.print_matlab("residual2.m");
    //sol[0].print_matlab("sol2.m");
    //soln.print_matlab("soln2.m");
        
    }
    
}

void FEproblem_Panfilov::computeJacobian(const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian){
    
    if(_howsolve){
        
    }
    
    
    else
        
        FEProblemBase::computeJacobian(soln, jacobian);
    

    /*SparseMatrix<Number> & S = _system.get_system_matrix();
    
    S.print_matlab("ciao2.m");
    
    jacobian.print_matlab("jacobian2.m");*/


    
}



