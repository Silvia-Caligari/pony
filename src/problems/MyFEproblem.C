#include "MyFEproblem.h"
#include "NonlinearSystemBase.h"
#include "DisplacedProblem.h"
#include "Assembly.h"

registerMooseObject("ponyApp",MyFEproblem);

defineLegacyParams(MyFEproblem);

InputParameters
MyFEproblem::validParams()
{
  InputParameters params = FEProblem::validParams();
//InputParameters params = SubProblem::validParams();
    params.addRequiredParam<Real>("alpha","splitting coefficient");
    params.addRequiredParam<Real>("theta","theta-method coefficient");
    params.addRequiredParam<Real>("c","FHN coefficient");
    params.addRequiredParam<Real>("u_rest","resting potential");
    params.addRequiredParam<Real>("u_thresh","threshold potential");
    params.addRequiredParam<Real>("u_dep","depolar potential");
    params.addRequiredParam<Real>("beta","w coefficient potential equation");
    params.addRequiredParam<Real>("delta","v coefficient in w equation");
    params.addRequiredParam<Real>("gamma","w coefficient in w equation");
    params.addRequiredParam<bool>("method","Newton or Halley method");
    params.addRequiredParam<bool>("how_solve","how_solve_system");
    params.addRequiredParam<std::string>("material_name", "material_name");

  return params;
}

MyFEproblem::MyFEproblem(const InputParameters & parameters):
FEProblem(parameters),
_system(_eq.get_system<TransientNonlinearImplicitSystem>("nl0")),
_solutionPointer(_system.solution),
_meshSystem(_system.get_mesh()),
_dof_map(_system.get_dof_map()),
_myAssembly(assembly(0)),
_qrule(_myAssembly.qRule()),
_material_name(getParam<std::string>("material_name")),
//_materialBase(*getMaterial(_material_name.c_str(),Moose::BLOCK_MATERIAL_DATA,0,true)),
//_epmaterials(dynamic_cast<EPmaterialsNew &>(_materialBase)),

//M(_system.get_system_matrix()),
//A(_system.get_system_matrix()),
_alpha(getParam<Real>("alpha")),
_theta(getParam<Real>("theta")),
_c(getParam<Real>("c")),
u_rest(getParam<Real>("u_rest")),
u_thresh(getParam<Real>("u_thresh")),
u_dep(getParam<Real>("u_dep")),
beta(getParam<Real>("beta")),
delta(getParam<Real>("delta")),
gamma(getParam<Real>("gamma")),
_method(getParam<bool>("method")),
_howsolve(getParam<bool>("how_solve"))
{}


void MyFEproblem::initialSetup(){

    FEProblem::initialSetup();
    
    NumericVector<Number> * sol=_solutionPointer.get();
    
    _w.assign(sol[0].last_local_index() - sol[0].first_local_index(), 0.0);

    
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



void MyFEproblem::solve()
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
    
 

    //////////////////////////////////////////////////////////////////////////
    solveL1(_alpha, _theta, _c, u_rest, u_thresh, u_dep, _w, beta, delta, gamma);
    
    
    //solveL1(_alpha, _theta, _c, u_rest, u_thresh, u_dep);
    
    
    if (_solve){
    _nl->solve();
    }
    Real alpha_bis = 1.0 - _alpha;
    
    solveL1(alpha_bis, _theta, _c, u_rest, u_thresh, u_dep, _w, beta, delta, gamma);
        //solveL1(alpha_bis, _theta, _c, u_rest, u_thresh, u_dep);
        
//}

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

void MyFEproblem::solveL1(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep, std::vector<Number> &_w, Real beta, Real delta, Real gamma){
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
        
        //std::cout<<"w.at:"<<_w.at(i - sol[0].first_local_index())<<std::endl;
        //std::cout<<"v:"<<v<<std::endl;
        
    Real w_new;
    
    //v += _alpha * _dt * _c * (v - u_rest) * (u_dep - v) * (v - u_thresh);
        //NEWTWON METHOD (Solve Non linear system L_1)
        
        w_new = (1.0/((1.0/_dt)+_alpha * gamma * _theta))*(((1.0/_dt) - _alpha * gamma * (1.0 - _theta))*w + _alpha * delta * v);
        
        int N = 100;
        int k = 0;
        Real tol = 1e-8;
        Real error = 1.0;
        Real v0 = v;
        Real v_new;
        
        while(k<=N & error >= tol){
            
           Real Fv = Function(_alpha, _theta, _c, v0, u_rest, u_thresh, u_dep, v, w_new);
            //Real Fv = Function(_alpha, _theta, _c, v0, u_rest, u_thresh, u_dep, v);
           Real DFv = DFunction(_alpha, _theta, _c, u_rest, u_thresh, u_dep, v);
           
            
            if(std::abs(DFv) == 0)
                break;
            
            else{
                
                if(_method)
                    
                    v_new = v - (Fv/DFv);
                
                else{
                    
                    Real DDFv = DDFunction(_alpha, _theta, _c, u_rest, u_thresh, u_dep, v);
                    
                    v_new = v - ((2*Fv*DFv)/(2*DFv*DFv - Fv*DDFv));
                }
                
                error = std::abs(v_new - v);
                v = v_new;
                k +=1;
            }
            
        }
    
    _w[i - sol[0].first_local_index()] = w_new;
        //std::cout<<"_w:"<<_w[i - sol[0].first_local_index()]<<std::endl;
    sol[0].set(i,v);
    
}
    
    sol[0].close();
    
    _system.update();
    
}

inline Real MyFEproblem::Function(Real _alpha, Real _theta, Real _c, Real v0, Real u_rest, Real u_thresh, Real u_dep, Real v, Real w){
    
    Real F;
    
    F = (1/_dt) * v + _alpha * _theta * _c * (v - u_rest) * (v - u_dep) * (v - u_thresh) - (1/_dt) * v0 + _alpha * (1 - _theta) * _c * (v0 - u_rest) * (v0 - u_dep) * (v0 - u_thresh) + _alpha*beta*w;
    //F = (1/_dt) * v + _alpha * _theta * _c * (v - u_rest) * (v - u_dep) * (v - u_thresh) - (1/_dt) * v0 + _alpha * (1 - _theta) * _c * (v0 - u_rest) * (v0 - u_dep) * (v0 - u_thresh);
    
    return F;
    
}
/*
inline Real MyFEproblem::Function(Real _alpha, Real _theta, Real _c, Real v0, Real u_rest, Real u_thresh, Real u_dep, Real v){
    
    Real F;
    
    //F = (1/_dt) * v + _alpha * _theta * _c * (v - u_rest) * (v - u_dep) * (v - u_thresh) - (1/_dt) * v0 + _alpha * (1 - _theta) * _c * (v0 - u_rest) * (v0 - u_dep) * (v0 - u_thresh) + _alpha*beta*w;
    F = (1/_dt) * v + _alpha * _theta * _c * (v - u_rest) * (v - u_dep) * (v - u_thresh) - (1/_dt) * v0 + _alpha * (1 - _theta) * _c * (v0 - u_rest) * (v0 - u_dep) * (v0 - u_thresh);
    
    return F;
    
}*/


inline Real MyFEproblem::DFunction(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep, Real v){
    
    Real dF;
    
    dF = 1/_dt + _alpha * _theta * _c * (3*v*v - 2*(u_dep + u_rest + u_thresh)*v + u_thresh*u_dep + u_thresh*u_rest + u_rest*u_dep);
    
    return dF;
    
}

inline Real MyFEproblem::DDFunction(Real _alpha, Real _theta, Real _c, Real u_rest, Real u_thresh, Real u_dep, Real v){
    
    Real ddF;
    
    ddF = _alpha * _theta * _c * (6*v - 2*(u_dep + u_rest + u_thresh));
    
    return ddF;
    
}


void MyFEproblem::computeResidual(const NumericVector< Number > & soln, NumericVector<Number> &     residual){
    
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
        
        
        residual.print_matlab("residual1.m");
        //sol[0].print_matlab("sol1.m");
        //soln.print_matlab("soln1.m");
    
        
    }
     
    else{
        
    FEProblemBase::computeResidual(soln, residual);
    
   
    
    residual.print_matlab("residual2.m");
    //sol[0].print_matlab("sol2.m");
    //soln.print_matlab("soln2.m");
        
    }
    
}

void MyFEproblem::computeJacobian(const NumericVector<Number> & soln, SparseMatrix<Number> & jacobian){
    
    if(_howsolve){
        
    }
    
    
    else
        
        FEProblemBase::computeJacobian(soln, jacobian);
    

    /*SparseMatrix<Number> & S = _system.get_system_matrix();
    
    S.print_matlab("ciao2.m");
    
    jacobian.print_matlab("jacobian2.m");*/


    
}

