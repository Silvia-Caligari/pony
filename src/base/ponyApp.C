#include "ponyApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
ponyApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

ponyApp::ponyApp(InputParameters parameters) : MooseApp(parameters)
{
  ponyApp::registerAll(_factory, _action_factory, _syntax);
}

ponyApp::~ponyApp() {}

void
ponyApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"ponyApp"});
  Registry::registerActionsTo(af, {"ponyApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
ponyApp::registerApps()
{
  registerApp(ponyApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ponyApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ponyApp::registerAll(f, af, s);
}
extern "C" void
ponyApp__registerApps()
{
  ponyApp::registerApps();
}
