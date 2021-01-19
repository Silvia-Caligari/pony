//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ponyTestApp.h"
#include "ponyApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
ponyTestApp::validParams()
{
  InputParameters params = ponyApp::validParams();
  return params;
}

ponyTestApp::ponyTestApp(InputParameters parameters) : MooseApp(parameters)
{
  ponyTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

ponyTestApp::~ponyTestApp() {}

void
ponyTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  ponyApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"ponyTestApp"});
    Registry::registerActionsTo(af, {"ponyTestApp"});
  }
}

void
ponyTestApp::registerApps()
{
  registerApp(ponyApp);
  registerApp(ponyTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ponyTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ponyTestApp::registerAll(f, af, s);
}
extern "C" void
ponyTestApp__registerApps()
{
  ponyTestApp::registerApps();
}
