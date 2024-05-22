//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "DiscoFluxMTestApp.h"
#include "DiscoFluxMApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
DiscoFluxMTestApp::validParams()
{
  InputParameters params = DiscoFluxMApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

DiscoFluxMTestApp::DiscoFluxMTestApp(InputParameters parameters) : MooseApp(parameters)
{
  DiscoFluxMTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

DiscoFluxMTestApp::~DiscoFluxMTestApp() {}

void
DiscoFluxMTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  DiscoFluxMApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"DiscoFluxMTestApp"});
    Registry::registerActionsTo(af, {"DiscoFluxMTestApp"});
  }
}

void
DiscoFluxMTestApp::registerApps()
{
  registerApp(DiscoFluxMApp);
  registerApp(DiscoFluxMTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
DiscoFluxMTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  DiscoFluxMTestApp::registerAll(f, af, s);
}
extern "C" void
DiscoFluxMTestApp__registerApps()
{
  DiscoFluxMTestApp::registerApps();
}
