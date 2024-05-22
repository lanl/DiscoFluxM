#include "DiscoFluxMApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
DiscoFluxMApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

DiscoFluxMApp::DiscoFluxMApp(InputParameters parameters) : MooseApp(parameters)
{
  DiscoFluxMApp::registerAll(_factory, _action_factory, _syntax);
}

DiscoFluxMApp::~DiscoFluxMApp() {}

void 
DiscoFluxMApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<DiscoFluxMApp>(f, af, s);
  Registry::registerObjectsTo(f, {"DiscoFluxMApp"});
  Registry::registerActionsTo(af, {"DiscoFluxMApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
DiscoFluxMApp::registerApps()
{
  registerApp(DiscoFluxMApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
DiscoFluxMApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  DiscoFluxMApp::registerAll(f, af, s);
}
extern "C" void
DiscoFluxMApp__registerApps()
{
  DiscoFluxMApp::registerApps();
}
