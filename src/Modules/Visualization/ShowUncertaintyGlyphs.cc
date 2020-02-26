#include <Modules/Visualization/ShowUncertaintyGlyphs.h>
#include <Core/Datatypes/Legacy/Field/Field.h>

using namespace SCIRun::Modules::Fields;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Dataflow::Networks;

MODULE_INFO_DEF(ShowUncertaintyGlyphs, NewField, SCIRun);


ShowUncertaintyGlyphs::ShowUncertaintyGlyphs() : Module(staticInfo_)
{
  INITIALIZE_PORT(InputField);
  INITIALIZE_PORT(OutputField);
}

void ShowUncertaintyGlyphs::setStateDefaults()
{
  auto state = get_state();
}

void ShowUncertaintyGlyphs::execute()
{
  auto field = getRequiredInput(InputField);
}
