#include <Modules/Visualization/ShowUncertaintyGlyphs.h>
#include <Core/Algorithms/Visualization/ShowUncertaintyGlyphsAlgorithm.h>
#include <Core/Datatypes/Legacy/Field/Field.h>

using namespace SCIRun;
using namespace Core::Datatypes;
using namespace Modules::Visualization;
using namespace Core::Algorithms::Visualization;

MODULE_INFO_DEF(ShowUncertaintyGlyphs, Visualization, SCIRun);

namespace SCIRun {
namespace Modules {
namespace Visualization {
ShowUncertaintyGlyphs::ShowUncertaintyGlyphs() : GeometryGeneratingModule(staticInfo_)
{
  INITIALIZE_PORT(InputFields);
  INITIALIZE_PORT(MeanTensorField);
  INITIALIZE_PORT(OutputGeom);
}

void ShowUncertaintyGlyphs::setStateDefaults()
{
  auto state = get_state();
}

void ShowUncertaintyGlyphs::execute()
{
  auto fields = getRequiredDynamicInputs(InputFields);

  if(needToExecute())
  {
    auto algorithm = ShowUncertaintyGlyphsAlgorithm();
    auto output = algorithm.run(*this, withInputData((InputFields, fields)));
    sendOutputFromAlgorithm(MeanTensorField, output);
    sendOutputFromAlgorithm(OutputGeom, output);
  }
}

}}}
