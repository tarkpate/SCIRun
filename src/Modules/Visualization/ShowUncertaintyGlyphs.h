#ifndef MODULES_FIELDS_SHOW_UNCERTAINTY_GLYPHS_H_
#define MODULES_FIELDS_SHOW_UNCERTAINTY_GLYPHS_H_

#include <Dataflow/Network/Module.h>
#include <Modules/Fields/share.h>

namespace SCIRun {
namespace Modules {
namespace Fields {
  class SCISHARE ShowUncertaintyGlyphs : public SCIRun::Dataflow::Networks::Module,
          public Has1InputPort<FieldPortTag>,
          public Has1OutputPort<FieldPortTag>
  {
  public:
    ShowUncertaintyGlyphs();
    virtual void execute();
    virtual void setStateDefaults();

    INPUT_PORT(0, InputField, Field);
    OUTPUT_PORT(0, OutputField, Field);

    MODULE_TRAITS_AND_INFO(ModuleHasUI);
  };
}}}

#endif
