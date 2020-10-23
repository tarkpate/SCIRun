#ifndef INTERFACE_MODULES_SHOW_UNCERTAINTY_GLYPHS_H
#define INTERFACE_MODULES_SHOW_UNCERTAINTY_GLYPHS_H

//This file is created from the ShowUncertaintyGlyphsDialog.ui in the module factory.
#include <Interface/Modules/Visualization/ui_ShowUncertaintyGlyphs.h>
#include <boost/shared_ptr.hpp>
#include <Interface/Modules/Base/ModuleDialogGeneric.h>
#include <Interface/Modules/Visualization/share.h>

namespace SCIRun {
namespace Gui {
  class SCISHARE ShowUncertaintyGlyphsDialog : public ModuleDialogGeneric,
                                               public Ui::ShowUncertaintyGlyphsDialog
  {
    Q_OBJECT

  public:
    ShowUncertaintyGlyphsDialog(const std::string& name,
                                SCIRun::Dataflow::Networks::ModuleStateHandle state,
                                QWidget* parent = 0);
  private:
    void setupTensorsTab();
    //this function would be from pulling data from module,
    // usually to change the UI.
    // virtual void pull() override;
  };
}}

#endif
