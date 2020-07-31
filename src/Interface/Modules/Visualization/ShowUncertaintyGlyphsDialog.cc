#include <Interface/Modules/Visualization/ShowUncertaintyGlyphsDialog.h>
#include <Dataflow/Network/ModuleStateInterface.h>
// #include <Core/Algorithms/Visualization/ShowUncertaintyGlyphsAlgo.h>

using namespace SCIRun::Gui;
using namespace SCIRun::Dataflow::Networks;

ShowUncertaintyGlyphsDialog::ShowUncertaintyGlyphsDialog(const std::string& name,
                                                         ModuleStateHandle state,
                                                         QWidget* parent)
: ModuleDialogGeneric(state, parent)
{
  setupUi(this);
  // setWindowTitle(QString::fromStdString(name));
  // fixSize();
}

// void ShowUncertaintyGlyphsDialog::pull()
// {
  // pull the code from the module and set in the dialog.
  // make changes necessary.
  // pull_newVersionToReplaceOld();
// }
