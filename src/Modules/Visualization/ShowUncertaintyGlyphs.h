#ifndef MODULES_VISUALIZATION_SHOW_UNCERTAINTY_GLYPHS_H_
#define MODULES_VISUALIZATION_SHOW_UNCERTAINTY_GLYPHS_H_

#include <Dataflow/Network/GeometryGeneratingModule.h>
#include <Modules/Visualization/share.h>

namespace SCIRun {
namespace Modules {
namespace Visualization {
  class SCISHARE ShowUncertaintyGlyphs : public SCIRun::Dataflow::Networks::GeometryGeneratingModule,
                                         public Has1InputPort<DynamicPortTag<FieldPortTag>>,
                                         public Has2OutputPorts<FieldPortTag, GeometryPortTag>
  {
  public:
    ShowUncertaintyGlyphs();
    virtual void execute();
    virtual void setStateDefaults();

    HAS_DYNAMIC_PORTS
    INPUT_PORT_DYNAMIC(0, InputFields, Field);
    OUTPUT_PORT(0, Mean, Field);
    OUTPUT_PORT(1, OutputGeom, GeometryObject);

    MODULE_TRAITS_AND_INFO(ModuleHasUIAndAlgorithm);

  private:
    boost::shared_ptr<class GlyphBuilder> builder_;

    int verifyData(const std::vector<FieldHandle> &fields,
                   const std::vector<std::vector<int>> &indices,
                   int fieldCount) const;
    Core::Geometry::Tensor computeMeanTensor(const std::vector<FieldHandle> &fields, int index) const;
    std::vector<Core::Geometry::Tensor> computeMeanTensors(const std::vector<FieldHandle> &fields,
                                                           int requiredSize) const;
    FieldHandle createOutputField(const std::vector<Core::Geometry::Point> &points,
                                  const std::vector<Core::Geometry::Tensor> &meanTensors,
                                  int requiredSize) const;
  };
}}}

#endif
