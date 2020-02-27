#ifndef CORE_ALGORITHMS_FIELDS_SHOW_UNCERTAINTY_GLYPHS_ALGORITHM_H
#define CORE_ALGORITHMS_FIELDS_SHOW_UNCERTAINTY_GLYPHS_ALGORITHM_H

#include <Core/Algorithms/Base/AlgorithmBase.h>
#include <Core/Algorithms/Visualization/share.h>
#include <Core/GeometryPrimitives/Tensor.h>
#include <Core/GeometryPrimitives/Point.h>

namespace SCIRun {
namespace Core {
namespace Algorithms {
namespace Visualization {
  class SCISHARE ShowUncertaintyGlyphsAlgorithm : public AlgorithmBase
  {
  public:
    ShowUncertaintyGlyphsAlgorithm();
    virtual AlgorithmOutput run(const AlgorithmInput& input);
  private:
    void computeMeanTensors();
    Core::Geometry::Tensor computeMeanTensor(int index) const;
    void verifyData(const FieldList &fields);
    void getPoints(const FieldList &fields);
    void getPointsForFields(FieldHandle field, std::vector<int> &indices,
                            std::vector<Core::Geometry::Point> &points);
    void getTensors(FieldList &fields);
    FieldHandle createOutputField() const;


    int fieldCount_ = 0;
    int fieldSize_ = 0;

    std::vector<std::vector<Core::Geometry::Point>> points_;
    std::vector<std::vector<int>> indices_;
    std::vector<std::vector<Core::Geometry::Tensor>> tensors_;
    std::vector<Core::Geometry::Tensor> meanTensors_;
  };
}}}}

#endif
