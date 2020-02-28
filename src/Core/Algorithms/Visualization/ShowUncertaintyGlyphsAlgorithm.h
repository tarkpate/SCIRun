#ifndef CORE_ALGORITHMS_FIELDS_SHOW_UNCERTAINTY_GLYPHS_ALGORITHM_H
#define CORE_ALGORITHMS_FIELDS_SHOW_UNCERTAINTY_GLYPHS_ALGORITHM_H

#include <Core/Algorithms/Base/AlgorithmBase.h>
#include <Core/Algorithms/Visualization/share.h>

namespace SCIRun {
namespace Core {
namespace Algorithms {
namespace Visualization {
  class SCISHARE ShowUncertaintyGlyphsAlgorithm : public AlgorithmBase
  {
  public:
    ShowUncertaintyGlyphsAlgorithm();
    virtual AlgorithmOutput run(const AlgorithmInput& input) const override;

    static const AlgorithmOutputName MeanTensorField;
  };
}}}}

#endif
