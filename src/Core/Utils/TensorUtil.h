#ifndef TENSORUTIL_H_
#define TENSORUTIL_H_

#include <Core/Datatypes/Dyadic3DTensor.h>
#include <Core/GeometryPrimitives/Tensor.h>

using namespace SCIRun;
using namespace Core;
using namespace Datatypes;
using namespace Geometry;

namespace SCIRun {
namespace Core {

Dyadic3DTensor scirunTensorToEigenTensor(const Tensor& t)
{
  Dyadic3DTensor newTensor(t.xx(), t.xy(), t.xz(), t.yy(), t.yz(), t.zz());
  return newTensor;
}

Tensor eigenTensorToScirunTensor(const Dyadic3DTensor& t)
{
  return Tensor(t(0, 0), t(1, 0), t(2, 0), t(1, 1), t(2, 1), t(2, 2));
}

}}

#endif // TENSORUTIL_H_
