/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2020 Scientific Computing and Imaging Institute,
   University of Utah.

   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/


#include <iostream>
#include <Core/Algorithms/Math/TensorInterpolation.h>
#include <Core/Algorithms/Base/AlgorithmVariableNames.h>
#include <Core/Algorithms/Base/AlgorithmPreconditions.h>
#include <unsupported/Eigen/MatrixFunctions>
// #include <Core/Datatypes/Legacy/Field/FieldInformation.h>
// #include <Core/Datatypes/Legacy/Field/VField.h>
// #include <Core/Datatypes/Legacy/Field/VMesh.h>
// #include <Core/Logging/Log.h>

using namespace SCIRun;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Algorithms::Fields;

TensorInterpolationAlgo::TensorInterpolationAlgo(Method method) {
    invariantMethod_ = method;
    orientMethod_ = method;

  // using namespace Parameters;
  // addParameter(InvariantMethod, Method::LOG_EUCLIDEAN);
  // addParameter(OrientMethod, Method::LOG_EUCLIDEAN);
}

TensorInterpolationAlgo::TensorInterpolationAlgo(Method invariantMethod, Method orientMethod) {
    invariantMethod_ = invariantMethod;
    orientMethod_ = orientMethod;
}

Dyadic3DTensor TensorInterpolationAlgo::matrixAverage(std::vector<Dyadic3DTensor>& tensors, std::optional<std::vector<float>> weights) {
  Dyadic3DTensor ret(0.0);
  double weight_total = 0.0;
  for (int i = 0; i < tensors.size(); ++i) {
      if (weights.has_value()) {
          ret += weights.value()[i] * tensors[i];
          weight_total += weights.value()[i];
      }
      else
          ret += tensors[i];
  }

  if (weights.has_value())
      ret = ret / weight_total;
  else
      ret = ret / static_cast<double>(tensors.size());

  return ret;
}

Dyadic3DTensor TensorInterpolationAlgo::logEuclidean(std::vector<Dyadic3DTensor>& tensors, std::optional<std::vector<float>> weights) {
  Eigen::Matrix3d sum = Eigen::Matrix3d::Zero();
  double weight_total = 0.0;
  for (size_t i = 0; i < tensors.size(); ++i)
  {
    tensors[i].setDescendingRHSOrder();
    Eigen::Vector3d eigvalsLog = tensors[i].getEigenvalues();
    Eigen::Matrix3d rotation = tensors[i].getEigenvectorsAsMatrix();
    for (size_t i = 0; i < eigvalsLog.size(); ++i)
      eigvalsLog[i] = std::log(eigvalsLog[i]);

    Eigen::Matrix3d sumAdd = rotation * eigvalsLog.asDiagonal() * rotation.transpose();
    if (weights.has_value()) {
        sum += weights.value()[i] * sumAdd;
        weight_total += weights.value()[i];
    }
    else
        sum += sumAdd;
  }

  if (weights.has_value())
      sum /= weight_total;
  else
      sum /= static_cast<double>(tensors.size());

  Eigen::Matrix3d mean = sum.exp();
  Dyadic3DTensor ret = Dyadic3DTensor(mean);
  return ret;
}

Dyadic3DTensor TensorInterpolationAlgo::linearInvariant(std::vector<Dyadic3DTensor>& tensors, std::optional<std::vector<float>> weights) {
  const static double oneThird = 1.0 / 3.0;
  const static double sqrtThreeHalves = sqrt(1.5);
  const static double threeSqrtSix = 3.0 * sqrt(6.0);

  double K1 = 0.0;
  double R2 = 0.0;
  double R3 = 0.0;

  for (size_t i = 0; i < tensors.size(); ++i)
  {
    const double trace = tensors[i].trace();
    const double fro = tensors[i].frobeniusNorm();

    Dyadic3DTensor anisotropicDeviation = tensors[i] - Dyadic3DTensor(trace * oneThird);
    anisotropicDeviation.setDescendingRHSOrder();
    const double anisotropicDeviationFro = anisotropicDeviation.frobeniusNorm();

    float R2_add = sqrtThreeHalves * anisotropicDeviationFro / fro;
    float R3_add = threeSqrtSix * (anisotropicDeviation / anisotropicDeviationFro).asMatrix().determinant();
    if (weights.has_value()) {
        K1 += weights.value()[i] * trace;
        R2 += weights.value()[i] * R2_add;
        R3 += weights.value()[i] * R3_add;
    } else {
        K1 += trace;
        R2 += R2_add;
        R3 += R3_add;
    }
  }

  // Equally weight all of the coeffecients
  double div = 1.0 / tensors.size();
  K1 *= div;
  R2 *= div;
  R3 *= div;

  // Clamp to avoid nan results with acos
  if (R3 > 1.0) R3 = 1.0;
  if (R3 < -1.0) R3 = -1.0;


  const double arccosR3 = std::acos(R3);
  Eigen::Vector3d eigvals;
  const double x = oneThird * K1;
  const double y = (2.0 * K1 * R2) / (3.0 * sqrt(3.0 - 2.0 * std::pow(R2, 2.0)));
  eigvals(0) = x + y * std::cos(arccosR3 * oneThird);
  eigvals(1) = x + y * std::cos((arccosR3 - 2.0 * M_PI) * oneThird);
  eigvals(2) = x + y * std::cos((arccosR3 + 2.0 * M_PI) * oneThird);

  // Using axis aligned orientation because this mean method is only for invariants
  return Dyadic3DTensor(eigvals);
}

Dyadic3DTensor TensorInterpolationAlgo::interpolate(std::vector<Dyadic3DTensor>& tensors) {
    return interpolate(tensors, std::nullopt);
}

Dyadic3DTensor TensorInterpolationAlgo::interpolate(std::vector<Dyadic3DTensor>& tensors, std::optional<std::vector<float>> weights) {
    Dyadic3DTensor invariant, orient;
    switch (invariantMethod_) {
        case Method::MATRIX_AVERAGE:
            invariant = matrixAverage(tensors, weights);
            break;
        case Method::LOG_EUCLIDEAN:
            invariant = logEuclidean(tensors, weights);
            break;
        default:
            invariant = linearInvariant(tensors, weights);
            break;
    }

    if (invariantMethod_ == orientMethod_)
        return invariant;

    switch (orientMethod_) {
        case Method::MATRIX_AVERAGE:
            orient = matrixAverage(tensors, weights);
            break;
        case Method::LOG_EUCLIDEAN:
            orient = logEuclidean(tensors, weights);
            break;
        default:
            orient = linearInvariant(tensors, weights);
            break;
    }

    return Dyadic3DTensor(orient.getEigenvectors(), invariant.getEigenvalues());
}

// AlgorithmOutput TensorInterpolationAlgo::run_generic(const AlgorithmInput& input) const
// {
//   auto inputField = input.get<Field>(Variables::InputField);

//   FieldHandle outputField(inputField->deep_clone());
//   double knob2 = get(Parameters::Knob2).toDouble();
//   if (get(Parameters::Knob1).getBool())
//   {
//     // do something
//   }
//   else
//   {
//     // do something else
//   }

//   AlgorithmOutput output;
//   output[Variables::OutputField] = outputField;
//   return output;
// }
