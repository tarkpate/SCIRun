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


#include <Core/Algorithms/Math/DownsampleTensorFieldAlgorithm.h>
#include <Core/Algorithms/Base/AlgorithmPreconditions.h>
#include <Core/Datatypes/Dyadic3DTensor.h>
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>
#include <Core/Datatypes/Mesh/MeshFacade.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <Core/GeometryPrimitives/Point.h>
#include <Core/GeometryPrimitives/Vector.h>
#include <boost/tuple/tuple.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <algorithm>

using namespace SCIRun;
using namespace Core;
using namespace Core::Datatypes;
using namespace Core::Algorithms;
using namespace Core::Algorithms::Math;
using namespace Core::Geometry;

AlgorithmOutputName DownsampleTensorFieldAlgorithm::MeanTensorField("MeanTensorField");
AlgorithmOutputName DownsampleTensorFieldAlgorithm::VariationField("VariationField");

ALGORITHM_PARAMETER_DEF(Math, MeanTensorInvariantMethod);
ALGORITHM_PARAMETER_DEF(Math, MeanTensorOrientationMethod);

// TODO move helpers
namespace helpers {
Dyadic3DTensor scirunTensorToEigenTensor(const Geometry::Tensor& t)
{
  Dyadic3DTensor newTensor(t.xx(), t.xy(), t.xz(), t.yy(), t.yz(), t.zz());
  return newTensor;
}

Tensor eigenTensorToScirunTensor(const Dyadic3DTensor& t)
{
  return Tensor(t(0, 0), t(1, 0), t(2, 0), t(1, 1), t(2, 1), t(2, 2));
}
}

// TODO move to separate algo
enum FieldDataType
{ node, edge, face, cell };

class DownsampleTensorFieldAlgorithmImpl
{
public:
  DownsampleTensorFieldAlgorithmImpl();
  void run(const FieldHandle field);
  void setInvariantMethod(AlgoOption method);
  void setOrientationMethod(AlgoOption method);
  void getPoints(const FieldHandle field);
  void getPointsForFields(FieldHandle field, std::vector<int>& indices, std::vector<Point>& points);
  FieldHandle getMeanTensors() const;
  FieldHandle getVariationField() const;
  DenseMatrixHandle getCovarianceMatrices() const;
private:
  size_t fieldCount_ = 0;
  size_t fieldSize_ = 0;
  size_t newFieldSize_ = 0;
  size_t resFactor_ = 2;
  VMesh::dimension_type dims_, newDims_;
  std::vector<Point> points_;
  std::vector<int> indices_;
  std::vector<Point> newPoints_;
  std::vector<double> FAVariation_;
  std::vector<Dyadic3DTensor> tensors_;
  std::vector<Dyadic3DTensor> meanTensors_;
  std::vector<DyadicTensor<6>> covarianceMatrices_;
  AlgoOption invariantMethod_;
  AlgoOption orientationMethod_;

  void verifyData(const FieldHandle field);
  void getTensors(const FieldHandle field);
  void computeCovarianceMatrices();
  void computeMeanTensors();
  void computeVariation();
  Dyadic3DTensor computeMeanMatrixAverage(int x, int y, int z) const;
  Dyadic3DTensor computeMeanLinearInvariant(int t) const;
  Dyadic3DTensor computeMeanLogEuclidean(int t) const;
  void computeMeanTensorsSameMethod();
  void computeMeanTensorsDifferentMethod();
};

void DownsampleTensorFieldAlgorithmImpl::setInvariantMethod(AlgoOption method)
{
  invariantMethod_ = method;
}

void DownsampleTensorFieldAlgorithmImpl::setOrientationMethod(AlgoOption method)
{
  orientationMethod_ = method;
}

FieldHandle DownsampleTensorFieldAlgorithmImpl::getMeanTensors() const
{
  FieldInformation ofinfo("PointCloudMesh", 0, "Tensor");
  auto ofield = CreateField(ofinfo);
  auto mesh = ofield->vmesh();
  auto field = ofield->vfield();

  std::vector<VMesh::index_type> meshIndices(newFieldSize_);
  for (size_t i = 0; i < newFieldSize_; ++i)
    meshIndices[i] = mesh->add_point(newPoints_[i]);

  field->resize_fdata();

  for (size_t i = 0; i < newFieldSize_; ++i)
    field->set_value(helpers::eigenTensorToScirunTensor(meanTensors_[i]), meshIndices[i]);

  return ofield;
}

FieldHandle DownsampleTensorFieldAlgorithmImpl::getVariationField() const
{
  FieldInformation ofinfo("PointCloudMesh", 0, "Scalar");
  auto ofield = CreateField(ofinfo);
  auto mesh = ofield->vmesh();
  auto field = ofield->vfield();

  std::vector<VMesh::index_type> meshIndices(newFieldSize_);
  for (size_t i = 0; i < newFieldSize_; ++i)
    meshIndices[i] = mesh->add_point(newPoints_[i]);

  field->resize_fdata();

  for (size_t i = 0; i < newFieldSize_; ++i)
    field->set_value(FAVariation_[i], meshIndices[i]);

  return ofield;
}

void DownsampleTensorFieldAlgorithmImpl::computeMeanTensorsSameMethod()
{
  if (orientationMethod_.option_ == "Matrix Average")
  {
    for (size_t z = 0; z < newDims_[2]; ++z)
      for (size_t y = 0; y < newDims_[1]; ++y)
        for (size_t x = 0; x < newDims_[0]; ++x)
          meanTensors_[z*newDims_[0]*newDims_[1]+y*newDims_[0]+x] = computeMeanMatrixAverage(x, y, z);
  }
  else if (orientationMethod_.option_ == "Log-Euclidean")
  {
    for (size_t t = 0; t < fieldSize_/4; ++t)
      meanTensors_[t] = computeMeanLogEuclidean(t);
  }
}

void DownsampleTensorFieldAlgorithmImpl::computeMeanTensorsDifferentMethod()
{
  std::vector<Dyadic3DTensor> orientationTensors(fieldSize_);
  std::vector<Dyadic3DTensor> invariantTensors(fieldSize_);
  if (orientationMethod_.option_ == "Matrix Average")
  {
    for (size_t t = 0; t < fieldSize_; ++t)
      orientationTensors[t] = computeMeanMatrixAverage(t, 0, 0);
  }
  else if (orientationMethod_.option_ == "Log-Euclidean")
  {
    for (size_t t = 0; t < fieldSize_; ++t)
      orientationTensors[t] = computeMeanLogEuclidean(t);
  }

  if (invariantMethod_.option_ == "Matrix Average")
  {
    for (size_t t = 0; t < fieldSize_; ++t)
      invariantTensors[t] = computeMeanMatrixAverage(t, 0, 0);
  }
  else if (invariantMethod_.option_ == "Log-Euclidean")
  {
    for (size_t t = 0; t < fieldSize_; ++t)
      invariantTensors[t] = computeMeanLogEuclidean(t);
  }
  else if (invariantMethod_.option_ == "Linear Invariant")
  {
    for (size_t t = 0; t < fieldSize_; ++t)
      invariantTensors[t] = computeMeanLinearInvariant(t);
  }

  for (size_t t = 0; t < fieldSize_; ++t)
  {
    invariantTensors[t].setDescendingRHSOrder();
    orientationTensors[t].setDescendingRHSOrder();
    meanTensors_[t] = Dyadic3DTensor(orientationTensors[t].getEigenvectors(),
                                     invariantTensors[t].getEigenvalues());
    meanTensors_[t].setDescendingRHSOrder();
  }
}

void DownsampleTensorFieldAlgorithmImpl::computeMeanTensors()
{
  meanTensors_ = std::vector<Dyadic3DTensor>(newFieldSize_);
  newPoints_ = std::vector<Point>(newFieldSize_);
  if (invariantMethod_.option_ == orientationMethod_.option_)
    computeMeanTensorsSameMethod();
  else
    computeMeanTensorsDifferentMethod();

  for (size_t z = 0; z < newDims_[2]; ++z)
    for (size_t y = 0; y < newDims_[1]; ++y)
      for (size_t x = 0; x < newDims_[0]; ++x)
      {
        double div = 0.0;
        Point p;
        for (int k = z * resFactor_; k < std::min((z+1)*resFactor_,(size_t)dims_[2]); ++k)
          for (int j = y * resFactor_; j < std::min((y+1)*resFactor_,(size_t)dims_[1]); ++j)
            for (int i = x * resFactor_; i < std::min((x+1)*resFactor_,(size_t)dims_[0]); ++i) {
              p += Point(points_[k*dims_[0]*dims_[1]+j*dims_[0]+i]);
              div++;
            }
        newPoints_[z*newDims_[0]*newDims_[1]+y*newDims_[0]+x] = p / div;
      }
}

void DownsampleTensorFieldAlgorithmImpl::computeVariation()
{
  FAVariation_ = std::vector<double>(newFieldSize_);
  for (size_t z = 0; z < newDims_[2]; ++z)
    for (size_t y = 0; y < newDims_[1]; ++y)
      for (size_t x = 0; x < newDims_[0]; ++x)
      {
        double div = 0.0;
        double v = 0.0;
        for (int k = z * resFactor_; k < std::min((z+1)*resFactor_,(size_t)dims_[2]); ++k)
          for (int j = y * resFactor_; j < std::min((y+1)*resFactor_,(size_t)dims_[1]); ++j)
            for (int i = x * resFactor_; i < std::min((x+1)*resFactor_,(size_t)dims_[0]); ++i) {
              if (abs(tensors_[k*dims_[0]*dims_[1]+j*dims_[0]+i].getEigenvalue(0)) > 2e-12) {
                v += tensors_[k*dims_[0]*dims_[1]+j*dims_[0]+i].fractionalAnisotropy();
                div++;
              }
            }

        double meanV = meanTensors_[z*newDims_[0]*newDims_[1]+y*newDims_[0]+x].fractionalAnisotropy();
        FAVariation_[z*newDims_[0]*newDims_[1]+y*newDims_[0]+x] = (div > 0.0) ? (v / div) - meanV : 0.0;
      }
}

Dyadic3DTensor DownsampleTensorFieldAlgorithmImpl::computeMeanMatrixAverage(int x, int y, int z) const
{
  Dyadic3DTensor sum = Dyadic3DTensor();
  double div = 0.0;
  for (int k = z * resFactor_; k < std::min((z+1)*resFactor_,(size_t)dims_[2]); ++k)
    for (int j = y * resFactor_; j < std::min((y+1)*resFactor_,(size_t)dims_[1]); ++j)
      for (int i = x * resFactor_; i < std::min((x+1)*resFactor_,(size_t)dims_[0]); ++i) {
        if (abs(tensors_[k*dims_[0]*dims_[1]+j*dims_[0]+i].getEigenvalue(0)) > 2e-12) {
          sum += tensors_[k*dims_[0]*dims_[1]+j*dims_[0]+i];
          div++;
        }
      }
  if (div > 0.0)
    sum = sum / div;
  sum.setDescendingRHSOrder();

  return sum;
}

// Mean calculation using Linear Invariant interpolation
Dyadic3DTensor DownsampleTensorFieldAlgorithmImpl::computeMeanLinearInvariant(int t) const
{
  const static double oneThird = 1.0 / 3.0;
  const static double sqrtThreeHalves = sqrt(1.5);
  const static double threeSqrtSix = 3.0 * sqrt(6.0);

  const static Dyadic3DTensor identity = Dyadic3DTensor(Eigen::Vector3d(1, 0, 0),
                                                        Eigen::Vector3d(0, 1, 0),
                                                        Eigen::Vector3d(0, 0, 1));
  const static Dyadic3DTensor identityThird = oneThird * identity;

  double K1 = 0.0;
  double R2 = 0.0;
  double R3 = 0.0;

  { // TODO put in interpolation
    const double trace = tensors_[t].trace();
    Dyadic3DTensor trThird = trace * identityThird;
    K1 += trace;

    const double fro = tensors_[t].frobeniusNorm();
    Dyadic3DTensor anisotropicDeviation = tensors_[t] - trThird;
    anisotropicDeviation.setDescendingRHSOrder();
    const double anisotropicDeviationFro = anisotropicDeviation.frobeniusNorm();

    R2 += sqrtThreeHalves * anisotropicDeviationFro / fro;
    R3 += threeSqrtSix * (anisotropicDeviation / anisotropicDeviationFro).asMatrix().determinant();
  }

  // Equally weight all of the coeffecients
  K1 /= static_cast<double>(fieldCount_);
  R2 /= static_cast<double>(fieldCount_);
  R3 /= static_cast<double>(fieldCount_);

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
  Dyadic3DTensor ret = Dyadic3DTensor({Eigen::Vector3d({1,0,0}), Eigen::Vector3d({0,1,0}), Eigen::Vector3d({0,0,1})}, eigvals);
  return ret;
}

Dyadic3DTensor DownsampleTensorFieldAlgorithmImpl::computeMeanLogEuclidean(int t) const
{
  Eigen::Matrix3d sum = Eigen::Matrix3d::Zero();
// TODO
  {
    Dyadic3DTensor tensor = tensors_[t];
    tensor.setDescendingRHSOrder();
    Eigen::Vector3d eigvalsLog = tensor.getEigenvalues();
    Eigen::Matrix3d rotation = tensor.getEigenvectorsAsMatrix();
    for (size_t i = 0; i < eigvalsLog.size(); ++i)
      eigvalsLog[i] = std::log(eigvalsLog[i]);
    sum += rotation * eigvalsLog.asDiagonal() * rotation.transpose();
  }
  sum /= (double)fieldCount_;
  Eigen::Matrix3d mean = sum.exp();
  Dyadic3DTensor ret = Dyadic3DTensor(mean);
  return ret;
}

DownsampleTensorFieldAlgorithmImpl::DownsampleTensorFieldAlgorithmImpl()
{}

void DownsampleTensorFieldAlgorithmImpl::computeCovarianceMatrices()
{
  covarianceMatrices_ = std::vector<DyadicTensor<6>>(fieldSize_);

  Eigen::Matrix<double, 6, 6> covarianceMat;
  for (size_t t = 0; t < fieldSize_; ++t)
  {
    covarianceMat.fill(0.0);
    // TODO
    {
      Eigen::Matrix<double, 6, 1> diffTensor = (tensors_[t] - meanTensors_[t]).mandel();
      covarianceMat += diffTensor * diffTensor.transpose();
    }
    covarianceMat /= static_cast<double>(fieldCount_);
    covarianceMatrices_[t] = DyadicTensor<6>(covarianceMat);
  }
}

void DownsampleTensorFieldAlgorithmImpl::getPointsForFields(
    FieldHandle field, std::vector<int>& indices, std::vector<Point>& points)
{
  // Collect indices and points from facades
  FieldDataType fieldLocation;
  FieldInformation finfo(field);
  if (finfo.is_point() || finfo.is_linear())
    fieldLocation = FieldDataType::node;
  else if (finfo.is_line())
    fieldLocation = FieldDataType::edge;
  else if (finfo.is_surface())
    fieldLocation = FieldDataType::face;
  else
    fieldLocation = FieldDataType::cell;

  auto mesh = field->vmesh();
  auto primaryFacade = field->mesh()->getFacade();
  switch (fieldLocation)
  {
  case FieldDataType::node:
    for (const auto& node : primaryFacade->nodes())
    {
      indices.push_back(node.index());
      Point p;
      mesh->get_center(p, node.index());
      points.push_back(p);
    }
    break;
  case FieldDataType::edge:
    for (const auto& edge : primaryFacade->edges())
    {
      indices.push_back(edge.index());
      Point p;
      mesh->get_center(p, edge.index());
      points.push_back(p);
    }
    break;
  case FieldDataType::face:
    for (const auto& face : primaryFacade->faces())
    {
      indices.push_back(face.index());
      Point p;
      mesh->get_center(p, face.index());
      points.push_back(p);
    }
    break;
  case FieldDataType::cell:
    for (const auto& cell : primaryFacade->cells())
    {
      indices.push_back(cell.index());
      Point p;
      mesh->get_center(p, cell.index());
      points.push_back(p);
    }
    break;
  }
}

void DownsampleTensorFieldAlgorithmImpl::getPoints(const FieldHandle field)
{
  // fieldCount_ = fields.size();

  indices_ = std::vector<int>();
  points_ = std::vector<Point>();

getPointsForFields(field, indices_, points_);
}

void DownsampleTensorFieldAlgorithmImpl::verifyData(const FieldHandle field)
{
  // Verify field is of tensors
    FieldInformation finfo(field);
    if (!finfo.is_tensor())
      THROW_ALGORITHM_INPUT_ERROR_SIMPLE("This module only supports tensor fields.");
}

void DownsampleTensorFieldAlgorithmImpl::getTensors(const FieldHandle field)
{
  tensors_ = std::vector<Dyadic3DTensor>(fieldSize_);
  auto vfield = field->vfield();
  for (size_t v = 0; v < fieldSize_; ++v)
  {
  Tensor temp;
    vfield->get_value(temp, v);
    // temp.build_eigens_from_mat();
    tensors_[v] = helpers::scirunTensorToEigenTensor(temp);
    // tensors_[v].makePositive();
    // tensors_[v].setDescendingRHSOrder();
    // std::cout << "ten " << v << ": " << tensors_[v] << std::endl;
  }
}

DenseMatrixHandle DownsampleTensorFieldAlgorithmImpl::getCovarianceMatrices() const
{
  auto m = std::make_shared<DenseMatrix>(21, fieldSize_);
  for (size_t i = 0; i < fieldSize_; ++i)
  {
    m->col(i) = covarianceMatrices_[i].mandel();
    Eigen::Matrix<double,21,1> recov_m = m->col(i);
    DyadicTensor<6> recov = DyadicTensor<6>(recov_m);
  }
  return m;
}

void DownsampleTensorFieldAlgorithmImpl::run(const FieldHandle field)
{
  VMesh* imesh = field->vmesh();
  imesh->get_dimensions(dims_);
  imesh->get_dimensions(newDims_);
  for (int i = 0; i < 3; ++i)
    if (dims_[i] == 0)
      dims_[i] = newDims_[i] = 1;
  for (int i = 0; i < 3; ++i)
    if (newDims_[i] > 1)
      newDims_[i] /= resFactor_;

  std::cout <<"dims: " << dims_[0] << ", "<< dims_[1] << ", "<< dims_[2] << "\n";
  std::cout <<"new dims: " << newDims_[0] << ", "<< newDims_[1] << ", "<< newDims_[2] << "\n";
  fieldSize_ = dims_[0]*dims_[1]*dims_[2];
  newFieldSize_ = newDims_[0]*newDims_[1]*newDims_[2];
  getPoints(field);
  verifyData(field);
  getTensors(field);
  computeMeanTensors();
  computeVariation();
  // computeCovarianceMatrices();
}

DownsampleTensorFieldAlgorithm::DownsampleTensorFieldAlgorithm()
{
  addOption(Parameters::MeanTensorInvariantMethod, "Linear Invariant",
            "Linear Invariant|Log-Euclidean|Matrix Average");
  addOption(Parameters::MeanTensorOrientationMethod, "Log-Euclidean",
            "Log-Euclidean|Matrix Average");
}

AlgorithmOutput DownsampleTensorFieldAlgorithm::run(const AlgorithmInput& input) const
{
  auto field = input.get<Field>(Variables::InputField);
  auto invariantMethod = get(Parameters::MeanTensorInvariantMethod).toOption();
  auto orientMethod = get(Parameters::MeanTensorOrientationMethod).toOption();
  auto data = runImpl(field, invariantMethod, orientMethod);

  AlgorithmOutput output;
  output[MeanTensorField] = data.get<0>();
  output[VariationField] = data.get<1>();
  return output;
}

boost::tuple<FieldHandle, FieldHandle> DownsampleTensorFieldAlgorithm::runImpl(const FieldHandle field, const AlgoOption& invariantMethod, const AlgoOption& orientationOption) const
{
  // if (field.empty())
    // THROW_ALGORITHM_INPUT_ERROR("No input fields given");

  auto impl = DownsampleTensorFieldAlgorithmImpl();
  impl.setInvariantMethod(invariantMethod);
  impl.setOrientationMethod(orientationOption);
  impl.run(field);

  // return boost::make_tuple(impl.getMeanTensors(), impl.getCovarianceMatrices());
  return boost::make_tuple(impl.getMeanTensors(), impl.getVariationField());
}
