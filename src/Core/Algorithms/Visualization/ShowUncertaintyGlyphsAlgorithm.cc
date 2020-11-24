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

#include <Core/Algorithms/Base/AlgorithmPreconditions.h>
#include <Core/Algorithms/Base/AlgorithmVariableNames.h>
#include <Core/Algorithms/Visualization/ShowUncertaintyGlyphsAlgorithm.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Datatypes/Dyadic3DTensor.h>
#include <Core/Datatypes/Geometry.h>
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/Mesh/MeshFacade.h>
#include <Core/Logging/Log.h>
#include <Graphics/Datatypes/GeometryImpl.h>
#include <Graphics/Glyphs/GlyphConstructor.h>
#include <Graphics/Glyphs/GlyphGeom.h>
#include <Graphics/Glyphs/GlyphGeomUtility.h>
#include <Graphics/Glyphs/TensorGlyphBuilder.h>

using namespace SCIRun;
using namespace Graphics;
using namespace Graphics::Datatypes;
using namespace Core;
using namespace Core::Datatypes;
using namespace Core::Algorithms;
using namespace Core::Geometry;
using namespace Core::Algorithms::Visualization;

using GGU = GlyphGeomUtility;

const AlgorithmOutputName ShowUncertaintyGlyphsAlgorithm::MeanTensorField("MeanTensorField");
const AlgorithmOutputName ShowUncertaintyGlyphsAlgorithm::OutputGeom("OutputGeom");

class ShowUncertaintyGlyphsImpl
{
 public:
  ShowUncertaintyGlyphsImpl();
  GeometryHandle run(const GeometryIDGenerator& idGen, const FieldList& fields);
  FieldHandle getMeanTensors() const;

 private:
  void computeMeanTensors();
  Dyadic3DTensor computeMeanTensor(int index) const;
  Dyadic3DTensor computeMeanLinearInvariant(int t) const;
  void computeCovarianceMatrices();
  void verifyData(const FieldList& fields);
  void getPoints(const FieldList& fields);
  void getPointsForFields(FieldHandle field, std::vector<int>& indices, std::vector<Point>& points);
  void getTensors(const FieldList& fields);
  void computeOffsetSurface();

  GlyphConstructor constructor_;
  int fieldCount_ = 0;
  int fieldSize_ = 0;

  std::vector<std::vector<Point>> points_;
  std::vector<std::vector<int>> indices_;
  std::vector<std::vector<Dyadic3DTensor>> tensors_;
  std::vector<Dyadic3DTensor> meanTensors_;
  std::vector<Eigen::Matrix<double, 6, 6>> covarianceMatrices_;

  // TODO make ui params
  double emphasis_ = 3.0;
  int resolution_ = 100;
};

ShowUncertaintyGlyphsAlgorithm::ShowUncertaintyGlyphsAlgorithm() {}

ShowUncertaintyGlyphsImpl::ShowUncertaintyGlyphsImpl()
{
  constructor_ = GlyphConstructor();
}

// TODO move to separate algo
enum FieldDataType
{
  node,
  edge,
  face,
  cell
};

GeometryHandle ShowUncertaintyGlyphsImpl::run(
    const GeometryIDGenerator& idgen, const FieldList& fields)
{
  getPoints(fields);
  verifyData(fields);
  getTensors(fields);
  computeMeanTensors();
  computeCovarianceMatrices();
  computeOffsetSurface();

  // Creates id
  std::string idname = "ShowUncertaintyGlyphs";
  // if(!state->getValue(ShowFieldGlyphs::FieldName).toString().empty())
  // idname += GeometryObject::delimiter + state->getValue(ShowFieldGlyphs::FieldName).toString() +
  // " (from " + moduleId_ +")";

  RenderState renState;
  renState.set(RenderState::USE_NORMALS, true);
  renState.set(RenderState::IS_ON, true);
  renState.set(RenderState::USE_TRANSPARENCY, true);
  renState.mGlyphType = RenderState::GlyphType::SUPERQUADRIC_TENSOR_GLYPH;
  renState.defaultColor = ColorRGB(1.0, 1.0, 1.0);
  renState.set(RenderState::USE_DEFAULT_COLOR, true);
  SpireIBO::PRIMITIVE primIn = SpireIBO::PRIMITIVE::TRIANGLES;
  auto vmesh = fields[0]->vmesh();

  auto geom(boost::make_shared<GeometryObjectSpire>(idgen, idname, true));
  constructor_.buildObject(*geom, geom->uniqueID(), true, 0.5, ColorScheme::COLOR_UNIFORM, renState,
      primIn, vmesh->get_bounding_box(), true, nullptr);

  return geom;
}

void ShowUncertaintyGlyphsImpl::computeCovarianceMatrices()
{
  covarianceMatrices_ = std::vector<Eigen::Matrix<double, 6, 6>>(fieldSize_);

  for (int t = 0; t < fieldSize_; ++t)
  {
    covarianceMatrices_[t] = Eigen::Matrix<double, 6, 6>();
    covarianceMatrices_[t].fill(0.0);

    for (int f = 0; f < fieldCount_; ++f)
    {
      Eigen::Matrix<double, 6, 1> diffTensor = (tensors_[f][t] - meanTensors_[t]).mandel();
      covarianceMatrices_[t] += (diffTensor * diffTensor.transpose());
    }

    covarianceMatrices_[t] /= fieldCount_;
  }
}

void ShowUncertaintyGlyphsImpl::computeOffsetSurface()
{
  for (int f = 0; f < fieldSize_; ++f)
  {
    UncertaintyTensorOffsetSurfaceBuilder builder(meanTensors_[f], points_[0][f], emphasis_);
    builder.setResolution(51);
    builder.generateOffsetSurface(constructor_, covarianceMatrices_[f]);
  }
}

void ShowUncertaintyGlyphsImpl::getPoints(const FieldList& fields)
{
  fieldCount_ = fields.size();

  indices_ = std::vector<std::vector<int>>(fieldCount_);
  points_ = std::vector<std::vector<Point>>(fieldCount_);

  for (int f = 0; f < fieldCount_; ++f)
    getPointsForFields(fields[f], indices_[f], points_[f]);
}

void ShowUncertaintyGlyphsImpl::getPointsForFields(
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

void ShowUncertaintyGlyphsImpl::verifyData(const FieldList& fields)
{
  fieldSize_ = indices_[0].size();
  for (int f = 1; f < fieldCount_; ++f)
    if (indices_[f].size() != fieldSize_)
      throw std::invalid_argument("All field inputs must have the same size.");

  // Verify all are tensors
  for (auto field : fields)
  {
    FieldInformation finfo(field);
    if (!finfo.is_tensor())
      throw std::invalid_argument("Currently, this module only supports tensor fields.");
  }
}

Dyadic3DTensor scirunTensorToEigenTensor(const Geometry::Tensor& t)
{
  Dyadic3DTensor newTensor(t.xx(), t.xy(), t.xz(), t.yy(), t.yz(), t.zz());
  return newTensor;
}

void ShowUncertaintyGlyphsImpl::getTensors(const FieldList& fields)
{
  tensors_ = std::vector<std::vector<Dyadic3DTensor>>(fieldCount_);
  for (int f = 0; f < fieldCount_; ++f)
  {
    tensors_[f] = std::vector<Dyadic3DTensor>(fieldSize_);
    auto vfield = fields[f]->vfield();
    Tensor temp;
    for (int v = 0; v < fieldSize_; ++v)
    {
      vfield->get_value(temp, v);
      tensors_[f][v] = scirunTensorToEigenTensor(temp);
      tensors_[f][v].setDescendingRHSOrder();
    }
  }
}

Tensor eigenTensorToScirunTensor(const Dyadic3DTensor& t)
{
  return Tensor(t(0, 0), t(1, 0), t(2, 0), t(1, 1), t(2, 1), t(2, 2));
}

void ShowUncertaintyGlyphsImpl::computeMeanTensors()
{
  meanTensors_ = std::vector<Dyadic3DTensor>(fieldSize_);
  for (int t = 0; t < fieldSize_; ++t)
    meanTensors_[t] = computeMeanLinearInvariant(t);
    // meanTensors_[t] = computeMeanTensor(t);

  // auto eigvecs = meanTensors_[0].getEigenvectors();
  // auto eigvals = meanTensors_[0].getEigenvalues();
  // std::cout << "e0 " << e0 << "\n";
  // std::cout << "e1 " << e1 << "\n";
  // std::cout << "e2 " << e2 << "\n";
  // for (int i = 0; i < 3; ++i) std::cout << "eigvec " << i << ": " << eigvecs[i] << "\n";
  // for (int i = 0; i < 3; ++i) std::cout << "eigval " << i << ": " << eigvals[i] << "\n";
  // std::cout << "\n\n\n\n\n\n";
  // Eigen::Matrix3d diagEig = eigvals.asDiagonal();
  // std::cout << "D: " << diagEig << "\n";
  // std::cout << "V: " << meanTensors_[0].getEigenvectorsAsMatrix() << "\n";
}

Dyadic3DTensor ShowUncertaintyGlyphsImpl::computeMeanTensor(int t) const
{
  auto sum = Dyadic3DTensor();
  for (int f = 0; f < fieldCount_; ++f)
    sum += tensors_[f][t];
  sum = sum / (double)fieldCount_;

  return sum;
}

double frob(const Dyadic3DTensor& t)
{
  return sqrt((t * t.transpose()).trace());
}

// Mean calculation using Linear Invariant interpolation
Dyadic3DTensor ShowUncertaintyGlyphsImpl::computeMeanLinearInvariant(int t) const
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

  for (int f = 0; f < fieldCount_; ++f)
  {
    const double trace = tensors_[f][t].trace();
    Dyadic3DTensor trThird = trace * identityThird;
    K1 += trace;

    const double fro = tensors_[f][t].frobeniusNorm();
    Dyadic3DTensor anisotropicDeviation = tensors_[f][t] - trThird;
    anisotropicDeviation.setDescendingRHSOrder();
    const double anisotropicDeviationFro = anisotropicDeviation.frobeniusNorm();

    R2 += sqrtThreeHalves * anisotropicDeviationFro / fro;
    R3 += threeSqrtSix * (anisotropicDeviation / anisotropicDeviationFro).asMatrix().determinant();
  }

  // Equally weight all of the coeffecients
  K1 /= fieldCount_;
  R2 /= fieldCount_;
  R3 /= fieldCount_;

  const double arccosR3 = std::acos(R3);

  Eigen::Vector3d eigvals;
  const double x = oneThird * K1;
  const double y = (2.0 * K1 * R2) / (3.0 * sqrt(3.0 - 2.0 * std::pow(R2, 2.0)));
  eigvals[0] = x + y * std::cos(arccosR3 * oneThird);
  eigvals[1] = x + y * std::cos((arccosR3 - 2.0 * M_PI) * oneThird);
  eigvals[2] = x + y * std::cos((arccosR3 + 2.0 * M_PI) * oneThird);

  std::vector<Eigen::Vector3d> eigvecs = computeMeanTensor(t).getEigenvectors();

  Dyadic3DTensor ret = Dyadic3DTensor(eigvecs, eigvals);
  ret.setDescendingRHSOrder();
  ret.setTensorValues();
  return ret;
}


FieldHandle ShowUncertaintyGlyphsImpl::getMeanTensors() const
{
  FieldInformation ofinfo("PointCloudMesh", 0, "Tensor");
  auto ofield = CreateField(ofinfo);
  auto mesh = ofield->vmesh();
  auto field = ofield->vfield();

  std::vector<VMesh::index_type> meshIndices(fieldSize_);
  for (int i = 0; i < fieldSize_; ++i)
    meshIndices[i] = mesh->add_point(points_[0][i]);

  field->resize_fdata();

  for (int i = 0; i < fieldSize_; ++i)
    field->set_value(eigenTensorToScirunTensor(meanTensors_[i]), meshIndices[i]);

  return ofield;
}

// main algorithm function
AlgorithmOutput ShowUncertaintyGlyphsAlgorithm::run(
    const GeometryIDGenerator& idGen, const AlgorithmInput& input) const
{
  auto fields = input.getList<Field>(Variables::InputFields);
  if (fields.empty()) THROW_ALGORITHM_INPUT_ERROR("No input fields given");

  auto impl = ShowUncertaintyGlyphsImpl();

  AlgorithmOutput output;
  output[OutputGeom] = impl.run(idGen, fields);
  output[MeanTensorField] = impl.getMeanTensors();
  return output;
}

AlgorithmOutput ShowUncertaintyGlyphsAlgorithm::run(const AlgorithmInput& input) const
{
  auto fields = input.getList<Field>(Variables::InputFields);
  if (fields.empty()) THROW_ALGORITHM_INPUT_ERROR("No input fields given");

  auto impl = ShowUncertaintyGlyphsImpl();

  AlgorithmOutput output;
  output[MeanTensorField] = impl.getMeanTensors();
  return output;
}
