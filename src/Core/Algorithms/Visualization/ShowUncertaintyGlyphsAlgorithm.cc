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


#include <Core/Algorithms/Visualization/ShowUncertaintyGlyphsAlgorithm.h>
#include <Core/Algorithms/Base/AlgorithmVariableNames.h>
#include <Core/Algorithms/Base/AlgorithmPreconditions.h>
#include <Core/Logging/Log.h>
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/Mesh/MeshFacade.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Graphics/Glyphs/GlyphGeom.h>
#include <Core/Datatypes/Geometry.h>
#include <Graphics/Datatypes/GeometryImpl.h>
#include <Graphics/Glyphs/TensorGlyphBuilder.h>
#include <Graphics/Glyphs/GlyphConstructor.h>
#include <Graphics/Glyphs/GlyphGeomUtility.h>

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
  GeometryHandle run(const GeometryIDGenerator& idGen, const FieldList &fields);
  FieldHandle getMeanTensors() const;
private:
  void computeMeanTensors();
  Tensor computeMeanTensor(int index) const;
  void computeCovarianceMatrices();
  void verifyData(const FieldList &fields);
  void getPoints(const FieldList &fields);
  void getPointsForFields(FieldHandle field, std::vector<int> &indices, std::vector<Point> &points);
  void getTensors(const FieldList &fields);
  void computeOffsetSurface();
  double evaluateSuperquadricImpl(bool linear, const Point& p, double A, double B);
  double evaluateSuperquadricImplLinear(const Point& p, double A, double B);
  double evaluateSuperquadricImplPlanar(const Point& p, double A, double B);
  double diffT(const TensorGlyphBuilder& builder, const DenseMatrix& mean, const
               DenseMatrix& finiteDiff, const Point& p);

  GlyphConstructor constructor_;
  int fieldCount_ = 0;
  int fieldSize_ = 0;

  std::vector<std::vector<Point>> points_;
  std::vector<std::vector<int>> indices_;
  std::vector<std::vector<Tensor>> tensors_;
  std::vector<Tensor> meanTensors_;
  std::vector<DenseMatrix> covarianceMatrices_;

  //TODO make ui params
  double emphasis_ = 3.0;
  int resolution_ = 100;
};

ShowUncertaintyGlyphsAlgorithm::ShowUncertaintyGlyphsAlgorithm()
{
}

ShowUncertaintyGlyphsImpl::ShowUncertaintyGlyphsImpl()
{
  constructor_ = GlyphConstructor();
}

//TODO move to separate algo
enum FieldDataType
{
  node,
  edge,
  face,
  cell
};

GeometryHandle ShowUncertaintyGlyphsImpl::run(const GeometryIDGenerator& idgen, const FieldList &fields)
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
    // idname += GeometryObject::delimiter + state->getValue(ShowFieldGlyphs::FieldName).toString() + " (from " + moduleId_ +")";

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
  constructor_.buildObject(*geom, geom->uniqueID(), true, 0.5,
                           ColorScheme::COLOR_UNIFORM, renState, primIn, vmesh->get_bounding_box(),
                           true, nullptr);

  return geom;
}

void ShowUncertaintyGlyphsImpl::computeCovarianceMatrices()
{
  covarianceMatrices_ = std::vector<DenseMatrix>(fieldSize_);

  for (int t = 0; t < fieldSize_; ++t)
  {
    covarianceMatrices_[t] = DenseMatrix(6, 6, 0.0);

    for (int f = 0; f < fieldCount_; ++f)
    {
      DenseColumnMatrix diffTensor = (tensors_[f][t] - meanTensors_[t]).mandel();
      covarianceMatrices_[t] += (diffTensor * diffTensor.transpose());
    }

    covarianceMatrices_[t] /= fieldCount_;
  }
}

double ShowUncertaintyGlyphsImpl::evaluateSuperquadricImpl(bool linear, const Point& p,
                                                           double A, double B)
{
  if (linear)
    return evaluateSuperquadricImplLinear(p, A, B);
  else
    return evaluateSuperquadricImplPlanar(p, A, B);
}

// Generate around x-axis
double ShowUncertaintyGlyphsImpl::evaluateSuperquadricImplLinear(const Point& p, double A, double B)
{
  double twoDivA = 2.0 / A;
  double twoDivB = 2.0 / B;
  return GGU::spow(GGU::spow(std::abs(p.y()), twoDivA) + GGU::spow(std::abs(p.z()), twoDivA), A/B)
    + GGU::spow(std::abs(p.x()), twoDivB) - 1;
}

// Generate around z-axis
double ShowUncertaintyGlyphsImpl::evaluateSuperquadricImplPlanar(const Point& p, double A, double B)
{
  double twoDivA = 2.0 / A;
  double twoDivB = 2.0 / B;
  return GGU::spow(GGU::spow(std::abs(p.x()), twoDivA) + GGU::spow(std::abs(p.y()), twoDivA), A/B)
    + GGU::spow(std::abs(p.z()), twoDivB) - 1;
}

double ShowUncertaintyGlyphsImpl::diffT(const TensorGlyphBuilder& builder, const DenseMatrix& t1,
                                        const DenseMatrix& t2, const Point& p)
{
  std::vector<double> dist(2);
  for (int i = 0; i < 2; ++i)
  {
    DenseColumnMatrix tMandel = (i == 0) ? t1 : t2;

    TensorGlyphBuilder newBuilder = builder;
    Tensor newT = symmetricTensorFromMandel(tMandel);
    newBuilder.setTensor(newT);
    // std::cout << "t " << symmetricTensorFromMandel(t) << "\n";
    newBuilder.makeTensorPositive(true, false);
    newBuilder.computeTransforms();
    newBuilder.postScaleTransorms();
    bool linear = newBuilder.isLinear();
    newBuilder.computeAAndB(emphasis_);
    double A = newBuilder.getA();
    double B = newBuilder.getB();
    std::vector<Vector> eigvecs(3);
    newBuilder.getTensor().get_eigenvectors(eigvecs[0], eigvecs[1], eigvecs[2]);
    const Transform rotate = Transform(Point(0,0,0), eigvecs[0], eigvecs[1], eigvecs[2]);//newBuilder.getRotate();
    const Transform scale = newBuilder.getScale();
    auto rotateInv = rotate;
    rotateInv.invert();
    auto scaleInv = scale;
    scaleInv.invert();
    // std::cout << "scale " << scale << "\n";
    // std::cout << "scaleInv " << scaleInv << "\n";

    Tensor t = newBuilder.getTensor();
    Vector pseudoInv = Vector();
    t.get_eigenvalues(pseudoInv[0], pseudoInv[1], pseudoInv[2]);
    for (int i = 0; i < 3; ++i)
      pseudoInv[i] = 1/pseudoInv[i];
    // pseudoInv /= t.magnitude();
    // std::cout << "pesudoInv " << pseudoInv << "\n";

    Vector newP = (pseudoInv * Vector(rotateInv * p));
    dist[i] = evaluateSuperquadricImpl(linear, Point(newP), A, B);
  }

  return dist[0] - dist[1];
}

void ShowUncertaintyGlyphsImpl::computeOffsetSurface()
{
  const static double h = 0.000001;
  const static double hHalf = 0.5 * h;
  Point origin = Point(0,0,0);

  for (int f = 0; f < fieldSize_; ++f)
  {
    TensorGlyphBuilder builder(meanTensors_[f], points_[0][f]);
    builder.setResolution(resolution_);
    builder.makeTensorPositive(true, false);
    builder.computeTransforms();
    builder.postScaleTransorms();
    builder.computeSinCosTable(false);
    bool linear = builder.isLinear();
    builder.computeAAndB(emphasis_);

    auto eigvecs = builder.getEigenVectors();
    Vector cross = Cross(eigvecs[0], eigvecs[1]);
    bool flipVertices = Dot(cross, eigvecs[2]) < 2e-12;

    // std::cout << "t " << builder.getTensor() << "\n";
    // Vector e1, e2, e3;
    // builder.getTensor().get_eigenvectors(e1, e2, e3);
    // std::cout << "e1 " << e1 << "\n";
    // std::cout << "e2 " << e2 << "\n";
    // std::cout << "e3 " << e3 << "\n";
    // std::cout << "tMandel\n " << builder.getTensor().mandel() << "\n";
    DenseColumnMatrix tMandel = builder.getTensor().mandel();

    const Transform trans = builder.getTrans();
    // const Transform rotate = builder.getRotate();
    const Transform rotate = Transform(Point(0,0,0), eigvecs[0], eigvecs[1], eigvecs[2]);
    const Transform scale = builder.getScale();
    const Transform scaleThenRotate = rotate * scale;
    Transform rotateInv = rotate;
    rotateInv.invert();
    std::cout << "rotate " << rotate << "\n";
    // std::cout << "rotateInv " << rotateInv << "\n";

    SuperquadricPointParams params;
    params.A = builder.getA();
    params.B = builder.getB();

    SuperquadricPointParams normalParams;
    normalParams.A = 2.0 - params.A;
    normalParams.B = 2.0 - params.B;

    DenseColumnMatrix finiteDiff(6);
    finiteDiff.fill(0.0);
    double fro = meanTensors_[f].magnitude();

    int nv = resolution_;
    int nu = resolution_ + 1;
    for (int v=0; v < nv-1; v++)
    {
      double sinPhi[2];
      sinPhi[0] = builder.computeSinPhi(v);
      sinPhi[1] = builder.computeSinPhi(v+1);

      double cosPhi[2];
      cosPhi[0] = builder.computeCosPhi(v);
      cosPhi[1] = builder.computeCosPhi(v+1);

      for (int u=0; u<nu; u++)
      {
        constructor_.setOffset();
        params.sinTheta = normalParams.sinTheta = builder.computeSinTheta(u);
        params.cosTheta = normalParams.cosTheta = builder.computeCosTheta(u);

        for(int i = 0; i < 2; ++i)
        {
          params.sinPhi = normalParams.sinPhi = sinPhi[i];
          params.cosPhi = normalParams.cosPhi = cosPhi[i];

          Vector p = scaleThenRotate * Vector(builder.evaluateSuperquadricPoint(linear, params));
          Vector pseudoInv = Vector();
          meanTensors_[f].get_eigenvalues(pseudoInv[0], pseudoInv[1], pseudoInv[2]);
          for (int i = 0; i < 3; ++i)
            pseudoInv[i] = 1/pseudoInv[i];
          // pseudoInv /= meanTensors_[f].magnitude();

          // Surface Derivative
          DenseColumnMatrix nn(3);
          nn.fill(0.0);
          Vector diff = Vector(0.0, 0.0, 0.0);
          for (int j = 0; j < 3; ++j)
          {
            diff[j] = hHalf;
            Point newP = Point(pseudoInv * (rotateInv * (p+diff)));
            double d1 = evaluateSuperquadricImpl(linear, newP, params.A, params.B);

            newP = Point(pseudoInv * (rotateInv * (p-diff)));
            double d2 = evaluateSuperquadricImpl(linear, newP, params.A, params.B);
            diff[j] = 0.0;
            nn(j) = (d1 - d2) / h;
          }

          DenseColumnMatrix qn(6);
          qn.fill(0.0);
          // Vector pNorm = rotate * Vector(builder.evaluateSuperquadricPoint(linear, params));
          for (int j = 0; j < 6; ++j)
          {
            finiteDiff(j) = hHalf;
            qn(j) = diffT(builder, tMandel + finiteDiff, tMandel - finiteDiff, Point(p));
            finiteDiff(j) = 0.0;
          }
          qn /= h;
          qn /= nn.norm();

          double q = std::sqrt(std::abs((qn.transpose().eval() * (covarianceMatrices_[f] * qn)
                                         .eval()).eval().value()));
          DenseColumnMatrix n = nn/nn.norm();
          Vector nVector = Vector(n(0), n(1), n(2));

          Vector normal = Vector(builder.evaluateSuperquadricPoint(linear, normalParams));
          normal = rotate * normal;
          normal.safe_normalize();

          Vector offsetP = p + q * normal;

          constructor_.addVertex(offsetP, nVector, ColorRGB(1.0, 1.0, 1.0));
        }

        if (flipVertices)
        {
          constructor_.addIndicesToOffset(0, 2, 1);
          constructor_.addIndicesToOffset(2, 3, 1);
        }
        else
        {
          constructor_.addIndicesToOffset(0, 1, 2);
          constructor_.addIndicesToOffset(2, 1, 3);
        }
      }
    }
    constructor_.popIndicesNTimes(6);
  }
}

void ShowUncertaintyGlyphsImpl::getPoints(const FieldList &fields)
{
  fieldCount_ = fields.size();

  indices_ = std::vector<std::vector<int>>(fieldCount_);
  points_  = std::vector<std::vector<Point>>(fieldCount_);

  for (int f = 0; f < fieldCount_; ++f)
    getPointsForFields(fields[f], indices_[f], points_[f]);
}

void ShowUncertaintyGlyphsImpl::getPointsForFields(FieldHandle field, std::vector<int> &indices,
                                                   std::vector<Point> &points)
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

  auto mesh =  field->vmesh();
  auto primaryFacade = field->mesh()->getFacade();
  switch(fieldLocation)
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

void ShowUncertaintyGlyphsImpl::verifyData(const FieldList &fields)
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

void ShowUncertaintyGlyphsImpl::getTensors(const FieldList &fields)
{
  tensors_ = std::vector<std::vector<Tensor>>(fieldCount_);
  for (int f = 0; f < fieldCount_; ++f)
  {
    tensors_[f] = std::vector<Tensor>(fieldSize_);
    auto vfield = fields[f]->vfield();
    for (int v = 0; v < fieldSize_; ++v)
      vfield->get_value(tensors_[f][v], v);
  }
}

void ShowUncertaintyGlyphsImpl::computeMeanTensors()
{
  meanTensors_ = std::vector<Tensor>(fieldSize_);
  for (int t = 0; t < fieldSize_; ++t)
    meanTensors_[t] = computeMeanTensor(t);
}

Tensor ShowUncertaintyGlyphsImpl::computeMeanTensor(int t) const
{
  Tensor sum = Tensor();
  for(int f = 0; f < fieldCount_; ++f)
    sum += tensors_[f][t];
  sum = sum / (double)fieldCount_;

  return sum;
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
    field->set_value(meanTensors_[i], meshIndices[i]);

  return ofield;
}

//main algorithm function
AlgorithmOutput ShowUncertaintyGlyphsAlgorithm::run(const GeometryIDGenerator& idGen, const AlgorithmInput &input) const
{
  auto fields = input.getList<Field>(Variables::InputFields);
  if (fields.empty())
    THROW_ALGORITHM_INPUT_ERROR("No input fields given");

  auto impl = ShowUncertaintyGlyphsImpl();

  AlgorithmOutput output;
  output[OutputGeom] = impl.run(idGen, fields);
  output[MeanTensorField] = impl.getMeanTensors();
  return output;
}

AlgorithmOutput ShowUncertaintyGlyphsAlgorithm::run(const AlgorithmInput &input) const
{
  auto fields = input.getList<Field>(Variables::InputFields);
  if (fields.empty())
    THROW_ALGORITHM_INPUT_ERROR("No input fields given");

  auto impl = ShowUncertaintyGlyphsImpl();

  AlgorithmOutput output;
  output[MeanTensorField] = impl.getMeanTensors();
  return output;
}
