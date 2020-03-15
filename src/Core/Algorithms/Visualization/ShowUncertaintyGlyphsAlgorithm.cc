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

using namespace SCIRun;
using namespace Graphics;
using namespace Core::Datatypes;
using namespace Core::Algorithms;
using namespace Core::Geometry;
using namespace Core::Algorithms::Visualization;

const AlgorithmOutputName ShowUncertaintyGlyphsAlgorithm::MeanTensorField("MeanTensorField");

class ShowUncertaintyGlyphsImpl
{
public:
  ShowUncertaintyGlyphsImpl();
  void run(const FieldList &fields);
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
  Point getSuperquadricTensorPoint(int u, int v, Transform& transform,
                                   SinCosTable& tab1, SinCosTable& tab2,
                                   double A, double B, bool linear);
  Vector getSuperquadricTensorNormal(int u, int v, Transform& rotate,
                                     SinCosTable& tab1, SinCosTable& tab2,
                                     double A, double B, bool linear);
  void makeTensorPositive(Tensor& t);

  int fieldCount_ = 0;
  int fieldSize_ = 0;

  std::vector<std::vector<Point>> points_;
  std::vector<std::vector<int>> indices_;
  std::vector<std::vector<Tensor>> tensors_;
  std::vector<Tensor> meanTensors_;
  std::vector<DenseMatrix> covarianceMatrices_;

  //TODO make ui params
  double emphasis_ = 3.0;
  int resolution_ = 10;
};

ShowUncertaintyGlyphsAlgorithm::ShowUncertaintyGlyphsAlgorithm()
{
}

ShowUncertaintyGlyphsImpl::ShowUncertaintyGlyphsImpl()
{
}

//TODO move to separate algo
enum FieldDataType
{
  node,
  edge,
  face,
  cell
};

void ShowUncertaintyGlyphsImpl::run(const FieldList &fields)
{
  getPoints(fields);
  verifyData(fields);
  getTensors(fields);
  computeMeanTensors();
  computeCovarianceMatrices();
  computeOffsetSurface();
}

void ShowUncertaintyGlyphsImpl::computeCovarianceMatrices()
{
  covarianceMatrices_ = std::vector<DenseMatrix>(fieldSize_);

  for(int t = 0; t < fieldSize_; ++t)
  {
    auto cov = DenseMatrix(6, 6, 0.0);

    for(int f = 0; f < fieldCount_; ++f)
    {
      auto diffTensor = (tensors_[f][t] - meanTensors_[t]).mandel();
      cov += diffTensor * diffTensor.transpose();
    }

    cov /= fieldCount_;
    covarianceMatrices_[t] = cov;
  }
}

void reorderTensor(std::vector<Vector>& eigvectors, std::vector<double>& eigvals)
{
  std::vector<int> indices = {0, 1, 2};
  std::vector<std::pair<double, Vector>> sortList =
    { std::make_pair(eigvals[0], eigvectors[0]),
      std::make_pair(eigvals[1], eigvectors[1]),
      std::make_pair(eigvals[2], eigvectors[2]) };

  std::sort(std::begin(sortList), std::end(sortList));

  for(int i = 0; i < 3; ++i)
  {
    eigvals[i] = sortList[i].first;
    eigvectors[i] = sortList[i].second;
  }
}

inline double spow(double e, double x)
{
  // This for round off of very small numbers.
  if( abs( e ) < 1.0e-6)
    e = 0.0;

  if (e < 0.0)
  {
    return (double)(pow(abs(e), x) * -1.0);
  }
  else
  {
    return (double)(pow(e, x));
  }
}

Point ShowUncertaintyGlyphsImpl::getSuperquadricTensorPoint(int u, int v, Transform& transform,
                                                            SinCosTable& tab1, SinCosTable& tab2,
                                                            double A, double B, bool linear)
{
  double nr = tab2.sin(v);
  double nz = tab2.cos(v);

  double nx = tab1.sin(u);
  double ny = tab1.cos(u);

  double x, y, z;
  if (linear) // Generate around x-axis
  {
    x =  spow(nz, B);
    y = -spow(nr, B) * spow(ny, A);
    z =  spow(nr, B) * spow(nx, A);
  }
  else       // Generate around z-axis
  {
    x = spow(nr, B) * spow(nx, A);
    y = spow(nr, B) * spow(ny, A);
    z = spow(nz, B);
  }
  return Point(transform * Point(x, y, z));
}

// TODO delete if not used later
Vector ShowUncertaintyGlyphsImpl::getSuperquadricTensorNormal(int u, int v, Transform& rotate,
                                                              SinCosTable& tab1, SinCosTable& tab2,
                                                              double A, double B, bool linear)
{
  double nr = tab2.sin(v);
  double nz = tab2.cos(v);

  double nx = tab1.sin(u);
  double ny = tab1.cos(u);

  double nnx, nny, nnz;
  if(linear)
  {
    nnx =  spow(nz, 2.0-B);
    nny = -spow(nr, 2.0-B) * spow(ny, 2.0-A);
    nnz =  spow(nr, 2.0-B) * spow(nx, 2.0-A);
  }
  else
  {
    nnx = spow(nr, 2.0-B) * spow(nx, 2.0-A);
    nny = spow(nr, 2.0-B) * spow(ny, 2.0-A);
    nnz = spow(nz, 2.0-B);
  }
  Vector normal = rotate * Vector(nnx, nny, nnz);
  normal.safe_normalize();
  return normal;
}

void ShowUncertaintyGlyphsImpl::makeTensorPositive(Tensor& t)
{
  static const double zeroThreshold = 0.000001;
  std::vector<Vector> eigvecs(3);
  t.get_eigenvectors(eigvecs[0], eigvecs[1], eigvecs[2]);

  double eigval1, eigval2, eigval3;
  t.get_eigenvalues(eigval1, eigval2, eigval3);
  std::vector<double> eigvals = {fabs(eigval1), fabs(eigval2), fabs(eigval3)};
  reorderTensor(eigvecs, eigvals);

  for (auto& e : eigvals)
    if(e <= zeroThreshold)
      e = 0;

  bool eig_x_0 = eigvals[0] == 0;
  bool eig_y_0 = eigvals[1] == 0;
  bool eig_z_0 = eigvals[2] == 0;

  bool flatTensor = (eig_x_0 + eig_y_0 + eig_z_0) >= 1;
  if (flatTensor)
  {
    // Check for zero eigenvectors
    if (eig_x_0)
      eigvecs[0] = Cross(eigvecs[1], eigvecs[2]);
    else if (eig_y_0)
      eigvecs[1] = Cross(eigvecs[0], eigvecs[2]);
    else if (eig_z_0)
      eigvecs[2] = Cross(eigvecs[0], eigvecs[1]);
  }

  t.set_outside_eigens(eigvecs[0], eigvecs[1], eigvecs[2],
                       eigvals[0], eigvals[1], eigvals[2]);
}

void ShowUncertaintyGlyphsImpl::diffT(TensorGlyphBuilder& builder, const DenseMatrix& mean, const
                                      DenseMatrix& finiteDiff, SuperquadricPointParams& params)
{
  auto t1 = symmetricTensorFromMandel(mean - finiteDiff);
  auto t2 = symmetricTensorFromMandel(mean + finiteDiff);

  std::vector<Point> points;
  for (const auto& t : {t1, t2})
  {
    builder.setTensor(t);
    bulider.makeTensorPositive();
    bool linear = builder.isLinear();
    builder.computeAAndB();
    params.A = builder.getA();
    params.B = builder.getB();
    points.push_back(buider.evaluateSuperquadricPoint(linear, params));
  }

  return points[0] = points[1];
}

void ShowUncertaintyGlyphsImpl::computeOffsetSurface()
{
  const static double h = 0.000001;
  const static double hHalf = 0.5 * h;

  for (int f = 0; f < fieldSize_; ++f)
  {
    TensorGlyphBuilder builder(meanTensors_[f], points_[0][f]);
    builder.setResolution(resolution_);
    builder.makeTensorPositive();
    builder.computeTransforms();
    builder.postScaleTransorms();
    builder.computeSinCosTable(false);

    auto trans = builder.getTrans();
    bool linear = builder.isLinear();
    SuperquadricPointParams params;
    params.A = builder.getA();
    params.B = builder.getB();

    Point point;

    std::vector<DenseMatrix> diffs(6);
    auto finiteDiff = DenseMatrix(6, 1, 0.0);

    for (int v=0; v < nv-1; v++)
    {
      params.sinPhi = builder.computeSinPhi(v);
      params.cosPhi = builder.computeCosPhi(v);
      for (int u=0; u<nu; u++)
      {
        params.sinTheta = builder.computeSinTheta(u);
        params.cosTheta = builder.computeCosTheta(u);
        point = builder.evaluateSuperquadricPoint(linear, params);
        point = trans * point;

        auto mean = meanTensors_[f].mandel();
        for (int i = 0; i < 6; ++i)
        {
          finiteDiff.put(i, 0, hHalf);
          diffs[i] = diffT(builder, mean, finiteDiff, params) / h;
          finiteDiff.put(i, 0, 0.0);
        }
      }
    }

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

  return sum / fieldCount_;
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
AlgorithmOutput ShowUncertaintyGlyphsAlgorithm::run(const AlgorithmInput &input) const
{
  auto fields = input.getList<Field>(Variables::InputFields);
  if (fields.empty())
    THROW_ALGORITHM_INPUT_ERROR("No input fields given");

  auto impl = ShowUncertaintyGlyphsImpl();

  AlgorithmOutput output;
  impl.run(fields);
  output[MeanTensorField] = impl.getMeanTensors();
  return output;
}
