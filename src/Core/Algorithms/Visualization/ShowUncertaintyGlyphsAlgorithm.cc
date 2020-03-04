#include <Core/Algorithms/Visualization/ShowUncertaintyGlyphsAlgorithm.h>
#include <Core/Algorithms/Base/AlgorithmVariableNames.h>
#include <Core/Algorithms/Base/AlgorithmPreconditions.h>
#include <Core/Logging/Log.h>
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/Mesh/MeshFacade.h>
#include <Core/Datatypes/DenseMatrix.h>

using namespace SCIRun;
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

  int fieldCount_ = 0;
  int fieldSize_ = 0;

  std::vector<std::vector<Point>> points_;
  std::vector<std::vector<int>> indices_;
  std::vector<std::vector<Tensor>> tensors_;
  std::vector<Tensor> meanTensors_;
  std::vector<DenseMatrix> covarianceMatrices_;
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
      auto diffTensorMatrix = DenseMatrix(6, 1);
      for(int i = 0; i < 6; ++i)
        diffTensorMatrix.put(i, 0, diffTensor[i]);
      cov += diffTensorMatrix * diffTensorMatrix.transpose();
    }

    cov /= fieldCount_;
    covarianceMatrices_[t] = cov;
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
