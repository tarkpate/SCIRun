#include <Core/Algorithms/Visualization/ShowUncertaintyGlyphsAlgorithm.h>
#include <Core/Algorithms/Base/AlgorithmVariableNames.h>
#include <Core/Algorithms/Base/AlgorithmPreconditions.h>
#include <Core/Logging/Log.h>
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/Mesh/MeshFacade.h>

using namespace SCIRun;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Geometry;
using namespace SCIRun::Core::Algorithms::Visualization;

//  this function is for setting defaults for state variables.
//  Mostly for UI variables.
ShowUncertaintyGlyphsAlgorithm::ShowUncertaintyGlyphsAlgorithm()
{
  // using namespace Parameters;
  // addParameter(Knob1, false);
  // addParameter(Knob2, 1.0);
}

//TODO move to separate algo
enum FieldDataType
{
  node,
  edge,
  face,
  cell
};

void ShowUncertaintyGlyphsAlgorithm::getPointsForFields(FieldHandle field,
                                                        std::vector<int> &indices,
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

void ShowUncertaintyGlyphsAlgorithm::verifyData(const FieldList &fields)
{
  for(int f = 1; f < fieldCount_; ++f)
    if(indices_[f].size() != fieldSize_)
      throw std::invalid_argument("All field inputs must have the same size.");

  // Verify all are tensors
  for(auto field : fields)
  {
    FieldInformation finfo(field);
    if(!finfo.is_tensor())
      throw std::invalid_argument("Currently, this module only supports tensor fields.");
  }
}

void ShowUncertaintyGlyphsAlgorithm::getPoints(const FieldList &fields)
{
  fieldCount_ = fields.size();
  fieldSize_ = indices_[0].size();

  indices_ = std::vector<std::vector<int>>(fieldCount_);
  points_  = std::vector<std::vector<Point>>(fieldCount_);

  for(int f = 0; f < fieldCount_; ++f)
    getPointsForFields(fields[f], indices_[f], points_[f]);
}

Tensor ShowUncertaintyGlyphsAlgorithm::computeMeanTensor(int index) const
{
  Tensor sum = Tensor();
  for(const auto tensor : tensors_[index])
    sum += tensor;

  // double div = (1.0 / fields.size());
  Tensor mean = sum / fieldCount_;
  return mean;
}

void ShowUncertaintyGlyphsAlgorithm::getTensors(FieldList &fields)
{
  tensors_ = std::vector<std::vector<Tensor>>(fieldCount_);
  for(int f = 0; f < fieldCount_; ++f)
  {
    tensors_[f] = std::vector<Tensor>(fieldSize_);
    auto vfield = fields[f]->vfield();
    for (int v = 0; v < fieldSize_; ++v)
      vfield->get_value(tensors_[f][v], v);
  }
}

void ShowUncertaintyGlyphsAlgorithm::computeMeanTensors()
{
  meanTensors_ = std::vector<Tensor>(fieldCount_);
  for (int f = 0; f < fieldCount_; ++f)
    meanTensors_[f] = computeMeanTensor(f);
}

FieldHandle ShowUncertaintyGlyphsAlgorithm::createOutputField() const
{
  FieldInformation ofinfo("PointCloudMesh", 0, "Tensor");
  auto ofield = CreateField(ofinfo);
  auto mesh = ofield->vmesh();
  auto field = ofield->vfield();

  std::vector<VMesh::index_type> meshIndices(fieldSize_);
  for(int i = 0; i < fieldSize_; ++i)
    meshIndices[i] = mesh->add_point(points_[0][i]);

  field->resize_fdata();

  for(int i = 0; i < fieldSize_; ++i)
    field->set_value(meanTensors_[i], meshIndices[i]);

  return ofield;
}

//main algorithm function
AlgorithmOutput ShowUncertaintyGlyphsAlgorithm::run(const AlgorithmInput &input)
{
  auto fields = input.getList<Field>(Variables::InputFields);
  if (fields.empty())
    THROW_ALGORITHM_INPUT_ERROR("No input fields given");

  getPoints(fields);
  verifyData(fields);
  getTensors(fields);
  computeMeanTensors();

  AlgorithmOutput output;
  output[Variables::InputFields] = createOutputField();
  return output;
}
