#include <Modules/Visualization/ShowUncertaintyGlyphs.h>
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/Legacy/Field/Mesh.h>
#include <Core/Datatypes/Mesh/MeshFacade.h>
#include <Core/Datatypes/Legacy/Field/Field.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>
#include <Core/GeometryPrimitives/Tensor.h>

using namespace SCIRun;
using namespace Core::Datatypes;
using namespace Dataflow::Networks;
using namespace Core::Geometry;
using namespace Modules::Visualization;

MODULE_INFO_DEF(ShowUncertaintyGlyphs, Visualization, SCIRun);

namespace SCIRun {
namespace Modules {
namespace Visualization {

class GlyphBuilder
{
public:
  GlyphBuilder(const std::string& moduleId) : moduleId_(moduleId){}
  void getPoints(FieldHandle field, std::vector<int>& indices, std::vector<Point>& points);
  Tensor getMeanTensor(std::vector<FieldHandle> fields, int index);
private:
  std::string moduleId_;
};

enum FieldDataType
{
  node,
  edge,
  face,
  cell
};

void GlyphBuilder::getPoints(FieldHandle field, std::vector<int>& indices, std::vector<Point>& points)
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

Tensor GlyphBuilder::getMeanTensor(std::vector<FieldHandle> fields, int index)
{
  std::vector<Tensor> tensors(fields.size());
  for(int i = 0; i < fields.size(); ++i)
  {
    auto vfield = fields[i]->vfield();
    vfield->get_value(tensors[i], index);
  }
  Tensor sum = Tensor();
  for(auto tensor : tensors)
    sum = sum + tensor;

  double div = (1.0 / fields.size());
  Tensor mean = sum * div;
  return mean;
}


ShowUncertaintyGlyphs::ShowUncertaintyGlyphs() : GeometryGeneratingModule(staticInfo_),
                                                 builder_(new GlyphBuilder(id().id_))
{
  INITIALIZE_PORT(InputFields);
  INITIALIZE_PORT(Mean);
  INITIALIZE_PORT(OutputGeom);
}

void ShowUncertaintyGlyphs::setStateDefaults()
{
  auto state = get_state();
}

void ShowUncertaintyGlyphs::execute()
{
  auto fields = getRequiredDynamicInputs(InputFields);

  if(needToExecute())
  {
    auto state = get_state();

    // Get data
    int fieldCount = fields.size();
    auto indices = std::vector<std::vector<int>>(fieldCount);
    auto points = std::vector<std::vector<Point>>(fieldCount);
    for (int i = 0; i < fieldCount; ++i)
      builder_->getPoints(fields[i], indices[i], points[i]);

    // Verify size is same
    int requiredSize = indices[0].size();
    for(int i = 1; i < fieldCount; ++i)
      if(indices[i].size() != requiredSize)
        throw std::invalid_argument("All field inputs must have the same size.");

    // Verify all are tensors
    for(auto field : fields)
    {
      FieldInformation finfo(field);
      if(!finfo.is_tensor())
        throw std::invalid_argument("This module only supports tensor fields.");
    }

    std::vector<Tensor> meanTensors(requiredSize);
    for (int i = 0; i < requiredSize; ++i)
      meanTensors[i] = builder_->getMeanTensor(fields, i);

    // Create output field
    FieldInformation firstFieldInfo(fields[0]);
    FieldInformation ofinfo("PointCloudMesh", 0, "Tensor");
    auto ofield = CreateField(ofinfo);
    auto mesh = ofield->vmesh();
    auto field = ofield->vfield();

    std::vector<VMesh::index_type> meshIndices(requiredSize);
    // Assumes all points are at same location
    for(int i = 0; i < requiredSize; ++i)
      meshIndices[i] = mesh->add_point(points[0][i]);
    field->resize_fdata();
    for(int i = 0; i < requiredSize; ++i)
      field->set_value(meanTensors[i], meshIndices[i]);

    sendOutput(Mean, ofield);
  }
}

}}}
