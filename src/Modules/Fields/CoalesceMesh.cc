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


#include <Core/Algorithms/Field/CoalesceMeshAlgorithm.h>
#include <Core/Datatypes/Legacy/Field/Field.h>
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <Core/Datatypes/Mesh/MeshFacade.h>
#include <Core/Datatypes/Scalar.h>
#include <Core/GeometryPrimitives/Tensor.h>
#include <Modules/Fields/CoalesceMesh.h>
#include <Dataflow/Network/ModuleStateInterface.h>
#include <Core/Datatypes/Legacy/Base/PropertyManager.h>


using namespace SCIRun::Modules::Fields;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Algorithms::Fields;
using namespace SCIRun::Core::Geometry;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Dataflow::Networks;
using namespace SCIRun;

MODULE_INFO_DEF(CoalesceMesh, ChangeFieldData, SCIRun)

CoalesceMesh::CoalesceMesh() : Module(staticInfo_)
{
  INITIALIZE_PORT(InputField);
  INITIALIZE_PORT(OutputField);
  INITIALIZE_PORT(IsoValueField);
}

void CoalesceMesh::setStateDefaults()
{
  setStateStringFromAlgoOption(Parameters::AddConstraints);
  setStateStringFromAlgoOption(Parameters::CoalesceMethod);
  setStateDoubleFromAlgo(Parameters::IsoValue);
}

// TODO move to separate algo
enum FieldDataType
{ node, edge, face, cell };

void get_dimensions(FieldHandle fh, std::vector<Mesh::index_type>& dims) {
  VMesh *mesh = fh->vmesh();
  mesh->get_dimensions(dims);
}

void getPointsForFields(
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

void CoalesceMesh::execute()
{
  // FieldList input = {getRequiredInput(InputField), getRequiredInput(IsoValueField)};
  auto input = getRequiredInput(InputField);
  // auto isovaluefield = getOptionalInput(IsoValueField);

  if (needToExecute()) {
    setAlgoOptionFromState(Parameters::AddConstraints);
    setAlgoOptionFromState(Parameters::CoalesceMethod);
    setAlgoDoubleFromState(Parameters::IsoValue);

    std::vector<Mesh::index_type> dims;
    get_dimensions(input, dims);

    auto vfield = input->vfield();
    Tensor temp;
    auto indices = std::vector<int>();
    auto points = std::vector<Point>();
    getPointsForFields(input, indices, points);
    auto fieldSize = indices.size();
    auto tensors = std::vector<Tensor>(fieldSize);
    for (size_t v = 0; v < fieldSize; ++v)
    {
      vfield->get_value(temp, v);
      tensors[v] = temp;
    }
    std::cout << "dims: " << dims[0] << ", " << dims[1] << ", " << dims[2] << std::endl;



    // auto output = algo().run(withInputData((InputField, input)));
    // sendOutputFromAlgorithm(OutputField, output);
  }
}
