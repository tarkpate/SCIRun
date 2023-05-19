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
#include <Core/Datatypes/Legacy/Field/VMesh.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
// For mapping matrices
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>
#include <Core/Algorithms/Base/AlgorithmPreconditions.h>
#include <Core/Algorithms/Base/AlgorithmVariableNames.h>

#include <Core/Datatypes/Dyadic3DTensor.h>
#include <Core/Utils/TensorUtil.h>

//STL classes needed
#include <algorithm>
#include <set>
#include "Core/Datatypes/Legacy/Field/LatVolMesh.h"
#include "Core/Datatypes/Legacy/Field/PointCloudMesh.h"
#include "Core/Datatypes/TensorFwd.h"

using namespace SCIRun;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Algorithms::Fields;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Geometry;
using namespace SCIRun::Core::Logging;


ALGORITHM_PARAMETER_DEF(Fields, AddConstraints);
ALGORITHM_PARAMETER_DEF(Fields, CoalesceMethod);
// ALGORITHM_PARAMETER_DEF(Fields, IsoValue);

enum class FieldDataType
{
  node,
  edge,
  face,
  cell
};

CoalesceMeshAlgo::CoalesceMeshAlgo()
{
		using namespace Parameters;
		addOption(AddConstraints,"all","all|greaterthan|unequal|lessthan|none");
		addOption(CoalesceMethod,"Default","Default|Expand coalescement volume to improve element quality");
		// addParameter(IsoValue,0.0);
}

AlgorithmOutput CoalesceMeshAlgo::run(const AlgorithmInput& input) const
{
	auto field = input.get<Field>(Variables::InputField);
	auto isoVals = input.get<Field>(IsoValueField);
  FieldHandle outputField;

  if (!runImpl(field, isoVals, outputField))
    THROW_ALGORITHM_PROCESSING_ERROR("False returned on legacy run call.");
	AlgorithmOutput output;
	output[Variables::OutputField] = outputField;
  return output;
}

VMesh::Node::index_type getArrayIndex(VMesh::dimension_type dim, int i, int j, int k) {
  return k*(dim[1]*dim[0]) + j*(dim[0]) + i;
}

Dyadic3DTensor getTensor(VField* vfield, VMesh::Node::index_type idx) {
  Tensor t;
  vfield->get_value(t, idx);
  return scirunTensorToEigenTensor(t);
}

Point getPoint(VMesh* vmesh, VMesh::Node::index_type idx) {
  Point p;
  vmesh->get_center(p, idx);
  return p;
}

float getUncertaintyValue(VField* vfield, VMesh::Node::index_type idx) {
  float f;
  vfield->get_value(f, idx);
  return f;
}

void addTensorToField(FieldHandle& field, Point& point, Dyadic3DTensor& tensor) {
  auto vmesh = field->vmesh();
  auto vfield = field->vfield();
  auto new_idx = vmesh->add_node(point);

  VMesh::dimension_type dim;
  vmesh->get_elem_dimensions(dim);

  VMesh::Node::size_type sz;
  vmesh->size(sz);

  vfield->resize_fdata();
  vfield->size(sz);
  vfield->set_value(eigenTensorToScirunTensor(tensor), new_idx);
}

// General access function
bool
CoalesceMeshAlgo::runImpl(FieldHandle& input, FieldHandle& isoValueField, FieldHandle& output) const
{
	REPORT_STATUS(CoalesceMesh);
  if (!input) {
    error("No input field");
    return false;
  }
  if (!isoValueField) {
    error("No isoValue field");
    return false;
  }

	FieldInformation fi(input);
	FieldInformation fo(input);
  fo.make_pointcloudmesh();

	const std::string rMethod = getOption(Parameters::CoalesceMethod);
	// const double isoVals = get(Parameters::IsoValue).toDouble();
	const std::string addCon = getOption(Parameters::AddConstraints);

  if (input->vfield()->num_values() == 0)
  {
    error("Input field has no data values. The CoalesceMesh algorithm requires input fields to contain data.");
    return false;
  }

  if (addCon == "none")
  {
    output = input;
    return true;
  }

  if (fi.is_pnt_element() || fi.is_prism_element())
  {
    error("This algorithm does not support point or prism meshes");
    return false;
  }

  if (!fi.is_tensor())
  {
    error("Field data needs to be of tensor type");
    return false;
  }

  if (!input) {
    // Error reporting:
    // we forward the specific message to the ProgressReporter and return a
    // false to indicate that an error has occured.
    error("Could not obtain input field");
    return false;
  }

  output = CreateField(fo);

  VField* ifield = input->vfield();
  VField* ofield = output->vfield();
  VMesh*  imesh = input->vmesh();
  VMesh*  omesh = output->vmesh();
  VField* unc_field = isoValueField->vfield();

  VMesh::dimension_type dim;
  imesh->get_dimensions(dim);
  for (int i = 0; i < DIM_SIZE_; ++i)
    dim[i] = std::max(static_cast<int>(dim[i]), 1);

  for (int i = 0; i < dim[0]-1; ++i)
    for (int j = 0; j < dim[1]-1; ++j)
      for (int k = 0; k < dim[2]-1; ++k) {
        // Get indices for block of 8 points
        std::vector<VMesh::Node::index_type> indices = {
          getArrayIndex(dim, i, j, k),
          getArrayIndex(dim, i, j, k+1),
          getArrayIndex(dim, i, j+1, k),
          getArrayIndex(dim, i, j+1, k+1),
          getArrayIndex(dim, i+1, j, k),
          getArrayIndex(dim, i+1, j, k+1),
          getArrayIndex(dim, i+1, j+1, k),
          getArrayIndex(dim, i+1, j+1, k+1)};

        std::vector<Point> points(BLOCK_SIZE_);
        std::vector<Dyadic3DTensor> tensors(BLOCK_SIZE_);
        float unc_val_max;
        auto avgPoint = Point();
        auto avgTensor = Dyadic3DTensor();
        for (int p = 0; p < BLOCK_SIZE_; ++p) {
          points[p] = getPoint(imesh, indices[p]);
          avgPoint += points[p];
          tensors[p] = getTensor(ifield, indices[p]);
          avgTensor += tensors[p];
          unc_val_max = std::max(unc_val_max, getUncertaintyValue(unc_field, indices[p]));
        }
        auto unc_threshold = 0.5;
        if (unc_val_max > unc_threshold) {
          for (int p = 0; p < BLOCK_SIZE_; ++p) {
            addTensorToField(output, points[p], tensors[p]);
          }
        } else {
          avgPoint /= static_cast<float>(BLOCK_SIZE_);
          avgTensor = avgTensor / static_cast<double>(BLOCK_SIZE_);
          addTensorToField(output, avgPoint, avgTensor);
        }
      }

  return true;
}
