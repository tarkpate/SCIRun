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
#include <Core/Algorithms/Math/TensorInterpolation.h>

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
ALGORITHM_PARAMETER_DEF(Fields, IsoValue);
ALGORITHM_PARAMETER_DEF(Fields, CoalesceCount);
ALGORITHM_PARAMETER_DEF(Fields, NeighborThreshold);
ALGORITHM_PARAMETER_DEF(Fields, BlockSize);
ALGORITHM_PARAMETER_DEF(Fields, OverlapSize);

enum class FieldDataType
{
  node,
  edge,
  face,
  cell
};

class CoalesceDataNode {
  bool coalesced;
  bool valid;
  int neighbors_coalesced;
  Dyadic3DTensor t;
  Point p;
  float v;
  public:
    CoalesceDataNode() {
      coalesced = false;
      valid = false;
      neighbors_coalesced = 0;
      t = Dyadic3DTensor();
      p = Point();
      v = 0;
    }
    CoalesceDataNode(const CoalesceDataNode& other) {
      coalesced = other.isCoalesced();
      valid = other.isValid();
      neighbors_coalesced = other.getNeighborsCoalesced();
      t = other.getTensor();
      p = other.getPoint();
      v = other.getUncertaintyVal();
    }
    // oalesceDataNode(const Point& p_in, const Dyadic3DTensor& t_in, bool coalesced_in, float v_in) {
    //   p = p_in;
    //   t = t_in;
    //   v = v_in;
    //   coalesced_in = coalesced;
    // }
    void resetNeighbors() { neighbors_coalesced=0; }
    void setValid(bool val) { valid=val; }
    void setCoalesced(bool val) { coalesced=val; }
    void addNeighborCoalesced() { neighbors_coalesced++; }
    void setTensor(Dyadic3DTensor val) { t=val; }
    void setPoint(Point val) { p=val; }
    void setUncertaintyVal(float val) { v=val; }
    bool isCoalesced() const { return coalesced; }
    bool isValid() const { return valid; }
    int getNeighborsCoalesced() const { return neighbors_coalesced; }
    const Dyadic3DTensor& getTensor() const { return t; }
    const Point& getPoint() const { return p; }
    float getUncertaintyVal() const { return v; }

    CoalesceDataNode& operator+(const CoalesceDataNode& other) {
      auto otherT = other.getTensor();
      t += otherT;
      // p = p + Vector(other.getPoint());
      for (int i = 0; i < 3; ++i)
        p[i] += other.getPoint()[i];
      v += other.getUncertaintyVal();
      return *this;
    }

    CoalesceDataNode& operator/(double div) {
      t = t / div;
      v = v / div;
      for (int i = 0; i < 3; ++i)
        p[i] /= div;
      return *this;
    }

    CoalesceDataNode& operator*(double mult) {
      t = t * mult;
      return *this;
    }
};

void addTensorToField(FieldHandle& field, const Point& point, const Dyadic3DTensor& tensor) {
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

class CoalesceData {
  public:
    typedef VMesh::Node::index_type data_index_type[3];
  private:
    std::vector<CoalesceDataNode> data;
    data_index_type dims;
    data_index_type data_dims; // need to maintain to keep track of indexing
    VMesh::Node::index_type getArrayIndex(int x, int y, int z) const {
        return z*(data_dims[1]*data_dims[0]) + y*(data_dims[0]) + x;
    }
    VMesh::Node::index_type getArrayIndex(data_index_type idx) const {
        return idx[2]*(data_dims[1]*data_dims[0]) + idx[1]*(data_dims[0]) + idx[0];
      // return idx[2]*(dim[1]*dim[0]) + idx[1]*(dim[0]) + idx[0];

        // return idx[2]*(dims[1]*dims[0]) + idx[1]*(dims[0]) + idx[0];
        // return idx[1]*(dims[2]*dims[0]) + idx[2]*(dims[0]) + idx[0];
        // return idx[2]*(dims[0]*dims[1]) + idx[0]*(dims[1]) + idx[1];
        // return idx[0]*(dims[2]*dims[1]) + idx[2]*(dims[1]) + idx[1];
        // return idx[1]*(dims[0]*dims[2]) + idx[0]*(dims[2]) + idx[2];
        // return idx[2]*(dims[1]*dims[2]) + idx[1]*(dims[2]) + idx[2];
    }

  public:
    CoalesceData(const data_index_type& in_dims) {
      // dims.resize(3);
      // for (int i = 0; i < in_dims.size(); ++i)
      // {
      //   dims[i] = in_dims[i];
      // }
      for (int i = 0; i < 3; ++i) {
        dims[i] = in_dims[i];
        data_dims[i] = in_dims[i];
      }
      data = std::vector<CoalesceDataNode>(data_dims[0]*data_dims[1]*data_dims[2]);
      // data = std::vector<std::vector<std::vector<CoalesceDataNode>>>(dims[0],
      // std::vector<std::vector<CoalesceDataNode>>(dims[1], std::vector<CoalesceDataNode>(dims[2])));
    }

    CoalesceDataNode* getNode(data_index_type idx) { return &data[getArrayIndex(idx)]; }
    CoalesceDataNode* getNode(int x, int y, int z) { return &data[getArrayIndex(x, y, z)]; }
    const CoalesceDataNode* getNodeConst(int x, int y, int z) const { return &data[getArrayIndex(x, y, z)]; }
    void setNode(data_index_type idx, const CoalesceDataNode& node) { data[getArrayIndex(idx)] = CoalesceDataNode(node); }
    void setNode(int x, int y, int z, const CoalesceDataNode& node) { data[getArrayIndex(x, y, z)] = CoalesceDataNode(node); }

    void addNodeIfValid(std::vector<CoalesceDataNode*>& vec, int i, int j, int k) {
      auto n = getNode(i,j,k);
      if (n->isValid()) vec.push_back(n);
    }

    void resetCoalesced() {
      for (int k = 0; k < dims[2]; k++)
        for (int j = 0; j < dims[1]; j++)
          for (int i = 0; i < dims[0]; i++) {
            auto n = getNode(i,j,k);
            n->setCoalesced(false);
            // n->setNeighborCoalesced(false);
          }
    }

    void addUncoalesced(FieldHandle& output, int neighborThres) {
      bool flat = dims[0]==1 || dims[1]==1 || dims[2]==1;
      for (int k = 0; k < dims[2]; k++)
        for (int j = 0; j < dims[1]; j++)
          for (int i = 0; i < dims[0]; i++) {
            // for (int i = 0; i < data.size(); ++i) {
            // auto n = data[i];
            auto n = getNode(i,j,k);
            if (n->isValid() && n->getNeighborsCoalesced() < neighborThres) {
              addTensorToField(output, n->getPoint(), n->getTensor());
            }
            continue;
            // if (!n.isCoalesced() && n.isValid() && !n.isNeighborCoalesced()) {
            // if (!n.isNeighborCoalesced()) {
            if (n->isValid() && n->getNeighborsCoalesced() > neighborThres) {
              if (!flat && (i==(dims[0]-1) || j==(dims[1]-1) || k==(dims[2]-1))) {
                  addTensorToField(output, n->getPoint(), n->getTensor());
                  n->setValid(false);
              } else if (!n->isCoalesced()) {
                addTensorToField(output, n->getPoint(), n->getTensor());
                n->setValid(false);
              }
            }
            // if (n->isValid() && ((!n->isNeighborCoalesced() && (i==(dims[0]-1) || j==(dims[1]-1) || k==(dims[2]-1))))) {
              // addTensorToField(output, n->getPoint(), n->getTensor());
              // n->setValid(false);
            // }
            // n->setCoalesced(false);
            // n->setNeighborCoalesced(false);
          }
    }

    void decrementDims() {
      for (int i = 0; i < 3; ++i)
        dims[i] = std::max(static_cast<int>(dims[i]-1), 1);
    }

    void makeNode(index_type idx, const Point& p, const Dyadic3DTensor& t, bool coalesced, bool valid, float unc) {
      // auto n = getDataNode(idx);
      // data[idx[0]][idx[1]][idx[2]] = CoalesceDataNode();
      // CoalesceDataNode& n = getDataNode(idx);
      CoalesceDataNode& n = data[idx];
      // n = CoalesceDataNode(p, t, coalesced, unc);
      n.setPoint(p);
      // std::cout << "make p: " << getDataNode(idx).getPoint() << "\n";
      n.setTensor(t);
      n.setValid(valid);
      n.setCoalesced(coalesced);
      n.setUncertaintyVal(unc);
    }

    const data_index_type& getDims() const { return dims; }
    // void setCoalesced(data_index_type idx, bool val) { getNode(idx)->setCoalesced(val); }
    // void setTensor(data_index_type idx, Dyadic3DTensor val) { getNode(idx)->setTensor(val); }
    // void setPoint(data_index_type idx, Point val) { getNode(idx)->setPoint(val); }
    // void setUncertaintyVal(data_index_type idx, float val) { getNode(idx)->setUncertaintyVal(val); }
    // bool getCoalesced(data_index_type idx) { return getNode(idx)->isCoalesced(); }
    // const Dyadic3DTensor& getTensor(data_index_type idx) { return getNode(idx)->getTensor(); }
    const Dyadic3DTensor& getTensor(int x, int y, int z) { return getNode(x, y, z)->getTensor(); }
    // const Point& getPoint(data_index_type idx) { return getNode(idx)->getPoint(); }
    const Point& getPoint(int x, int y, int z) { return getNode(x, y, z)->getPoint(); }
    // float getUncertaintyVal(data_index_type idx) { return getNode(idx)->getUncertaintyVal(); }
};

// VMesh::Node::index_type getArrayIndex(CoalesceData::data_index_type dim, int i, int j, int k) {
  // return k*(dim[1]*dim[0]) + j*(dim[0]) + i;
// }

CoalesceMeshAlgo::CoalesceMeshAlgo()
{
		using namespace Parameters;
		addOption(AddConstraints,"all","all|greaterthan|unequal|lessthan|none");
		addOption(CoalesceMethod,"Default","Default|Expand coalescement volume to improve element quality");
		addParameter(IsoValue,0.5);
		addParameter(CoalesceCount,1);
		addParameter(NeighborThreshold,2);
		addParameter(BlockSize,3);
		addParameter(OverlapSize,1);
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
  return std::abs(f);
}

void addRemaining(CoalesceData& coalesceData, FieldHandle& output) {
  const CoalesceData::data_index_type& dims = coalesceData.getDims();
  for (int k = 0; k < dims[2]; k++)
    for (int j = 0; j < dims[1]; j++)
      for (int i = 0; i < dims[0]; i++) {
        auto n = coalesceData.getNode(i,j,k);
        if (n->isValid())
          addTensorToField(output, n->getPoint(), n->getTensor());
      }
}
void furtherCoalesce(CoalesceData& coalesceData, CoalesceData& newData, float unc_threshold, int blockSize, int overlapSize, int borderSize) {

  TensorInterpolationAlgo interpAlgo(TensorInterpolationAlgo::Method::LINEAR_INVARIANT,
                                     TensorInterpolationAlgo::Method::LOG_EUCLIDEAN);

  const CoalesceData::data_index_type& dims = coalesceData.getDims();
  const CoalesceData::data_index_type& newDims = newData.getDims();
  bool x_flat = dims[0] <= 1;
  bool y_flat = dims[1] <= 1;
  bool z_flat = dims[2] <= 1;
  int incSize = blockSize - overlapSize;
  for (int k = 0; k < newDims[2]; k++) // TODO change increment
    for (int j = 0; j < newDims[1]; j++)
      for (int i = 0; i < newDims[0]; i++) {
        int sampleX = x_flat ? 0 : i*incSize + borderSize;
        int sampleY = y_flat ? 0 : j*incSize + borderSize;
        int sampleZ = z_flat ? 0 : k*incSize + borderSize;

        std::vector<CoalesceDataNode*> nodes = {coalesceData.getNode(sampleX,sampleY,sampleZ)};
        if (!nodes[0]->isValid()) continue;
        for (int d = -blockSize/2; d <= blockSize/2; ++d) {
          if (!x_flat)
            nodes.push_back(coalesceData.getNode(sampleX+d,sampleY,sampleZ));
          if (!y_flat)
            nodes.push_back(coalesceData.getNode(sampleX,sampleY+d,sampleZ));
          if (!z_flat)
            nodes.push_back(coalesceData.getNode(sampleX,sampleY,sampleZ+d));
          if (!x_flat && !y_flat)
            nodes.push_back(coalesceData.getNode(sampleX+d,sampleY+d,sampleZ));
          if (!x_flat && !z_flat)
            nodes.push_back(coalesceData.getNode(sampleX+d,sampleY,sampleZ+d));
          if (!y_flat && !z_flat)
            nodes.push_back(coalesceData.getNode(sampleX,sampleY+d,sampleZ+d));
          if (!x_flat && !y_flat && !z_flat)
            nodes.push_back(coalesceData.getNode(sampleX+d,sampleY+d,sampleZ+d));
        }
        bool valid = true;
        for (auto n : nodes)
          if (!n->isValid()) valid = false;
        if (!valid) continue;

        std::vector<Dyadic3DTensor> tensors(nodes.size());
        Point avgPoint = Point(0, 0, 0);
        float maxUncVal = 0;
        for (int n = 0; n < nodes.size(); ++n) {
          tensors[n] = nodes[n]->getTensor();
          avgPoint += nodes[n]->getPoint();
          maxUncVal = std::max(maxUncVal, nodes[n]->getUncertaintyVal());
        }
        avgPoint /= nodes.size();

        const static double scale_factor = 2.0;
        Dyadic3DTensor avgTensor = interpAlgo.interpolate(tensors);
        avgTensor = avgTensor * scale_factor;

        // CoalesceDataNode avg(*nodes[0]);
        // for (int n = 1; n < nodes.size(); ++n) {
          // if (!nodes[n]->isValid()) continue;
          // avg = avg + *nodes[n];
        // }
        CoalesceDataNode avgNode;
        avgNode.setTensor(avgTensor);
        avgNode.setUncertaintyVal(maxUncVal);
        if (avgNode.getUncertaintyVal() < unc_threshold) {
          avgNode.setValid(true);
          avgNode.resetNeighbors();
          avgNode.setPoint(avgPoint);
          newData.setNode(i,j,k, avgNode);
          coalesceData.getNode(sampleX,sampleY,sampleZ)->setValid(false);
          for (auto n : nodes) n->addNeighborCoalesced();
        }
    }
  // std::cout << "actual x dim: " << count << "\n";
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
	const double uncThreshold = get(Parameters::IsoValue).toDouble();
	const int blockSize = get(Parameters::BlockSize).toInt();
	const int overlapSize = get(Parameters::OverlapSize).toInt();
	const int coalesceCount = get(Parameters::CoalesceCount).toInt();
	const std::string addCon = getOption(Parameters::AddConstraints);
	const int neighborThres = get(Parameters::NeighborThreshold).toInt();

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

  VMesh::dimension_type mesh_dims;
  imesh->get_dimensions(mesh_dims);
  CoalesceData::data_index_type dims;
  for (int i = 0; i < DIM_SIZE_; ++i)
    dims[i] = std::max(static_cast<int>(mesh_dims[i]), 1);

  std::vector<VMesh::Node::index_type> indices;
  if (imesh->is_imagemesh()) {
    dims[2] = 1;
  }

  CoalesceData coalesceData(dims);
  CoalesceData::data_index_type idx;
  for (int i = 0; i < dims[0]*dims[1]*dims[2]; ++i)
        coalesceData.makeNode(i, getPoint(imesh, i),
                              getTensor(ifield, i),
                              false,
                              true,
                              getUncertaintyValue(unc_field, i));
  bool x_flat = dims[0] <= 1;
  bool y_flat = dims[1] <= 1;
  bool z_flat = dims[2] <= 1;
  VMesh::index_type borderSize = blockSize / 2; // Remainder ignored
  printf("borderSize %i\n", borderSize);
  // int x_start = x_flat ? 0 : borderSize;
  // int y_start = y_flat ? 0 : borderSize;
  // int z_start = z_flat ? 0 : borderSize;
  CoalesceData::data_index_type newDims;

  for (int c = 0; c < coalesceCount; ++c) {
    // int x_end = x_flat ? 1 : static_cast<int>(dims[0])-borderSize;
    // int y_end = y_flat ? 1 : static_cast<int>(dims[1])-borderSize;
    // int z_end = z_flat ? 1 : static_cast<int>(dims[2])-borderSize;
    int incSize = blockSize - overlapSize;

    for (int i = 0; i < 3; ++i)
      newDims[i] = std::max(1, static_cast<int>(VMesh::index_type(dims[i] - 2*borderSize) / incSize + 1));
    // std::cout << "new dims: " << newDims[0] << ", " << newDims[1] << ", " << newDims[2] << "\n";
    CoalesceData newData(newDims);

    furtherCoalesce(coalesceData, newData, uncThreshold, blockSize, overlapSize, borderSize);
    coalesceData.addUncoalesced(output, neighborThres);
    coalesceData = newData;
    for (int i = 0; i < 3; ++i)
      dims[i] = newDims[i];
  }

  addRemaining(coalesceData, output);

  return true;
}
