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
ALGORITHM_PARAMETER_DEF(Fields, IsoValue);
ALGORITHM_PARAMETER_DEF(Fields, CoalesceCount);
ALGORITHM_PARAMETER_DEF(Fields, NeighborThreshold);

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
            if (n->getNeighborsCoalesced() < neighborThres) {
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
void furtherCoalesce(CoalesceData& coalesceData, CoalesceData& newData, float unc_threshold) {
  // coalesceData.resetCoalesced();
  // coalesceData.decrementDims();
  const CoalesceData::data_index_type& dims = coalesceData.getDims();
  // CoalesceData::data_index_type idx;
  std::cout << "dims " << dims[0] << ", "<< dims[1] << ", "<< dims[2] << "\n";
  // for (idx[0] = 0; idx[0] < dims[0]; idx[0]++)
  //   for (idx[1] = 0; idx[1] < dims[1]; idx[1]++)
  //     for (idx[2] = 0; idx[2] < dims[2]; idx[2]++) {
  // addTensorToField(output, coalesceData.getPoint(0,0,0), coalesceData.getTensor(0,0,0));
  // addTensorToField(output, coalesceData.getPoint(1,0,0), coalesceData.getTensor(1,0,0));
  // addTensorToField(output, coalesceData.getPoint(0,1,0), coalesceData.getTensor(0,1,0));
  // addTensorToField(output, coalesceData.getPoint(0,0,1), coalesceData.getTensor(0,0,1));
  // return;
  bool x_flat = dims[0] <= 1;
  bool y_flat = dims[1] <= 1;
  bool z_flat = dims[2] <= 1;
  for (int k = 0; k < std::max(1, static_cast<int>(dims[2]-1)); k++)
    for (int j = 0; j < std::max(1, static_cast<int>(dims[1]-1)); j++)
      for (int i = 0; i < std::max(1, static_cast<int>(dims[0]-1)); i++) {
        // idx[0] = i;
        // idx[1] = j;
        // idx[2] = k;
        // std::cout << "add idx: " <<  idx[0] << ", "<<  idx[1] << ", "<<  idx[2] << "\n";
        // auto p = coalesceData.getPoint(idx);
        // auto t = coalesceData.getTensor(idx);
        // std::cout << "p: " <<  p << "\n";
        // std::cout << "t: " <<  t << "\n";
        std::vector<const CoalesceDataNode*> nodes = {coalesceData.getNodeConst(i,j,k)};
        // nodes[0].setCoalesced(false);
        if (!nodes[0]->isValid()) continue;
        if (!x_flat)
          nodes.push_back(coalesceData.getNodeConst(i+1,j,k));
        if (!y_flat)
          nodes.push_back(coalesceData.getNodeConst(i,j+1,k));
        if (!z_flat)
          nodes.push_back(coalesceData.getNodeConst(i,j,k+1));
        if (!x_flat && !y_flat)
          nodes.push_back(coalesceData.getNodeConst(i+1,j+1,k));
        if (!x_flat && !z_flat)
          nodes.push_back(coalesceData.getNodeConst(i+1,j,k+1));
        if (!y_flat && !z_flat)
          nodes.push_back(coalesceData.getNodeConst(i,j+1,k+1));
        if (!x_flat && !y_flat && !z_flat)
          nodes.push_back(coalesceData.getNodeConst(i+1,j+1,k+1));
        bool valid = true;
        for (auto n : nodes)
          if (!n->isValid()) valid = false;
        if (!valid) continue;
        // std::cout << "nodes size: "<< nodes.size() << "\n";
        // std::cout << "idx: " <<  i << ", "<<  j << ", "<<  k << "\n";
        // if (!x_flat)
        //   coalesceData.addNodeIfValid(nodes, i+1,j,k);
        // if (!y_flat)
        //   coalesceData.addNodeIfValid(nodes, i,j+1,k);
        // if (!z_flat)
        //   coalesceData.addNodeIfValid(nodes, i,j,k+1);
        // if (!x_flat && !y_flat)
        //   coalesceData.addNodeIfValid(nodes, i+1,j+1,k);
        // if (!x_flat && !z_flat)
        //   coalesceData.addNodeIfValid(nodes, i+1,j,k+1);
        // if (!y_flat && !z_flat)
        //   coalesceData.addNodeIfValid(nodes, i,j+1,k+1);
        // if (!x_flat && !y_flat && !z_flat)
        //   coalesceData.addNodeIfValid(nodes, i+1,j+1,k+1);

        // for (int n = 0; n < nodes.size(); ++n)
            // nodes[n]->setNeighborCoalesced(true);

        CoalesceDataNode avg(*nodes[0]);
        for (int n = 1; n < nodes.size(); ++n) {
          if (!nodes[n]->isValid()) continue;
          avg = avg + *nodes[n];
        }
        double scale_factor = 1.2;
        avg = avg / nodes.size() * scale_factor;
        if (avg.getUncertaintyVal() < unc_threshold) {
          // avg.setCoalesced(true);
          avg.setValid(true);
          avg.resetNeighbors();
          newData.setNode(i,j,k, avg);
          // printf("coalesced\n");
          // for (int n = 1; n < nodes.size(); ++n)
            // nodes[n]->setNeighborCoalesced(true);
          // newData.getNode(i,j,k)->setCoalesced(true);
          if (!x_flat) {
            coalesceData.getNode(i+1,j,k)->addNeighborCoalesced();
            if (i > 0)
              coalesceData.getNode(i-1,j,k)->addNeighborCoalesced();
          }
          if (!y_flat) {
            coalesceData.getNode(i,j+1,k)->addNeighborCoalesced();
            if (j > 0)
              coalesceData.getNode(i,j-1,k)->addNeighborCoalesced();
          }
          if (!z_flat) {
            coalesceData.getNode(i,j,k+1)->addNeighborCoalesced();
            if (k > 0)
              coalesceData.getNode(i,j,k-1)->addNeighborCoalesced();
          }
          if (!x_flat || !z_flat) {
            coalesceData.getNode(i+1,j,k+1)->addNeighborCoalesced();
            if (i > 0 && k > 0)
              coalesceData.getNode(i-1,j,k-1)->addNeighborCoalesced();
          }
          if (!x_flat || !y_flat) {
            coalesceData.getNode(i+1,j+1,k)->addNeighborCoalesced();
            if (i > 0 && j > 0)
              coalesceData.getNode(i-1,j-1,k)->addNeighborCoalesced();
          }
          if (!y_flat || !z_flat) {
            coalesceData.getNode(i,j+1,k+1)->addNeighborCoalesced();
            if (j > 0 && k > 0)
              coalesceData.getNode(i,j-1,k-1)->addNeighborCoalesced();
          }
          if (!y_flat || !z_flat) {
            coalesceData.getNode(i+1,j+1,k+1)->addNeighborCoalesced();
            if (i > 0 && j > 0 && k > 0)
              coalesceData.getNode(i-1,j-1,k-1)->addNeighborCoalesced();
          }
          // if (x_flat || y_flat || z_flat)
            // for (auto n : nodes) n.setCoalesced(true);
          // addTensorToField(output, avg.getPoint(), avg.getTensor());
        // } else {
          // for (int n = 0; n < nodes.size(); ++n) {
            // nodes[n].setCoalesced(false);
            // if (nodes.size() >= 8)
              // addTensorToField(output, nodes[n].getPoint(), nodes[n].getTensor());
          // }
        }
        // addTensorToField(output, coalesceData.getPoint(idx), Dyadic3DTensor(idx[2]*(dims[1]*dims[0]) + idx[1]*(dims[0]) + idx[0]));
        // addTensorToField(output, coalesceData.getPoint(idx), Dyadic3DTensor(idx[1]*(dims[2]*dims[0]) + idx[2]*(dims[0]) + idx[0]));
        // addTensorToField(output, coalesceData.getPoint(idx), Dyadic3DTensor(idx[2]*(dims[0]*dims[1]) + idx[0]*(dims[1]) + idx[1]));
        // addTensorToField(output, coalesceData.getPoint(idx), Dyadic3DTensor(idx[0]*(dims[2]*dims[1]) + idx[2]*(dims[1]) + idx[1]));
        // addTensorToField(output, coalesceData.getPoint(idx), Dyadic3DTensor(idx[1]*(dims[0]*dims[2]) + idx[0]*(dims[2]) + idx[2]));
        // addTensorToField(output, coalesceData.getPoint(idx), Dyadic3DTensor(idx[2]*(dims[1]*dims[2]) + idx[1]*(dims[2]) + idx[2]));
      }
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
	const double unc_threshold = get(Parameters::IsoValue).toDouble();
	const int coalesce_count = get(Parameters::CoalesceCount).toInt();
	const std::string addCon = getOption(Parameters::AddConstraints);
	const int neighborThres = get(Parameters::NeighborThreshold).toInt();
  std::cout << "neigh thres: " << neighborThres << "\n";

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

  // std::vector<std::vector<std::vector<CoalesceData>>> coalesceData(dim[0]-1,
  //   std::vector<std::vector<CoalesceData>>(dim[1]-1, std::vector<CoalesceData>(std::max(static_cast<int>(dim[2]-1),1))));
  CoalesceData coalesceData(dims);
  CoalesceData::data_index_type idx;
  // for (idx[2] = 0; idx[2] < dims[2]; idx[2]++)
  //   for (idx[1] = 0; idx[1] < dims[1]; idx[1]++)
  //     for (idx[0] = 0; idx[0] < dims[0]; idx[0]++) {
  for (int i = 0; i < dims[0]*dims[1]*dims[2]; ++i) {
        // auto array_index = getArrayIndex(dims, idx[0], idx[1], idx[2]);
        // std::cout << "make p: " << getPoint(imesh, array_index) << "\n";
        coalesceData.makeNode(i, getPoint(imesh, i),
                              getTensor(ifield, i),
                              false,
                              true,
                              getUncertaintyValue(unc_field, i));
        // std::cout << "make p: " << coalesceData.getPoint(idx) << "\n";
        // coalesceData.setPoint(idx, getPoint(imesh, array_index));
        // coalesceData.setTensor(idx, getTensor(ifield, array_index));
        // coalesceData.setCoalesced(idx, false);
        // coalesceData.setUncertaintyVal(idx, getUncertaintyValue(unc_field, array_index));
      }
  printf("1\n");


  // for (int i = 0; i < dims[0]; i++)
    // for (int j = 0; j < dims[1]; j++)
      // for (int k = 0; k < std::max(static_cast<int>(dims[2]-1), 1); k++) {
        // VMesh::Node::index_type index;
        // if (imesh->is_imagemesh())
          // index = getArrayIndex(dims, i, j, 0);
        // else if (imesh->is_latvolmesh())
          // index = getArrayIndex(dims, i, j, k);
        // coalesceData.setPoint(i, j, k, getPoint(imesh, indices[p]));
        // coalesceData.setTensor(i, j, k, getTensor(imesh, indices[p]));
        // coalesceData.setCoalsced(i, j, k, false);
        // coalesceData.setUncertaintyValue(i, j, k, getUncertaintyValue(unc_field, indices[p]));
  // }

  // std::vector<std::vector<std::vector<bool>>> coalesced(dim[0],
    // std::vector<std::vector<bool>>(dim[1], std::vector<bool>(std::max(static_cast<int>(dim[2]),1))));
  // std::vector<std::vector<std::vector<bool>>> coalesced(dim[0]-1,
    // std::vector<std::vector<bool>>(dim[1]-1, std::vector<bool>(std::max(static_cast<int>(dim[2]-1),1))));
  // std::vector<std::vector<std::vector<Dyadic3DTensor>>> coalescedTensors(dim[0],
    // std::vector<std::vector<Dyadic3DTensor>>(dim[1], std::vector<Dyadic3DTensor>(std::max(static_cast<int>(dim[2]-1), 1))));
  // std::cout << "new dim size " << coalescedTensors.size() << ", " << coalescedTensors[0].size() << ", " << coalescedTensors[0][0].size() << "\n";
  // std::vector<std::vector<std::vector<Point>>> coalescedPoints(dim[0],
    // std::vector<std::vector<Point>>(dim[1], std::vector<Point>(std::max(static_cast<int>(dim[2]-1), 1))));
  // std::vector<std::vector<std::vector<float>>> coalescedUncertaintyVals(dim[0],
    // std::vector<std::vector<float>>(dim[1], std::vector<float>(std::max(static_cast<int>(dim[2]-1), 1))));
  // for (int i = 0; i < dims[0]-1; i++)
  //   for (int j = 0; j < dims[1]-1; j++)
  //     for (int k = 0; k < std::max(static_cast<int>(dims[2]-1), 1); k++) {
  //       // Make list of indices for group
  //       if (imesh->is_imagemesh()) {
  //         indices = {
  //         getArrayIndex(dims, i, j, 0),
  //         getArrayIndex(dims, i, j+1, 0),
  //         getArrayIndex(dims, i+1, j, 0),
  //         getArrayIndex(dims, i+1, j+1, 0)};
  //       }
  //       if (imesh->is_latvolmesh()) {
  //         indices = {
  //           getArrayIndex(dims, i, j, k),
  //           getArrayIndex(dims, i, j, k+1),
  //           getArrayIndex(dims, i, j+1, k),
  //           getArrayIndex(dims, i, j+1, k+1),
  //           getArrayIndex(dims, i+1, j, k),
  //           getArrayIndex(dims, i+1, j, k+1),
  //           getArrayIndex(dims, i+1, j+1, k),
  //           getArrayIndex(dims, i+1, j+1, k+1)};
  //       }
  //       std::vector<Point> points(indices.size());
  //       std::vector<Dyadic3DTensor> tensors(indices.size());
  //       auto avgPoint = Point();
  //       auto avgTensor = Dyadic3DTensor();
  //       float unc_val_max = 0;
  //       // Average the points together
  //       for (int p = 0; p < indices.size(); ++p) {
  //         points[p] = getPoint(imesh, indices[p]);
  //         avgPoint += points[p];
  //         tensors[p] = getTensor(ifield, indices[p]);
  //         avgTensor += tensors[p];
  //         unc_val_max = std::max(unc_val_max, getUncertaintyValue(unc_field, indices[p]));
  //       }
  //       // If threshold is met, coalesce tensors
  //       if (unc_val_max > unc_threshold) {
  //         coalesceData[i][j][k].setCoalesced(false);
  //         for (int p = 0; p < indices.size(); ++p) {
  //           addTensorToField(output, points[p], tensors[p]);
  //         }
  //       } else {
  //         coalesceData[i][j][k].setCoalesced(true);
  //         avgPoint /= static_cast<float>(indices.size());
  //         avgTensor = avgTensor / static_cast<double>(indices.size());
  //         avgTensor = avgTensor * 1.5;
  //         coalesceData[i][j][k].setTensor(avgTensor);
  //         coalesceData[i][j][k].setPoint(avgPoint);
  //         coalesceData[i][j][k].setUncertaintyVal(unc_val_max);
  //         // coalescedTensors[i][j][k] = avgTensor;
  //         // coalescedPoints[i][j][k] = avgPoint;
  //         // coalescedUncertaintyVals[i][j][k] = unc_val_max;
  //         // addTensorToField(output, avgPoint, avgTensor);
  //       }
  //     }

  for (int i = 0; i < coalesce_count; ++i) {
      for (int i = 0; i < 3; ++i)
        dims[i] = std::max(static_cast<int>(dims[i]-1), 1);
    CoalesceData newData(dims);

    furtherCoalesce(coalesceData, newData, unc_threshold);
    coalesceData.addUncoalesced(output, neighborThres);
    // coalesceData.decrementDims();
    coalesceData = newData;
  }

  // furtherCoalesce(coalesceData, unc_threshold);
  // coalesceData.addUncoalesced(output);
  // coalesceData.decrementDims();
  addRemaining(coalesceData, output);


  // std::vector<std::vector<int>> newIndices;
  // for (int i = 0; i < dim[0]-2; ++i)
  //   for (int j = 0; j < dim[1]-2; ++j)
  //     for (int k = 0; k < std::max(static_cast<int>(dim[2]-2), 1); ++k) {
  //       if (imesh->is_imagemesh()) {
  //         newIndices = {
  //         {i, j, 0},
  //         {i, j+1, 0},
  //         {i+1, j, 0},
  //         {i+1, j+1, 0}};
  //       }
  //       if (imesh->is_latvolmesh()) {
  //         newIndices = {
  //           {i, j, k},
  //           {i, j, k+1},
  //           {i, j+1, k},
  //           {i, j+1, k+1},
  //           {i+1, j, k},
  //           {i+1, j, k+1},
  //           {i+1, j+1, k},
  //           {i+1, j+1, k+1}};
  //       }
  //       std::cout << "new in size " << newIndices.size() << "\n";
  //       std::cout << "new in vec size " << newIndices[0].size() << "\n";
  //       std::vector<Point> points(newIndices.size());
  //       std::vector<Dyadic3DTensor> tensors(newIndices.size());
  //       auto avgPoint = Point();
  //       auto avgTensor = Dyadic3DTensor();
  //       float unc_val_max = 0;
  //       printf("1\n");
  //       for (int p = 0; p < newIndices.size(); ++p) {
  //         points[p] = coalescedPoints[newIndices[p][i]][newIndices[p][j]][newIndices[p][k]];
  //       printf("2\n");
  //         avgPoint += points[p];
  //         // tensors[p] = getTensor(ifield, indices[p]);
  //         std::cout << "new in " << newIndices[p][i] << ", " << newIndices[p][j] << ", " << newIndices[p][k] << "\n";
  //         tensors[p] = coalescedTensors[newIndices[p][i]][newIndices[p][j]][newIndices[p][k]];
  //       printf("3\n");
  //         avgTensor += tensors[p];
  //         unc_val_max = std::max(unc_val_max, coalescedUncertaintyVals[newIndices[p][i]][newIndices[p][j]][newIndices[p][k]]);
  //       }
  //       printf("4\n");
  //       if (unc_val_max > unc_threshold) {
  //         // coalesced[i][j][k] = false;
  //         for (int p = 0; p < indices.size(); ++p) {
  //           addTensorToField(output, points[p], tensors[p]);
  //         }
  //       } else {
  //         // coalesced[i][j][k] = true;
  //         avgPoint /= static_cast<float>(indices.size());
  //         avgTensor = avgTensor / static_cast<double>(indices.size());
  //         avgTensor = avgTensor * 1.5;
  //         coalescedTensors[i][j][k] = avgTensor;
  //         addTensorToField(output, avgPoint, avgTensor);
  //       }
  // }

  return true;
}
