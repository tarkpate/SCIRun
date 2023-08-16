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
  bool valid;
  int neighbors_coalesced;
  Dyadic3DTensor t;
  std::vector<Dyadic3DTensor> t_samples;
  Point p;
  public:
    CoalesceDataNode() {
      valid = false;
      neighbors_coalesced = 0;
      t = Dyadic3DTensor();
      p = Point();
      t_samples = {t};
    }
    void resetNeighbors() { neighbors_coalesced=0; }
    void setValid(bool val) { valid=val; }
    void addNeighborCoalesced() { neighbors_coalesced++; }
    void setTensor(Dyadic3DTensor val) { t=val; }
    void setTensors(std::vector<Dyadic3DTensor> val) {
      t_samples.resize(val.size());
      for(int i = 0; i < val.size(); ++i)
        t_samples[i]=val[i];
    }
    void setPoint(Point val) { p=val; }
    bool isValid() const { return valid; }
    int getNeighborsCoalesced() const { return neighbors_coalesced; }
    size_t getNumSamples() const { return t_samples.size(); }
    const Dyadic3DTensor& getTensor() const { return t; }
    Dyadic3DTensor getTensorSample(int i) const { return t_samples[i]; }
    const Point& getPoint() const { return p; }

    CoalesceDataNode& operator+(const CoalesceDataNode& other) {
      auto otherT = other.getTensor();
      t += otherT;
      for (int i = 0; i < 3; ++i)
        p[i] += other.getPoint()[i];
      return *this;
    }

    CoalesceDataNode& operator/(double div) {
      t = t / div;
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
    }

  public:
    CoalesceData(const data_index_type& in_dims) {
      for (int i = 0; i < 3; ++i) {
        dims[i] = in_dims[i];
        data_dims[i] = in_dims[i];
      }
      data = std::vector<CoalesceDataNode>(data_dims[0]*data_dims[1]*data_dims[2]);
      for (int i = 0; i < data.size(); ++i)
        data[i].setValid(false);
    }

    CoalesceDataNode* getNode(int x, int y, int z) { return &data[getArrayIndex(x, y, z)]; }

    void addUncoalesced(FieldHandle& output, int neighborThres) {
      for (int k = 0; k < dims[2]; k++)
        for (int j = 0; j < dims[1]; j++)
          for (int i = 0; i < dims[0]; i++) {
            auto n = getNode(i,j,k);
            if (n->isValid() && n->getNeighborsCoalesced() < neighborThres) {
              addTensorToField(output, n->getPoint(), n->getTensor());
            }
          }
    }

    void makeNode(index_type idx, const Point& p, const Dyadic3DTensor& t, bool valid) {
      CoalesceDataNode& n = data[idx];
      n.setPoint(p);
      n.setTensor(t);
      n.setValid(valid);
      n.setTensors({t});
    }

    const data_index_type& getDims() const { return dims; }
    const Dyadic3DTensor& getTensor(int x, int y, int z) { return getNode(x, y, z)->getTensor(); }
    const Point& getPoint(int x, int y, int z) { return getNode(x, y, z)->getPoint(); }
};

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

Dyadic3DTensor getTensor(const VField* vfield, VMesh::Node::index_type idx) {
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

void addRemaining(CoalesceData* coalesceData, FieldHandle& output) {
  const CoalesceData::data_index_type& dims = coalesceData->getDims();
  for (int k = 0; k < dims[2]; k++)
    for (int j = 0; j < dims[1]; j++)
      for (int i = 0; i < dims[0]; i++) {
        auto n = coalesceData->getNode(i,j,k);
        if (n->isValid())
          addTensorToField(output, n->getPoint(), n->getTensor());
      }
}

void furtherCoalesce(CoalesceData* coalesceData, CoalesceData* newData, float unc_threshold, int blockSize, int overlapSize, int borderSize, int coalesceCount) {


  TensorInterpolationAlgo interpAlgo(TensorInterpolationAlgo::Method::LINEAR_INVARIANT,
                                     TensorInterpolationAlgo::Method::LOG_EUCLIDEAN);

  const CoalesceData::data_index_type& dims = coalesceData->getDims();
  const CoalesceData::data_index_type& newDims = newData->getDims();

  bool x_flat = dims[0] <= 1;
  bool y_flat = dims[1] <= 1;
  bool z_flat = dims[2] <= 1;
  int incSize = blockSize - overlapSize;

  for (int k = 0; k < newDims[2]; k++)
    for (int j = 0; j < newDims[1]; j++)
      for (int i = 0; i < newDims[0]; i++) {
        int sampleX = x_flat ? 0 : i*incSize + borderSize;
        int sampleY = y_flat ? 0 : j*incSize + borderSize;
        int sampleZ = z_flat ? 0 : k*incSize + borderSize;

        bool valid = true;
        std::vector<CoalesceDataNode*> nodes;
        for (int blockX = x_flat ? 0 : -blockSize/2; x_flat ? blockX < 1 : blockX <= blockSize/2; ++blockX)
          for (int blockY = y_flat ? 0 : -blockSize/2; y_flat ? blockY < 1 : blockY <= blockSize/2; ++blockY)
            for (int blockZ = z_flat ? 0 : -blockSize/2; z_flat ? blockZ < 1 : blockZ <= blockSize/2; ++blockZ) {
              auto n = coalesceData->getNode(sampleX+blockX, sampleY+blockY, sampleZ+blockZ);
              valid = valid && n->isValid();
              nodes.push_back(n);
          }

        if (valid) {
        int numSamples = nodes[0]->getNumSamples();
        std::vector<Dyadic3DTensor> tensors(nodes.size()*numSamples);
        std::vector<Dyadic3DTensor::VectorType> eigvecs(nodes.size()*numSamples);
        Point avgPoint = Point(0, 0, 0);
        for (int n = 0; n < nodes.size(); ++n) {
          auto node = nodes[n];
          for (int s = 0; s < numSamples; ++s) {
            tensors[n*numSamples+s] = node->getTensorSample(s);
            eigvecs[n*numSamples+s] = tensors[n*numSamples+s].getEigenvector(0);
          }
          avgPoint += node->getPoint();
        }
        avgPoint /= nodes.size();

        Dyadic3DTensor::VectorType avgVec(0);
        for (auto e : eigvecs)
          avgVec += e.normalized();
        avgVec.normalize();

        std::vector<float> angleDiff(eigvecs.size());
        for (int c = 0; c < eigvecs.size(); ++c)
          angleDiff[c] = 1.0-abs(avgVec.dot(eigvecs[c]));

        std::vector<float> uncVals(tensors.size());
        float avgUncVal = 0;
        for (int c = 0; c < tensors.size(); ++c) {
          auto fa = tensors[c].fractionalAnisotropy();
          auto angle = abs(angleDiff[c]);
          if (fa < 0.2) // highly isotropic
            uncVals[c] = fa;
          else
            uncVals[c] = angle;
          avgUncVal += uncVals[c];
        }
        avgUncVal /= tensors.size();
        Dyadic3DTensor avgTensor = interpAlgo.interpolate(tensors);
        avgTensor = avgTensor * std::pow(2, coalesceCount+1);

        CoalesceDataNode* avgNode = newData->getNode(i,j,k);
        if (avgUncVal < unc_threshold) {
          avgNode->setValid(true);
          avgNode->setTensor(avgTensor);
          avgNode->setPoint(avgPoint);
          avgNode->setTensors(tensors);
          coalesceData->getNode(sampleX,sampleY,sampleZ)->setValid(false);
          for (auto n : nodes) n->addNeighborCoalesced();
        }
        else{
          avgNode->setValid(false);
        }
        }
    }
}

// General access function
bool
CoalesceMeshAlgo::runImpl(const FieldHandle& input, const FieldHandle& isoValueField, FieldHandle& output) const
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

	const FieldInformation fi(input);
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

  const VField* ifield = input->vfield();
  VMesh*  imesh = input->vmesh();

  VMesh::dimension_type mesh_dims;
  imesh->get_dimensions(mesh_dims);
  CoalesceData::data_index_type dims;
  for (int i = 0; i < DIM_SIZE_; ++i)
    dims[i] = std::max(static_cast<int>(mesh_dims[i]), 1);

  std::vector<VMesh::Node::index_type> indices;
  if (imesh->is_imagemesh())
    dims[2] = 1;

  CoalesceData* currData = new CoalesceData(dims);
  CoalesceData::data_index_type idx;
  for (int i = 0; i < dims[0]*dims[1]*dims[2]; ++i)
        currData->makeNode(i, getPoint(imesh, i),
                           getTensor(ifield, i),
                           true);
  VMesh::index_type borderSize = blockSize / 2; // Remainder ignored
  CoalesceData::data_index_type newDims;

  CoalesceData* newData;
  for (int c = 0; c < coalesceCount; ++c) {
    int incSize = blockSize - overlapSize;

    for (int i = 0; i < 3; ++i)
      newDims[i] = std::max(1, static_cast<int>((dims[i]-borderSize)/incSize));
    newData = new CoalesceData(newDims);

    furtherCoalesce(currData, newData, uncThreshold, blockSize, overlapSize, borderSize, c);
    currData->addUncoalesced(output, neighborThres);

    currData = newData;
    for (int i = 0; i < 3; ++i)
      dims[i] = newDims[i];

  }

  addRemaining(currData, output);
  delete currData;

  return true;
}
