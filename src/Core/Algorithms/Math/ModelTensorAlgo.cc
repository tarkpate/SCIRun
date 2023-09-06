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


#include <iostream>
#include <cmath>
#include <Core/Algorithms/Math/ModelTensorAlgo.h>
#include <Core/Datatypes/Matrix.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/Algorithms/Base/AlgorithmPreconditions.h>
#include <Core/Datatypes/Legacy/Field/FieldInformation.h>
#include <Core/Datatypes/Mesh/MeshFacade.h>
#include <Core/GeometryPrimitives/Point.h>
#include <Core/GeometryPrimitives/Tensor.h>
#include <Core/Datatypes/Legacy/Field/Field.h>
#include <Core/Datatypes/Legacy/Field/VField.h>
#include <math.h>
#include <unsupported/Eigen/MatrixFunctions>

using namespace SCIRun;
using namespace Core;
using namespace Core::Algorithms;
using namespace Core::Algorithms::Math;
using namespace Core::Datatypes;
using namespace Core::Geometry;

class TensorModeler {
  public:
    enum FitMethod { OLS, WLS };
  private:
    MatrixHandle dwi_, bvecs_, bvals_;
    FitMethod fit_;
    size_t numGradients_;
    size_t numVoxels_;

    Eigen::MatrixXd computeLnSignal() {
      Eigen::MatrixXd y(numGradients_, numVoxels_);

      for (int i = 0; i < numGradients_; ++i)
          for (int j = 0; j < numVoxels_; ++j) {
            auto f = dwi_->get(j,i);
            y(i,j) = f > 0 ? -std::log(f) : 0;
          }

      return y;
    }

    Eigen::MatrixXd computeDesignMatrix() {
      Eigen::MatrixXd design(numGradients_,7);
      for (int i = 0; i < numGradients_; ++i) {
        auto bval = bvals_->get(0,i);
        auto bvec = Eigen::Vector3d(bvecs_->get(0,i), bvecs_->get(1,i), bvecs_->get(2,i));
        design(i, 0) = bval*bvec[0]*bvec[0];
        design(i, 1) = bval*bvec[1]*bvec[1];
        design(i, 2) = bval*bvec[2]*bvec[2];
        design(i, 3) = 2*bval*bvec[0]*bvec[1];
        design(i, 4) = 2*bval*bvec[0]*bvec[2];
        design(i, 5) = 2*bval*bvec[1]*bvec[2];
        design(i, 6) = -1;
      }
      return design;
    }

    Eigen::MatrixXd computeOLS() {
      auto x = computeDesignMatrix();
      auto y = computeLnSignal();
      return x.completeOrthogonalDecomposition().pseudoInverse() * y;
    }

  public:
    TensorModeler(MatrixHandle dwi, MatrixHandle bvecs, MatrixHandle bvals, FitMethod fit) : dwi_(dwi), bvecs_(bvecs), bvals_(bvals), fit_(fit), numGradients_(dwi->ncols()), numVoxels_(dwi->nrows()) {
    }

    FieldHandle getTensors() {
      FieldInformation fi("LatVolMesh", "linear", "tensor");
      Eigen::MatrixXd ols = computeOLS();
      // TODO move to user input
      std::vector<int> shape = {81, 106, 76};

      Point min(0,0,0);
      Point max(shape[0],shape[1],shape[2]);
      MeshHandle mesh = CreateMesh(fi, shape[0], shape[1], shape[2], min, max);
      FieldHandle output = CreateField(fi, mesh);

      for (int k = 0; k < shape[2]; ++k)
        for (int j = 0; j < shape[1]; ++j)
          for (int i = 0; i < shape[0]; ++i) {
            auto ols_index = i*shape[1]*shape[2] + j*shape[2] + k;
            auto out_index = k*shape[1]*shape[0] + j*shape[0] + i;
            Tensor t(ols(0,ols_index), ols(3,ols_index), ols(4,ols_index),
                     ols(1,ols_index), ols(5,ols_index), ols(2,ols_index));
            output->vfield()->set_value(t, out_index);
      }

      return output;
    }

    MatrixHandle getResiduals() { return nullptr; }
};

ModelTensorAlgo::ModelTensorAlgo()
{
}

AlgorithmOutput ModelTensorAlgo::run(const AlgorithmInput& input) const
{
  const auto dwi = input.get<Matrix>(Variables::FirstMatrix);
  const auto bvecs = input.get<Matrix>(Variables::SecondMatrix);
  const auto bvals = input.get<Matrix>(Variables::ThirdMatrix);
  TensorModeler modeler(dwi, bvecs, bvals, TensorModeler::FitMethod::WLS);
  auto tensors = modeler.getTensors();
  auto residuals = modeler.getResiduals();

  AlgorithmOutput output;
  output[Variables::OutputField] = tensors;
  return output;
}
