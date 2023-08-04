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

using namespace SCIRun;
using namespace Core;
using namespace Core::Datatypes;
using namespace Core::Algorithms;
using namespace Core::Algorithms::Math;

class TensorModeler {
  public:
    enum FitMethod { OLS, WLS };
  private:
    MatrixHandle dwi_, bvecs_, bvals_;
    FitMethod fit_;
    size_t numGradients_;
    size_t numVoxels_;

    Eigen::VectorXd computeLnSignal0() {
      Eigen::VectorXd s0(numVoxels_);
      int num = 0;
      for (int i = 0; i < numGradients_; ++i)
        if (bvals_->get(i,0) == 0) {
          num++;
          for (int j = 0; j < numVoxels_; ++j)
            s0(j) += dwi_->get(i,j);
        }
      return s0/num; // Returns average
    }

    Eigen::MatrixXd computeLnSignal() {
      Eigen::VectorXd s0 = computeLnSignal0();
      Eigen::MatrixXd y(numGradients_, numVoxels_);

      for (int i = 0; i < numGradients_; ++i) {
        if (bvals_->get(i,0) != 0) {
          for (int j = 0; j < numVoxels_; ++j) {
            y(i,j) = std::log(dwi_->get(i,j));
          }
        }
      }
      return y;
    }

    Eigen::MatrixXd computeDesignMatrix() {
      Eigen::MatrixXd design(7,numGradients_);
      for (int i = 0; i < numGradients_; ++i) {
        auto bval = bvals_->get(i,0);
        auto bvec = Eigen::Vector3d(bvecs_->get(i,0), bvecs_->get(i,1), bvecs_->get(i,2));
        design(0, i) = -bval*bvec[0]*bvec[0];
        design(1, i) = -bval*bvec[1]*bvec[1];
        design(2, i) = -bval*bvec[2]*bvec[2];
        design(3, i) = -2*bval*bvec[0]*bvec[1];
        design(4, i) = -2*bval*bvec[1]*bvec[2];
        design(5, i) = -2*bval*bvec[0]*bvec[2];
        design(6, i) = -1;
      }
      return design;
    }

    void computeOLS() {
      auto x = computeDesignMatrix();
      auto y = computeLnSignal();
      Eigen::MatrixXd ols = (x.transpose()*x).inverse() * x.transpose() * y;
      std::cout << "ols dim: " << ols.cols() << ", "<<ols.rows()<<"\n";
    }

  public:
    TensorModeler(MatrixHandle dwi, MatrixHandle bvecs, MatrixHandle bvals, FitMethod fit) : dwi_(dwi), bvecs_(bvecs), bvals_(bvals), fit_(fit), numGradients_(dwi->ncols()), numVoxels_(dwi->nrows()) {
      computeDesignMatrix();
    }

    FieldHandle getTensors() {
      std::cout << "numVoxels: " << numVoxels_ << "\n";
      std::cout << "numGradients: " << numGradients_ << "\n";
      // computeLnSignal();
      computeOLS();
      return nullptr;
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
  std::cout << "in algo run\n";
  TensorModeler modeler(dwi, bvecs, bvals, TensorModeler::FitMethod::WLS);
  auto tensors = modeler.getTensors();
  auto residuals = modeler.getResiduals();

  AlgorithmOutput output;
  return output;
}
