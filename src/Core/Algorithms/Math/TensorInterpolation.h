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


#ifndef CORE_ALGORITHMS_FIELDS_TENSORINTERPOLATION_H
#define CORE_ALGORITHMS_FIELDS_TENSORINTERPOLATION_H

#include <Core/Datatypes/Dyadic3DTensor.h>
#include <Core/Algorithms/Base/AlgorithmBase.h>
#include <Core/Algorithms/Field/share.h>

namespace SCIRun {
  namespace Core {
    namespace Algorithms {
      namespace Fields {

        // ALGORITHM_PARAMETER_DECL(InvariantMethod);
        // ALGORITHM_PARAMETER_DECL(OrientMethod);

        class SCISHARE TensorInterpolationAlgo// : public AlgorithmBase
        {
        public:
          enum Method {
            MATRIX_AVERAGE,
            LOG_EUCLIDEAN,
            LINEAR_INVARIANT
          };

          private:
          Method invariantMethod_, orientMethod_;
          Core::Datatypes::Dyadic3DTensor matrixAverage(std::vector<Core::Datatypes::Dyadic3DTensor>& tensors, std::optional<std::vector<float>> weights);
          Core::Datatypes::Dyadic3DTensor linearInvariant(std::vector<Core::Datatypes::Dyadic3DTensor>& tensors, std::optional<std::vector<float>> weights);
          Core::Datatypes::Dyadic3DTensor logEuclidean(std::vector<Core::Datatypes::Dyadic3DTensor>& tensors, std::optional<std::vector<float>> weights);

          public:
          TensorInterpolationAlgo(Method method);
          TensorInterpolationAlgo(Method invariantMethod, Method orientMethod);
          Core::Datatypes::Dyadic3DTensor interpolate(std::vector<Core::Datatypes::Dyadic3DTensor>& tensors);
          Core::Datatypes::Dyadic3DTensor interpolate(std::vector<Core::Datatypes::Dyadic3DTensor>& tensors, std::optional<std::vector<float>> weights);

          // virtual AlgorithmOutput run_generic(const AlgorithmInput& input) const;
        };

      }}}}

#endif
