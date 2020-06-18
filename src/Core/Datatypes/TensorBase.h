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

#ifndef CORE_DATATYPES_TENSOR_BASE_H
#define CORE_DATATYPES_TENSOR_BASE_H

#include <unsupported/Eigen/CXX11/TensorSymmetry>
#include "unsupported/Eigen/CXX11/src/TensorSymmetry/Symmetry.h"
#include <Core/Datatypes/share.h>

namespace SCIRun {
namespace Core {
  namespace Datatypes {
    template <typename T, int Rank>
    class TensorBase : public Eigen::Tensor<T, Rank>
    {
      typedef Eigen::Tensor<T, Rank> parent;

     public:
      using parent::Tensor;  // adding parent constructors

      bool operator==(const TensorBase<T, Rank>& t) { return dimensionsEqual(t) && valuesEqual(t); }
      bool operator!=(const TensorBase<T, Rank>& t) { return !(*this == t); }

      bool dimensionsEqual(const TensorBase<T, Rank>& t)
      {
        if (parent::NumDimensions != t.NumDimensions) return false;

        auto thisDim = parent::dimensions();
        auto thatDim = t.dimensions();

        for (int i = 0; i < parent::NumDimensions; ++i)
          if (thisDim[i] != thatDim[i]) return false;

        return true;
      }

      bool valuesEqual(const TensorBase<T, Rank>& t)
      {
        auto dim = parent::dimensions();
        auto index = std::vector<int>(parent::NumDimensions, 0);
        return checkForEquals(index, dim, t);
      }

      bool checkForEquals(
          std::vector<int>& index, Eigen::DSizes<long int, Rank>& dim, const TensorBase<T, Rank>& t)
      {
        if ((*this)(index) != t(index)) return false;

        if (incrementIndex(index, dim))
          return checkForEquals(index, dim, t);
        else
          return true;
      }

      bool incrementIndex(std::vector<int>& index, Eigen::DSizes<long int, Rank>& dim)
      {
        int d = index.size() - 1;

        while (d >= 0)
        {
          auto dIndex = index[d] + 1;
          if (dIndex >= dim[d])
            --d;
          else
          {
            index[d] = dIndex;
            break;
          }
        }

        if (d == 0 && (index[0] + 1 >= dim[0])) return false;
        if (d != index.size() - 1)
        {
          for (int i = d + 1; i < index.size(); ++i)
            index[i] = 0;
        }
        return true;
      }
    };
  }
}
}

#endif
