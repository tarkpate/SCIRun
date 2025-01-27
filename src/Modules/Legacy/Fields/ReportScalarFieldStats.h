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


#ifndef MODULES_LEGACY_FIELDS_REPORTSCALARFIELDSTATS_H__
#define MODULES_LEGACY_FIELDS_REPORTSCALARFIELDSTATS_H__

#include <Dataflow/Network/Module.h>
#include <Modules/Legacy/Fields/share.h>

namespace SCIRun {
  namespace Core
  {
    namespace Algorithms
    {
      namespace Fields
      {
        ALGORITHM_PARAMETER_DECL(AutoRangeEnabled);
        ALGORITHM_PARAMETER_DECL(Mean);
        ALGORITHM_PARAMETER_DECL(Median);
        ALGORITHM_PARAMETER_DECL(StandardDeviation);
        ALGORITHM_PARAMETER_DECL(HistogramBinCount);
      }
    }
  }
  namespace Modules {
    namespace Fields {

      class SCISHARE ReportScalarFieldStats : public Dataflow::Networks::Module,
        public Has1InputPort<FieldPortTag>,
        public Has1OutputPort<MatrixPortTag>
      {
      public:
        ReportScalarFieldStats();

        void execute() override;
        void setStateDefaults() override;

        INPUT_PORT(0, InputField, Field);
        OUTPUT_PORT(0, HistogramData, Matrix);

        MODULE_TRAITS_AND_INFO(ModuleFlags::ModuleHasUI)
      };
    }
  }
}

#endif
