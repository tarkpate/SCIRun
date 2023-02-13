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


#include <Modules/Math/DownsampleTensorField.h>
#include <Core/Datatypes/Legacy/Field/Field.h>
#include <Core/Algorithms/Math/DownsampleTensorFieldAlgorithm.h>
#include <Core/Datatypes/DenseMatrix.h>

using namespace SCIRun;
using namespace Core::Algorithms::Math;
using namespace Core::Datatypes;
using namespace Modules::Math;
using namespace Dataflow::Networks;

MODULE_INFO_DEF(DownsampleTensorField, Math, SCIRun)

DownsampleTensorField::DownsampleTensorField() : Module(staticInfo_)
{
  INITIALIZE_PORT(InputField);
  INITIALIZE_PORT(MeanTensorField);
  INITIALIZE_PORT(VariationField);
}

void DownsampleTensorField::setStateDefaults()
{
  auto state = get_state();
  setStateStringFromAlgoOption(Parameters::MeanTensorInvariantMethod);
  setStateStringFromAlgoOption(Parameters::MeanTensorOrientationMethod);
}

void DownsampleTensorField::execute()
{
  std::cout << "exec\n";
  auto field = getRequiredInput(InputField);

  if (needToExecute())
  {
    setAlgoOptionFromState(Parameters::MeanTensorInvariantMethod);
    setAlgoOptionFromState(Parameters::MeanTensorOrientationMethod);
    std::cout << "run algo\n";
    auto output = algo().run(withInputData((InputField, field)));
    sendOutputFromAlgorithm(MeanTensorField, output);
    sendOutputFromAlgorithm(VariationField, output);
  }
}
