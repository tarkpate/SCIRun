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

//STL classes needed
#include <algorithm>
#include <set>

using namespace SCIRun;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Algorithms::Fields;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Core::Geometry;
using namespace SCIRun::Core::Logging;


ALGORITHM_PARAMETER_DEF(Fields, AddConstraints);
ALGORITHM_PARAMETER_DEF(Fields, CoalesceMethod);
ALGORITHM_PARAMETER_DEF(Fields, IsoValue);

CoalesceMeshAlgo::CoalesceMeshAlgo()
{
		using namespace Parameters;
		addOption(AddConstraints,"all","all|greaterthan|unequal|lessthan|none");
		addOption(CoalesceMethod,"Default","Default|Expand coalescement volume to improve element quality");
		addParameter(IsoValue,0.0);
}

AlgorithmOutput CoalesceMeshAlgo::run(const AlgorithmInput& input) const
{
	auto field = input.get<Field>(Variables::InputField);
  FieldHandle outputField;

  if (!runImpl(field, outputField))
    THROW_ALGORITHM_PROCESSING_ERROR("False returned on legacy run call.");
	AlgorithmOutput output;
	output[Variables::OutputField] = outputField;
  return output;
}

// General access function

bool
CoalesceMeshAlgo::runImpl(FieldHandle input, FieldHandle& output) const
{
	REPORT_STATUS(CoalesceMesh);
  if (!input)
  {
    error("No input field");
    return (false);
  }
	FieldInformation fi(input);
	FieldInformation fo(output);

	const std::string rMethod = getOption(Parameters::CoalesceMethod);
	const double isoVal = get(Parameters::IsoValue).toDouble();
	const std::string addCon = getOption(Parameters::AddConstraints);

  if (input->vfield()->num_values() == 0)
  {
    error("Input field has no data values. The CoalesceMesh algorithm requires input fields to contain data.");
    return (false);
  }

  if (addCon == "none")
  {
    output = input;
    return (true);
  }

  if (fi.is_pnt_element() || fi.is_prism_element())
  {
    error("This algorithm does not support point or prism meshes");
    return(false);
  }

  if ((!(fi.is_scalar())) && (addCon != "all"))
  {
    error("Field data needs to be of scalar type");
    return (false);
  }

  error("No coalescement method has been implemented for this type of mesh");
  return (false);
}
