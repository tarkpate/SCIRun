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


#include <Modules/Python/LoopEnd.h>
#include <Modules/Python/LoopIncrement.h>
#include <Modules/Python/PythonObjectForwarder.h>
#ifdef BUILD_WITH_PYTHON
#include <Modules/Python/PythonInterfaceParser.h>
#include <Core/Python/PythonInterpreter.h>
#include <Core/Logging/Log.h>
#include <Core/Datatypes/MetadataObject.h>
#include <Dataflow/Engine/Python/NetworkEditorPythonAPI.h>
#endif

using namespace SCIRun;
using namespace SCIRun::Modules::Python;
using namespace SCIRun::Core::Datatypes;
using namespace SCIRun::Dataflow::Networks;
using namespace SCIRun::Core;
using namespace SCIRun::Core::Thread;
using namespace SCIRun::Core::Algorithms;
using namespace SCIRun::Core::Algorithms::Python;

ALGORITHM_PARAMETER_DEF(Python, LoopEndCode);
ALGORITHM_PARAMETER_DEF(Python, LoopWhileCondition);

MODULE_INFO_DEF(LoopEnd, Python, SCIRun)

//Mutex InterfaceWithPython::lock_("InterfaceWithPython");
//bool InterfaceWithPython::matlabInitialized_{ false };

LoopEnd::LoopEnd() : Module(staticInfo_)
{
  INITIALIZE_PORT(LoopEndCodeObject);
  INITIALIZE_PORT(InputMatrix);
  INITIALIZE_PORT(InputField);
  INITIALIZE_PORT(InputString);

#ifdef BUILD_WITH_PYTHON
  translator_.reset(new InterfaceWithPythonCodeTranslatorImpl([this]() { return id().id_; }, get_state(), { Parameters::LoopWhileCondition }));
#endif
}

void LoopEnd::setStateDefaults()
{
  auto state = get_state();

  state->setValue(Parameters::LoopEndCode, std::string("# Insert your Python code here. The SCIRun API package is automatically imported."));
  state->setValue(Parameters::PollingIntervalMilliseconds, 10);
  state->setValue(Parameters::NumberOfRetries, 5);
  state->setValue(Parameters::LoopWhileCondition, std::string("loopWhileCondition"));
}

void LoopEnd::postStateChangeInternalSignalHookup()
{
  setProgrammableInputPortEnabled(true);
}

//
// std::vector<AlgorithmParameterName> InterfaceWithPython::outputNameParameters()
// {
//   return { Parameters::PythonOutputMatrix1Name, Parameters::PythonOutputMatrix2Name, Parameters::PythonOutputMatrix3Name,
//     Parameters::PythonOutputField1Name, Parameters::PythonOutputField2Name, Parameters::PythonOutputField3Name,
//     Parameters::PythonOutputString1Name, Parameters::PythonOutputString2Name, Parameters::PythonOutputString3Name };
// }
//
// std::vector<std::string> InterfaceWithPython::connectedPortIds() const
// {
//   std::vector<std::string> ids;
//   for (const auto& port : inputPorts())
//   {
//     if (port->nconnections() > 0)
//     {
//       ids.push_back(port->id().toString());
//     }
//   }
//   return ids;
// }

void LoopEnd::execute()
{
#ifdef BUILD_WITH_PYTHON
  auto matrices = getOptionalDynamicInputs(InputMatrix);
  auto fields = getOptionalDynamicInputs(InputField);
  auto strings = getOptionalDynamicInputs(InputString);
  if (needToExecute())
  {
    auto code = get_state()->getValue(Parameters::LoopEndCode).toString();
    remark(code);
    auto convertedCode = translator_->translate(code);
    PythonInterpreter::Instance().run_script(convertedCode.code);
    PythonObjectForwarderImpl<LoopEnd> impl(*this);
    DatatypeHandle whileCondition;
    if (oport_connected(LoopEndCodeObject))
      whileCondition = impl.waitForOutputFromTransientState(get_state()->getValue(Parameters::LoopWhileCondition).toString(), DummyPortName(), DummyPortName(), DummyPortName());

    if (whileCondition)
    {
      auto dense = dynamic_cast<DenseMatrix*>(whileCondition.get());
      auto endVal = dense->get(0,0);
      logCritical("Generated end condition, value = {}", endVal);
      if (endVal)
      {
        logCritical("Attempting to enqueue execute again...");
        enqueueExecuteAgain(true);
      }
    }
    else
      logCritical("End condition not generated.");
  }
#else
  error("This module does nothing, turn on BUILD_WITH_PYTHON to enable.");
#endif
}

bool LoopEnd::checkForVirtualConnection(const ModuleInterface& downstream) const
{
  return dynamic_cast<const LoopIncrement*>(&downstream) != nullptr;
}
