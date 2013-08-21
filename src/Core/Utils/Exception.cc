/*
 For more information, please see: http://software.sci.utah.edu

 The MIT License

 Copyright (c) 2012 Scientific Computing and Imaging Institute,
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

#include <Core/Utils/Exception.h>
#include <boost/regex.hpp>

using namespace SCIRun::Core;

const char* ExceptionBase::what() const throw()
{
  try
  {
    if (auto msg = boost::get_error_info<Core::ErrorMessage>(*this))
      return msg->c_str();
    else
      return "";
  }
  catch (...)
  {
    return "";
  }
}

std::string ExceptionBase::typeName() const
{
  try
  {
    //TODO very hacky.
    std::string name = typeid(*this).name();

    static boost::regex r(".*class.*<.*struct (.+)>");
    boost::smatch match;
    if (boost::regex_match(name, match, r))
      return (std::string)match[1];
    return name;
  }
  catch (...)
  {
    return "<Unknown>";
  }
}