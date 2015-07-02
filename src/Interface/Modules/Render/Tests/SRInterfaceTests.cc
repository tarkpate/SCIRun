/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2015 Scientific Computing and Imaging Institute,
   University of Utah.

   License for the specific language governing rights and limitations under
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

#if _WIN32
#include <GL/glew.h>
#endif
#include <gtest/gtest.h>
#include <Interface/Modules/Render/ES/SRInterface.h>

using namespace SCIRun;
using namespace SCIRun::Render;
using namespace SCIRun::Gui;

TEST(SRInterfaceTest, CanInstantiateSRInterface)
{
  std::shared_ptr<GLContext> context;
  std::vector<std::string> shaderDirs;
  SRInterface srinterface(context, shaderDirs);
}

class DummyGLContext : public GLContext
{
public:
  DummyGLContext() : GLContext(nullptr) {}
  virtual void makeCurrent() override 
  { 
    std::cout << "DummyGLContext::makeCurrent called" << std::endl;
  }
  virtual void swapBuffers() override
  {
    std::cout << "DummyGLContext::swapBuffers called" << std::endl;
  }
};


TEST(SRInterfaceTest, CanRenderEmptyFrame)
{
  std::shared_ptr<GLContext> context(new DummyGLContext);
  std::vector<std::string> shaderDirs;
  SRInterface srinterface(context, shaderDirs);

  srinterface.doFrame(0, 50);
}

