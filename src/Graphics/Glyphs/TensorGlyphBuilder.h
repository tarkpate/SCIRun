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

#ifndef Graphics_Glyphs_TENSOR_GLYPH_BUILDER_H
#define Graphics_Glyphs_TENSOR_GLYPH_BUILDER_H

#include <Core/Datatypes/Color.h>
#include <Core/Datatypes/Dyadic3DTensor.h>
#include <Core/GeometryPrimitives/Point.h>
#include <Core/GeometryPrimitives/Transform.h>
#include <Core/GeometryPrimitives/Vector.h>
#include <Core/Math/TrigTable.h>
#include <Graphics/Glyphs/GlyphConstructor.h>
#include <Graphics/Glyphs/share.h>

namespace SCIRun {
namespace Graphics {
  struct SCISHARE EllipsoidPointParams
  {
    double sinPhi;
    double cosPhi;
    double sinTheta;
    double cosTheta;
  };

  struct SCISHARE SuperquadricPointParams
  {
    double sinPhi;
    double cosPhi;
    double sinTheta;
    double cosTheta;
    double A;
    double B;
  };

  class SCISHARE TensorGlyphBuilder
  {
   public:
    TensorGlyphBuilder(
        const Core::Datatypes::Dyadic3DTensor& t, const Core::Geometry::Point& center);
    void scaleTensor(double scale);
    void makeTensorPositive(bool reorder = false, bool makeGlyph = true);
    void normalizeTensor();
    void setColor(const Core::Datatypes::ColorRGB& color);
    void setResolution(double resolution);
    void generateSuperquadricTensor(GlyphConstructor& constructor, double emphasis);
    void generateEllipsoid(GlyphConstructor& constructor, bool half);
    void generateBox(GlyphConstructor& constructor);

    Core::Geometry::Point evaluateSuperquadricPointLinear(const SuperquadricPointParams& params);
    Core::Geometry::Point evaluateSuperquadricPointPlanar(const SuperquadricPointParams& params);
    Core::Geometry::Point evaluateSuperquadricPoint(
        bool linear, const SuperquadricPointParams& params);
    Core::Geometry::Point evaluateEllipsoidPoint(EllipsoidPointParams& params);

    Core::Geometry::Transform getScale();
    void setTensor(const Core::Datatypes::Dyadic3DTensor& t);
    Core::Datatypes::Dyadic3DTensor getTensor() const;

   private:
    void generateBoxSide(GlyphConstructor& constructor,
        const std::vector<Core::Geometry::Vector>& points, const Core::Geometry::Vector& normal);
    std::vector<Core::Geometry::Vector> generateBoxPoints();
    std::pair<double, double> getAAndB(double emphasis);
    bool isLinear();
    void computeSinCosTable(bool half);
    void computeTransforms();
    void postScaleTransorms();

    const static int DIMENSIONS_ = 3;
    const static int BOX_FACE_POINTS_ = 4;
    Core::Datatypes::Dyadic3DTensor t_;
    Core::Geometry::Point center_;
    Core::Geometry::Transform trans_, rotate_;
    Core::Datatypes::ColorRGB color_ = {1.0, 1.0, 1.0};
    int resolution_ = 10;
    bool flatTensor_ = false;
    Core::Geometry::Vector zeroNorm_ = {0, 0, 0};
    SinCosTable tab1_, tab2_;
    int nv_ = 0;
    int nu_ = 0;
    double cl_, cp_;
  };
}
}

#endif
