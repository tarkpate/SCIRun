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

#include <Core/GeometryPrimitives/Transform.h>
#include <Core/Math/MiscMath.h>
#include <Graphics/Glyphs/GlyphGeomUtility.h>
#include <Graphics/Glyphs/TensorGlyphBuilder.h>

using namespace SCIRun;
using namespace Graphics;
using namespace Core::Geometry;
using namespace Core::Datatypes;

Vector makeSCIRunVector(const Eigen::Vector3d& vec)
{
  return Vector(vec[0], vec[1], vec[2]);
}

TensorGlyphBuilder::TensorGlyphBuilder(const Dyadic3DTensor& t, const Point& center)
{
  t_ = t;
  center_ = Point(center);
  trans_ = Transform();
  rotate_ = Transform();
}

void TensorGlyphBuilder::setTensor(const Dyadic3DTensor& t)
{
  t_ = t;
}

Dyadic3DTensor TensorGlyphBuilder::getTensor() const
{
  return t_;
}

void TensorGlyphBuilder::scaleTensor(double scale)
{
  t_ = t_ * scale;
}

void TensorGlyphBuilder::normalizeTensor()
{
  t_.normalize();
}

void TensorGlyphBuilder::setColor(const ColorRGB& color)
{
  color_ = color;
}

void TensorGlyphBuilder::setResolution(double resolution)
{
  resolution_ = resolution;
}

void TensorGlyphBuilder::makeTensorPositive(bool reorder, bool makeGlyph)
{
  static const double zeroThreshold = 0.000001;

  auto eigvals = t_.getEigenvalues();
  auto eigvecs = t_.getEigenvectors();
  for (auto& e : eigvals)
  {
    e = fabs(e);
    if (e <= zeroThreshold) e = 0;
  }

  // These are exactly zero after thresholding
  flatTensor_ = eigvals[0] == 0 || eigvals[1] == 0 || eigvals[2] == 0;

  if (makeGlyph)
  {
    auto cross = eigvecs[0].cross(eigvecs[1]);
    if (cross.dot(eigvecs[2]) < 2e-12) eigvecs[2] = cross;
  }

  for (int d = 0; d < DIMENSIONS_; ++d)
    if (eigvals[d] == 0)
    {
      auto cross = eigvecs[(d + 1) % DIMENSIONS_].cross(eigvecs[(d + 2) % DIMENSIONS_]);
      cross /= cross.norm();
      eigvecs[d] = cross;
      zeroNorm_ = makeSCIRunVector(cross);
      break;
    }

  t_.setEigens(eigvecs, eigvals);
}

void TensorGlyphBuilder::computeSinCosTable(bool half)
{
  nv_ = resolution_;
  nu_ = resolution_ + 1;

  // Should only happen when doing half ellipsoids.
  if (half) nv_ /= 2;
  if (nv_ < 2) nv_ = 2;

  double end = half ? M_PI : 2.0 * M_PI;
  tab1_ = SinCosTable(nu_, 0, end);
  tab2_ = SinCosTable(nv_, 0, M_PI);
}

void TensorGlyphBuilder::computeTransforms()
{
  auto eigvecs = t_.getEigenvectors();
  for (auto& v : eigvecs)
    v.normalize();
  GlyphGeomUtility::generateTransforms(center_, makeSCIRunVector(eigvecs[0]),
      makeSCIRunVector(eigvecs[1]), makeSCIRunVector(eigvecs[2]), trans_, rotate_);
}

void TensorGlyphBuilder::postScaleTransorms()
{
  auto eigvals = t_.getEigenvalues();
  Vector eigvalsVector(eigvals[0], eigvals[1], eigvals[2]);

  trans_.post_scale(eigvalsVector);
}

void TensorGlyphBuilder::generateEllipsoid(GlyphConstructor& constructor, bool half)
{
  makeTensorPositive(true);
  computeTransforms();
  postScaleTransorms();
  computeSinCosTable(half);

  auto eigvals = t_.getEigenvalues();
  Vector eigvalsVector(eigvals[0], eigvals[1], eigvals[2]);
  Transform rotateThenInvScale = rotate_;
  rotateThenInvScale.post_scale(Vector(1.0, 1.0, 1.0) / eigvalsVector);

  for (int v = 0; v < nv_ - 1; ++v)
  {
    double sinPhi[2] = {tab2_.sin(v), tab2_.sin(v + 1)};
    double cosPhi[2] = {tab2_.cos(v), tab2_.cos(v + 1)};

    for (int u = 0; u < nu_; ++u)
    {
      EllipsoidPointParams params;
      params.sinTheta = tab1_.sin(u);
      params.cosTheta = tab1_.cos(u);

      // Transorm points and add to points list
      constructor.setOffset();
      for (int i = 0; i < 2; ++i)
      {
        params.sinPhi = sinPhi[i];
        params.cosPhi = cosPhi[i];
        Point point = evaluateEllipsoidPoint(params);
        Vector pVector = Vector(trans_ * point);

        Vector normal;
        if (flatTensor_)
        {
          // Avoids recalculating norm vector and prevents vectors with infinite length
          bool first_half = v < nv_ / 2;
          normal = first_half ? zeroNorm_ : -zeroNorm_;
        }
        else
        {
          normal = rotateThenInvScale * Vector(point);
          normal.safe_normalize();
        }

        constructor.addVertex(pVector, normal, color_);
      }

      constructor.addIndicesToOffset(0, 1, 2);
      constructor.addIndicesToOffset(2, 1, 3);
    }
    constructor.popIndicesNTimes(6);
  }
}

Point TensorGlyphBuilder::evaluateEllipsoidPoint(EllipsoidPointParams& params)
{
  double x, y, z;
  x = params.sinPhi * params.sinTheta;
  y = params.sinPhi * params.cosTheta;
  z = params.cosPhi;
  return Point(x, y, z);
}

bool TensorGlyphBuilder::isLinear()
{
  cl_ = t_.linearCertainty();
  cp_ = t_.planarCertainty();
  return cl_ >= cp_;
}

std::pair<double, double> TensorGlyphBuilder::getAAndB(double emphasis)
{
  double pPower = GlyphGeomUtility::spow((1.0 - cp_), emphasis);
  double lPower = GlyphGeomUtility::spow((1.0 - cl_), emphasis);
  if (isLinear())
    return std::make_pair(pPower, lPower);
  else
    return std::make_pair(lPower, pPower);
}

void TensorGlyphBuilder::generateSuperquadricTensor(GlyphConstructor& constructor, double emphasis)
{
  makeTensorPositive(true);
  computeTransforms();
  postScaleTransorms();
  computeSinCosTable(false);

  bool linear = isLinear();
  auto AAndB = getAAndB(emphasis);
  SuperquadricPointParams params;
  params.A = AAndB.first;
  params.B = AAndB.second;

  SuperquadricPointParams normalParams;
  normalParams.A = 2.0 - params.A;
  normalParams.B = 2.0 - params.B;

  for (int v = 0; v < nv_ - 1; ++v)
  {
    double sinPhi[2] = {tab2_.sin(v), tab2_.sin(v + 1)};
    double cosPhi[2] = {tab2_.cos(v), tab2_.cos(v + 1)};

    for (int u = 0; u < nu_; ++u)
    {
      constructor.setOffset();
      params.sinTheta = normalParams.sinTheta = tab1_.sin(u);
      params.cosTheta = normalParams.cosTheta = tab1_.cos(u);

      for (int i = 0; i < 2; ++i)
      {
        params.sinPhi = normalParams.sinPhi = sinPhi[i];
        params.cosPhi = normalParams.cosPhi = cosPhi[i];
        // Transorm points and add to points list
        Point p = evaluateSuperquadricPoint(linear, params);
        Vector pVector = Vector(trans_ * p);

        Vector normal;
        if (flatTensor_)
        {
          // Avoids recalculating norm vector and prevents vectors with infinite length
          bool first_half = v < nv_ / 2;
          normal = first_half ? zeroNorm_ : -zeroNorm_;
        }
        else
        {
          normal = Vector(evaluateSuperquadricPoint(linear, normalParams));
          normal = rotate_ * normal;
          normal.safe_normalize();
        }

        constructor.addVertex(pVector, normal, color_);
      }

      constructor.addIndicesToOffset(0, 1, 2);
      constructor.addIndicesToOffset(2, 1, 3);
    }
  }
  constructor.popIndicesNTimes(6);
}

Point TensorGlyphBuilder::evaluateSuperquadricPoint(
    bool linear, const SuperquadricPointParams& params)
{
  if (linear)
    return evaluateSuperquadricPointLinear(params);
  else
    return evaluateSuperquadricPointPlanar(params);
}

// Generate around x-axis
Point TensorGlyphBuilder::evaluateSuperquadricPointLinear(const SuperquadricPointParams& params)
{
  double x, y, z;
  x = GlyphGeomUtility::spow(params.cosPhi, params.B);
  y = -GlyphGeomUtility::spow(params.sinPhi, params.B) *
      GlyphGeomUtility::spow(params.sinTheta, params.A);
  z = GlyphGeomUtility::spow(params.sinPhi, params.B) *
      GlyphGeomUtility::spow(params.cosTheta, params.A);
  return Point(x, y, z);
}

// Generate around z-axis
Point TensorGlyphBuilder::evaluateSuperquadricPointPlanar(const SuperquadricPointParams& params)
{
  double x, y, z;
  x = GlyphGeomUtility::spow(params.sinPhi, params.B) *
      GlyphGeomUtility::spow(params.cosTheta, params.A);
  y = GlyphGeomUtility::spow(params.sinPhi, params.B) *
      GlyphGeomUtility::spow(params.sinTheta, params.A);
  z = GlyphGeomUtility::spow(params.cosPhi, params.B);
  return Point(x, y, z);
}

void TensorGlyphBuilder::generateBox(GlyphConstructor& constructor)
{
  computeTransforms();

  std::vector<Vector> points = generateBoxPoints();
  std::vector<Vector> normals = rotate_.get_rotation();
  if (flatTensor_)
    for (int d = 0; d < DIMENSIONS_; ++d)
      normals[d] = zeroNorm_;

  generateBoxSide(constructor, {points[5], points[4], points[7], points[6]}, normals[0]);
  generateBoxSide(constructor, {points[7], points[6], points[3], points[2]}, normals[1]);
  generateBoxSide(constructor, {points[1], points[5], points[3], points[7]}, normals[2]);
  generateBoxSide(constructor, {points[3], points[2], points[1], points[0]}, -normals[0]);
  generateBoxSide(constructor, {points[1], points[0], points[5], points[4]}, -normals[1]);
  generateBoxSide(constructor, {points[2], points[6], points[0], points[4]}, -normals[2]);
}

void TensorGlyphBuilder::generateBoxSide(
    GlyphConstructor& constructor, const std::vector<Point>& points, const Vector& normal)
{
  constructor.setOffset();
  for (auto& p : points)
    constructor.addVertex(p, normal, color_);
  constructor.addIndicesToOffset(2, 0, 3);
  constructor.addIndicesToOffset(1, 3, 0);
}

std::vector<Vector> TensorGlyphBuilder::generateBoxPoints()
{
  auto eigvals = t_.getEigenvalues();
  std::vector<Vector> boxPoints;

  for (int x : {-1, 1})
    for (int y : {-1, 1})
      for (int z : {-1, 1})
        boxPoints.emplace_back(trans_ * Point(x * eigvals[0], y * eigvals[1], z * eigvals[2]));

  return boxPoints;
}

Transform TensorGlyphBuilder::getScale()
{
  Transform scale = Transform();
  auto eigvals = t_.getEigenvalues();
  for (int i = 0; i < eigvals.size(); ++i)
    scale.set_mat_val(i, i, eigvals[i]);
  return scale;
}
