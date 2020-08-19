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
#include <chrono>

using namespace SCIRun;
using namespace Graphics;
using namespace Core::Geometry;
using namespace Core::Datatypes;

using GGU = GlyphGeomUtility;

Vector makeSCIRunVector(const Eigen::Vector3d& vec)
{
  return Vector(vec[0], vec[1], vec[2]);
}

UncertaintyTensorOffsetSurfaceBuilder::UncertaintyTensorOffsetSurfaceBuilder(
    const Core::Datatypes::Dyadic3DTensor& t, const Core::Geometry::Point& center, double emphasis)
    : TensorGlyphBuilder(t, center), emphasis_(emphasis)
{}

void UncertaintyTensorOffsetSurfaceBuilder::generateOffsetSurface(
    GlyphConstructor& constructor, const Eigen::Matrix<double, 6, 6>& covarianceMatrix)
{
  auto start = std::chrono::system_clock::now();
  long time = 0;
  const static double h = 0.000001;
  const static double hHalf = 0.5 * h;
  Point origin = Point(0, 0, 0);

  t_.makePositive(true, false);
  // t_.reorderTensorValues();
  computeTransforms();
  postScaleTransforms();
  computeSinCosTable(false);

  auto eigvecs = t_.getEigenvectors();
  auto eigvals = t_.getEigenvalues();
  auto cross = eigvecs[0].cross(eigvecs[1]);
  bool flipVertices = cross.dot(eigvecs[2]) < 2e-12;

  MandelVector tMandel = t_.mandel();
  std::cout << "tMandel " << tMandel << "\n";

  const Transform rotate = Transform(Point(0, 0, 0), makeSCIRunVector(eigvecs[0]),
      makeSCIRunVector(eigvecs[1]), makeSCIRunVector(eigvecs[2]));
  Transform scale = getScale();
  // auto eigvalsVector = Vector(eigvals[0], eigvals[1], eigvals[2]);
  // scale.post_scale(eigvalsVector);
  const Transform scaleThenRotate = rotate * scale;
  Transform rotateInv = rotate;
  rotateInv.invert();

  bool linear = t_.linearCertainty() > t_.planarCertainty();
  auto AAndB = getAAndB(linear, emphasis_);
  SuperquadricPointParams params;
  params.A = AAndB.first;
  params.B = AAndB.second;

  SuperquadricPointParams normalParams;
  normalParams.A = 2.0 - params.A;
  normalParams.B = 2.0 - params.B;

  MandelVector finiteDiff;
  finiteDiff.fill(0.0);

  auto pseudoInv = t_.getEigenvalues();
  Vector pseudoInvVector;
  for (int i = 0; i < 3; ++i)
    pseudoInvVector[i] = 1 / pseudoInv[i];

  int nv = resolution_;
  int nu = resolution_ + 1;
  long avg_time = 0;
  long avg_time_difft = 0;
  for (int v = 0; v < nv - 1; ++v)
  {
    double sinPhi[2] = {tab2_.sin(v), tab2_.sin(v + 1)};
    double cosPhi[2] = {tab2_.cos(v), tab2_.cos(v + 1)};

    for (int u = 0; u < nu; ++u)
    {
      start = std::chrono::system_clock::now();
      constructor.setOffset();
      params.sinTheta = normalParams.sinTheta = tab1_.sin(u);
      params.cosTheta = normalParams.cosTheta = tab1_.cos(u);

      time = std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::system_clock::now() - start)
                 .count();
      // std::cout << "time to calculate angles " << time << "\n";
      for (int i = 0; i < 2; ++i)
      {
        start = std::chrono::system_clock::now();
        params.sinPhi = normalParams.sinPhi = sinPhi[i];
        params.cosPhi = normalParams.cosPhi = cosPhi[i];

        Vector p = scaleThenRotate * Vector(evaluateSuperquadricPoint(linear, params));
        // pseudoInv /= meanTensors_[f].magnitude();

        // Surface Derivative
        Eigen::Vector3d nn;
        nn.fill(0.0);
        Vector diff = Vector(0.0, 0.0, 0.0);
        for (int j = 0; j < 3; ++j)
        {
          diff[j] = hHalf;
          Point newP = Point(pseudoInvVector * (rotateInv * (p + diff)));
          double d1 = evaluateSuperquadricImpl(linear, newP, params.A, params.B);

          newP = Point(pseudoInvVector * (rotateInv * (p - diff)));
          double d2 = evaluateSuperquadricImpl(linear, newP, params.A, params.B);
          diff[j] = 0.0;
          nn(j) = (d1 - d2) / h;
        }

        auto difft_start = std::chrono::system_clock::now();
        MandelVector qn;
        // qn.fill(0.0);
        // Vector pNorm = rotate * Vector(builder.evaluateSuperquadricPoint(linear, params));
        for (int j = 0; j < 6; ++j)
        {
          finiteDiff(j) = hHalf;
          qn(j) = diffT(tMandel + finiteDiff, tMandel - finiteDiff, Point(p));
          // std::cout << "qn " << j << ": " << qn(j) << "\n";
          finiteDiff(j) = 0.0;
          // std::cout << "qn " << j << " " <<qn(j) << "\n";
        }
        auto time_difft = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::system_clock::now() - difft_start)
                              .count();
        avg_time_difft += time_difft;
        // std::cout << "nn norm " << nn.norm() << "\n";
        qn /= h;
        qn /= nn.norm();

        // std::cout << "qn " << qn << "\n";
        // std::cout << "covarianceMatrix " << covarianceMatrix << "\n";
        double q = std::sqrt(
            std::abs((qn.transpose().eval() * (covarianceMatrix * qn).eval()).eval().value()));
        auto n = nn / nn.norm();
        Vector nVector = Vector(n(0), n(1), n(2));

        Vector normal = Vector(evaluateSuperquadricPoint(linear, normalParams));
        normal = rotate * normal;
        normal.safe_normalize();

        Vector offsetP = p + q * normal;

        constructor.addVertex(offsetP + Vector(center_), nVector, ColorRGB(1.0, 1.0, 1.0));
        time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::system_clock::now() - start)
                   .count();
        avg_time += time;
        // std::cout << "time to add vertex " << time << "\n";
      }

      if (flipVertices)
      {
        // std::cout << "vert flip\n";
        constructor.addIndicesToOffset(0, 2, 1);
        constructor.addIndicesToOffset(2, 3, 1);
      }
      else
      {
        // std::cout << "not vert flip\n";
        constructor.addIndicesToOffset(0, 1, 2);
        constructor.addIndicesToOffset(2, 1, 3);
      }
    }
  }
  avg_time /= ((nv-1) * nu * 2);
  avg_time_difft /= ((nv - 1) * nu * 2);
  std::cout << "avg time for new vert " << avg_time << "\n";
  std::cout << "avg time for difft " << avg_time_difft << "\n";
  constructor.popIndicesNTimes(6);
  // avg_time += std::chrono::duration_cast<std::chrono::milliseconds>(
      // std::chrono::system_clock::now() - start)
                  // .count();
}

double UncertaintyTensorOffsetSurfaceBuilder::diffT(
    const MandelVector& t1, const MandelVector& t2, const Point& p)
{
  std::vector<double> dist(2);
  for (int i = 0; i < 2; ++i)
  {
    const MandelVector& tMandel = (i == 0) ? t1 : t2;

    auto start = std::chrono::system_clock::now();
    auto t = symmetricTensorFromMandel(tMandel);
    // std::cout << "t " << t << "\n";
    // t.makePositive(true, false);
    // t.reorderTensorValues();
    bool linear = t.linearCertainty() > t.planarCertainty();
    auto eigvecs = t.getEigenvectors();
    // for (int i = 0 ; i < 3 ; ++i){
      // std::cout << "eigvec " << i << ": ";
      // std::cout << "[" << eigvecs[i][0] << " " << eigvecs[i][1] << " " << eigvecs[i][2] << "]\n";
    // }

    auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::system_clock::now() - start)
                    .count();
    std::cout << "this point " << time << "\n";
    TensorGlyphBuilder builder(t, p);
    builder.computeTransforms();
    builder.postScaleTransforms();
    auto AAndB = builder.getAAndB(linear, emphasis_);
    double A = AAndB.first;
    double B = AAndB.second;
    // Transform trans, rotate, scale;
    // GGU::generateTransforms(p, makeSCIRunVector(eigvecs[0]), makeSCIRunVector(eigvecs[1]), makeSCIRunVector(eigvecs[2]), trans, rotate);
    const Transform rotate = Transform(Point(0,0,0), makeSCIRunVector(eigvecs[0]), makeSCIRunVector(eigvecs[1]), makeSCIRunVector(eigvecs[2]));
    const Transform scale = builder.getScale();
    // const Transform rotate = Transform(Point(0, 0, 0), makeSCIRunVector(eigvecs[0]),
        // makeSCIRunVector(eigvecs[1]), makeSCIRunVector(eigvecs[2]));
    // Transform scale;
    auto eigvals = t.getEigenvalues();
    auto eigvalsVector = Vector(eigvals[0], eigvals[1], eigvals[2]);
    // for (int i = 0; i < eigvals.size(); ++i)
      // scale.set_mat_val(i, i, eigvals[i]);
    // scale.post_scale(eigvalsVector);
    auto rotateInv = rotate;
    rotateInv.invert();
    // std::cout << "rotateInv " << rotateInv << "\n";
    auto scaleInv = scale;
    scaleInv.invert();
    // std::cout << "scaleInv " << scaleInv << "\n";

    // auto t = newBuilder.getTensor();
    Vector pseudoInv = eigvalsVector;
    for (int i = 0; i < 3; ++i)
      pseudoInv[i] = 1 / pseudoInv[i];
    // std::cout << "pseudoInv " << pseudoInv << "\n";

    Vector newP = (pseudoInv * Vector(rotateInv * p));
    // std::cout << "newP " << newP << "\n";
    dist[i] = evaluateSuperquadricImpl(linear, Point(newP), A, B);
    // std::cout << "dist " << i << " " << dist[i] <<"\n";
  }

  // std::cout << "dist diff " << dist[0] - dist[1] << "\n";
  return dist[0] - dist[1];
}

double UncertaintyTensorOffsetSurfaceBuilder::evaluateSuperquadricImpl(
    bool linear, const Point& p, double A, double B)
{
  if (linear)
    return evaluateSuperquadricImplLinear(p, A, B);
  else
    return evaluateSuperquadricImplPlanar(p, A, B);
}

// Generate around x-axis
double UncertaintyTensorOffsetSurfaceBuilder::evaluateSuperquadricImplLinear(
    const Point& p, double A, double B)
{
  double twoDivA = 2.0 / A;
  double twoDivB = 2.0 / B;
  return GGU::spow(
             GGU::spow(std::abs(p.y()), twoDivA) + GGU::spow(std::abs(p.z()), twoDivA), A / B) +
         GGU::spow(std::abs(p.x()), twoDivB) - 1;
}

// Generate around z-axis
double UncertaintyTensorOffsetSurfaceBuilder::evaluateSuperquadricImplPlanar(
    const Point& p, double A, double B)
{
  double twoDivA = 2.0 / A;
  double twoDivB = 2.0 / B;
  return GGU::spow(
             GGU::spow(std::abs(p.x()), twoDivA) + GGU::spow(std::abs(p.y()), twoDivA), A / B) +
         GGU::spow(std::abs(p.z()), twoDivB) - 1;
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

void TensorGlyphBuilder::setResolution(int resolution)
{
  resolution_ = resolution;
}

void TensorGlyphBuilder::makeTensorPositive(bool reorder, bool makeGlyph)
{
  t_.makePositive(reorder, makeGlyph);

  auto eigvals = t_.getEigenvalues();
  // This is exactly zero after thresholding
  flatTensor_ = t_.getEigenvalue(2) == 0;
  zeroNorm_ = makeSCIRunVector(t_.getEigenvector(2));
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

void TensorGlyphBuilder::postScaleTransforms()
{
  auto eigvals = t_.getEigenvalues();
  Vector eigvalsVector(eigvals[0], eigvals[1], eigvals[2]);

  trans_.post_scale(eigvalsVector);
}

void TensorGlyphBuilder::generateEllipsoid(GlyphConstructor& constructor, bool half)
{
  makeTensorPositive(true);
  t_.reorderTensorValues();
  computeTransforms();
  postScaleTransforms();
  computeSinCosTable(half);

  auto eigvals = t_.getEigenvalues();
  Vector eigvalsVector(eigvals[0], eigvals[1], eigvals[2]);
  Transform rotateThenInvScale = rotate_;
  rotateThenInvScale.post_scale(Vector(1.0, 1.0, 1.0) / eigvalsVector);

  for (int v = 0; v < nv_ - 1; ++v)
  {
    double sinPhi[2] = {tab2_.sin(v + 1), tab2_.sin(v)};
    double cosPhi[2] = {tab2_.cos(v + 1), tab2_.cos(v)};

    for (int u = 0; u < nu_; ++u)
    {
      EllipsoidPointParams params;
      params.sinTheta = tab1_.sin(u);
      params.cosTheta = tab1_.cos(u);

      // Transform points and add to points list
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

std::pair<double, double> TensorGlyphBuilder::getAAndB(bool linear, double emphasis)
{
  double pPower = GlyphGeomUtility::spow((1.0 - t_.planarCertainty()), emphasis);
  double lPower = GlyphGeomUtility::spow((1.0 - t_.linearCertainty()), emphasis);
  if (linear)
    return std::make_pair(pPower, lPower);
  else
    return std::make_pair(lPower, pPower);
}

void TensorGlyphBuilder::generateSuperquadricTensor(GlyphConstructor& constructor, double emphasis)
{
  makeTensorPositive(true);
  t_.reorderTensorValues();
  computeTransforms();
  postScaleTransforms();
  computeSinCosTable(false);

  bool linear = t_.linearCertainty() > t_.planarCertainty();
  auto AAndB = getAAndB(linear, emphasis);
  SuperquadricPointParams params;
  params.A = AAndB.first;
  params.B = AAndB.second;

  SuperquadricPointParams normalParams;
  normalParams.A = 2.0 - params.A;
  normalParams.B = 2.0 - params.B;
  auto eigvecs = t_.getEigenvectors();
  auto cross = eigvecs[0].cross(eigvecs[1]);
  bool flipVertices = cross.dot(eigvecs[2]) < 2e-12;

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
        // Transform points and add to points list
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

      if (flipVertices)
      {
        constructor.addIndicesToOffset(0, 2, 1);
        constructor.addIndicesToOffset(2, 3, 1);
      }
      else
      {
        constructor.addIndicesToOffset(0, 1, 2);
        constructor.addIndicesToOffset(2, 1, 3);
      }
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
  makeTensorPositive(true);
  computeTransforms();

  std::vector<Vector> points = generateBoxPoints();
  std::vector<Vector> normals = rotate_.get_rotation();
  // std::cout << "normals\n";
  for (auto& v : normals)
    // std::cout << v << "\n";
  if (flatTensor_) {
    // std::cout << "is flat!\n";
    for (int d = 0; d < DIMENSIONS_; ++d)
      normals[d] = zeroNorm_;
  }

  generateBoxSide(constructor, {points[5], points[4], points[7], points[6]}, normals[0]);
  generateBoxSide(constructor, {points[7], points[6], points[3], points[2]}, normals[1]);
  generateBoxSide(constructor, {points[1], points[5], points[3], points[7]}, normals[2]);
  generateBoxSide(constructor, {points[3], points[2], points[1], points[0]}, -normals[0]);
  generateBoxSide(constructor, {points[1], points[0], points[5], points[4]}, -normals[1]);
  generateBoxSide(constructor, {points[2], points[6], points[0], points[4]}, -normals[2]);
}

void TensorGlyphBuilder::generateBoxSide(
    GlyphConstructor& constructor, const std::vector<Vector>& points, const Vector& normal)
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
