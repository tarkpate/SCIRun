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
#include <ios>
#include "Core/Datatypes/TensorFwd.h"

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
{
  h_ = 0.000001;// * t.frobeniusNorm();
  hHalf_ = h_ * 0.5;
}

/* TODO remove. this is for testing purposes
double UncertaintyTensorOffsetSurfaceBuilder::diffT(const MandelVector& s1,
                                                    const MandelVector& s2,
                                                    const Eigen::Vector3d& p, double emphasis)
{

  auto t1 = symmetricTensorFromMandel(s1);
  auto t2 = symmetricTensorFromMandel(s2);

  // std::cout << "T1: " << t1 <<"\n";
  // std::cout << "T2: " << t2 <<"\n";
  // std::cout << "t1 eigvals: " << t1.getEigenvalues() <<"\n";
  // std::cout << "t2 eigvals: " << t2.getEigenvalues() <<"\n";

  auto cl1 = t1.linearCertainty();
  auto cp1 = t1.planarCertainty();
  // std::cout << "cl1: " << std::setprecision(10) << cl1 <<"\n";
  // std::cout << "cp1: " << std::setprecision(10) << cp1 <<"\n";

  auto cl2 = t2.linearCertainty();
  auto cp2 = t2.planarCertainty();
  // std::cout << "cl2: " << std::setprecision(10) << cl2 <<"\n";
  // std::cout << "cp2: " << std::setprecision(10) << cp2 <<"\n";

  double A1, B1, A2, B2;

  bool linear1 = (cl1 > cp1) || ((std::abs(cl1) - std::abs(cp1) > 0.00001));
  bool linear2 = (cl2 > cp2) || ((std::abs(cl2) - std::abs(cp2) > 0.00001));
  if(linear1)  {
    A1 = GGU::spow((1.0-cp1), emphasis);
    B1 = GGU::spow((1.0-cl1), emphasis);
  }
  else
  {
    A1 = GGU::spow((1.0-cl1), emphasis);
    B1 = GGU::spow((1.0-cp1), emphasis);
  }
  if(linear2)
  {
    A2 = GGU::spow((1.0-cp2), emphasis);
    B2 = GGU::spow((1.0-cl2), emphasis);
  }
  else
  {
    A2 = GGU::spow((1.0-cl2), emphasis);
    B2 = GGU::spow((1.0-cp2), emphasis);
  }


  // std::cout << std::scientific;
  const Eigen::Vector3d newP1 = t1.getEigenvalues().asDiagonal().inverse() * t1.getEigenvectorsAsMatrix().transpose() * p;
  const Eigen::Vector3d newP2 = t2.getEigenvalues().asDiagonal().inverse() * t2.getEigenvectorsAsMatrix().transpose() * p;
  // std::cout << "newP1: " << newP1 << "\n";
  // std::cout << "A1: " << A1 << "\n";
  // std::cout << "B1: " << B1 << "\n";
  // std::cout << "newP2: " << newP2 << "\n";
  // std::cout << "newP2.x(): " << newP2.x() << "\n";
  // std::cout << "abs(newP2.x()): " << std::abs(newP2.x()) << "\n";
  // std::cout << "A2: " << A2 << "\n";
  // std::cout << "B2: " << B2 << "\n";

  auto g1 = evaluateSuperquadricImpl(linear1, newP1, A1, B1);
  auto g2 = evaluateSuperquadricImpl(linear2, newP2, A2, B2);

  // std::cout << "newP2: " << newP2 << "\n";
  // std::cout << "newP2.x(): " << newP2.x() << "\n";
  // std::cout << "abs(newP2.x()): " << std::abs(newP2.x()) << "\n";
  // auto g2 = GGU::spow(GGU::spow(abs(newP2.x()), 2.0/A2) + GGU::spow(abs(newP2.y()), 2.0/A2), A2/B2)
    // + GGU::spow(abs(newP2.z()), 2.0/B2) - 1.0;

  // std::cout << "newP2: " << newP2 << "\n";
  // std::cout << "newP2.x(): " << newP2.x() << "\n";
  // std::cout << "abs(newP2.x()): " << std::abs(newP2.x()) << "\n";
  // std::cout << "TEST p1 with newP2: " << std::pow(std::abs(newP2.x()), 2.0 / A2) << "\n";
  // std::cout << "TEST p1 with 0: " << std::pow(std::abs(0), 2.0 / A2) << "\n";
  // std::cout << "TEST p1 1: " << std::pow(std::abs(newP2.y()), 2.0 / A2) << "\n";
  // std::cout << "A2/B2: " << A2/B2 << "\n";
  // std::cout << "TEST p1: "
            // << std::pow(std::pow(std::abs(newP2.x()), 2.0 / A2) +
                            // std::pow(std::abs(newP2.y()), 2.0 / A2),
                   // A2 / B2)
            // << "\n";
  // std::cout << "TEST p2: " << std::pow(std::abs(newP2.z()), 2 / B2) << "\n";

  // std::cout << "TEST: " << std::pow(std::pow(std::abs(newP2.x()), 2.0/A2) + std::pow(std::abs(newP2.y()), 2.0/A2), A2/B2) + std::pow(std::abs(newP2.z()), 2/B2) - 1.0 << "\n";
  // std::cout << "g1: " << g1 << "\n";
  // std::cout << "g2: " << g2 << "\n";
  // std::cout << "g1-g2: " << g1-g2 << "\n";

  return g1-g2;
}
*/

void UncertaintyTensorOffsetSurfaceBuilder::generateOffsetSurface(
    GlyphConstructor& constructor, const Eigen::Matrix<double, 6, 6>& covarianceMatrix)
{
  // auto start = std::chrono::system_clock::now();
  // long time = 0;
  // Point origin = Point(0, 0, 0);
  Vector centerVector = Vector(center_);

  t_.makePositive(true, false);
  computeTransforms();
  postScaleTransforms();
  computeSinCosTable(false);

  std::vector<Eigen::Vector3d> eigvecs = t_.getEigenvectors();
  // std::cout << "V: " << t_.getEigenvectorsAsMatrix() << "\n";
  // std::cout << "V[0]: " << t_.getEigenvector(0) << "\n";
  // auto cross = eigvecs[0].cross(eigvecs[1]);

  MandelVector tMandel = t_.mandel();

  Eigen::Vector3d a(1,2,3);
  Eigen::Vector3d b = a.asDiagonal().inverse().diagonal();

  // const Transform rotate = Transform(Point(0, 0, 0), makeSCIRunVector(eigvecs[0]),
      // makeSCIRunVector(eigvecs[1]), makeSCIRunVector(eigvecs[2]));
  Eigen::Matrix3d rotate = t_.getEigenvectorsAsMatrix();
  // std::cout << "rotate: " << rotate << "\n";
  Eigen::Matrix3d rotateTranspose = rotate.transpose();
  // Transform scale = getScale();
  Eigen::Vector3d eigvals = t_.getEigenvalues();
  Eigen::Vector3d eigvalsInv = eigvals.asDiagonal().inverse().diagonal();
  Eigen::Matrix3d scale = eigvals.asDiagonal();
  Eigen::Matrix3d scaleInv = eigvalsInv.asDiagonal();
  // std::cout << "scale: " << scale << "\n";

  const auto scaleThenRotate = rotate * scale;
  const auto undoScaleAndRotate = scaleInv * rotateTranspose;

  double cl = t_.linearCertainty();
  double cp = t_.planarCertainty();
  // std::cout << "cl: " << cl <<"\n";
  // std::cout << "cp: " << cp <<"\n";
  bool linear = cl >= cp;
  linear = false;
  std::cout << "linear: " << linear << "\n";
  auto AAndB = getAAndB(cl, cp, linear, emphasis_);
  SuperquadricPointParams params;
  params.A = AAndB.first;
  params.B = AAndB.second;

  SuperquadricPointParams normalParams;
  normalParams.A = 2.0 - params.A;
  normalParams.B = 2.0 - params.B;

  // auto pseudoInv = t_.getEigenvalues();
  // Vector pseudoInvVector;
  MandelVector qn;
  // for (int i = 0; i < 3; ++i)
    // pseudoInvVector[i] = 1.0 / pseudoInv[i];

  DifftValues difftVals;
  precalculateDifftValues(difftVals, tMandel);
  int nv = resolution_;
  int nu = resolution_ + 1;

  for (int v = 0; v < nv - 1; ++v)
  {
    double sinPhi[2] = {tab2_.sin(v), tab2_.sin(v + 1)};
    double cosPhi[2] = {tab2_.cos(v), tab2_.cos(v + 1)};
    for (int u = 0; u < nu; ++u)
    {
      constructor.setOffset();
      params.sinTheta = normalParams.sinTheta = tab1_.sin(u);
      params.cosTheta = normalParams.cosTheta = tab1_.cos(u);
      for (int i = 0; i < 2; ++i)
      {
        params.sinPhi = normalParams.sinPhi = sinPhi[i];
        params.cosPhi = normalParams.cosPhi = cosPhi[i];

        // std::cout << "sinPhi: " << params.sinPhi << "\n";
        // std::cout << "cosPhi: " << params.cosPhi << "\n";
        // std::cout << "sinTheta: " << params.sinTheta << "\n";
        // std::cout << "cosTheta: " << params.cosTheta << "\n";
        // std::cout << "eval sup: " << evaluateSuperquadricPoint(linear, params) << "\n";
        // std::cout << "scale*eval: " << scale*evaluateSuperquadricPoint(linear, params) << "\n";
        Vector pVector = evaluateSuperquadricPoint(linear, params);
        Eigen::Vector3d p = scaleThenRotate * Eigen::Vector3d(pVector.x(), pVector.y(), pVector.z());
        // std::cout << "p: " << p << "\n";
        Vector normalVector = evaluateSuperquadricPoint(linear, normalParams);
        Eigen::Vector3d normal = Eigen::Vector3d(normalVector.x(), normalVector.y(), normalVector.z());

        // Surface Derivative
        Eigen::Vector3d nn;
        Eigen::Vector3d diff(0.0, 0.0, 0.0);
        for (int j = 0; j < 3; ++j)
        {
          diff[j] = hHalf_;
          Eigen::Vector3d newP = undoScaleAndRotate * (p + diff);
          double d1 = evaluateSuperquadricImpl(linear, newP, params.A, params.B);

          newP = undoScaleAndRotate * (p - diff);
          double d2 = evaluateSuperquadricImpl(linear, newP, params.A, params.B);
          diff[j] = 0.0;
          nn(j) = (d1 - d2) / h_;
        }
        // std::cout << "nn: " << nn << "\n";
        auto qn = getQn(difftVals, p);
        // std::cout << "nn.norm: " << nn.norm() << "\n";

        // MandelVector diffs;
        // MandelVector difference;
        // difference << 0, 0, 0, 0, 0, 0;
        // for (int i = 0; i < 6; ++i)
        // {
          // difference(i) = hHalf_;
          // diffs(i) = diffT(tMandel + difference, tMandel - difference, p, emphasis_);
          // std::cout << "diffs ret: " << diffs(i) << "\n";
          // difference(i) = 0;
        // }

        // std::cout << "diffs before h: " << diffs.transpose() << "\n";
        // diffs /= h_;
        qn /= h_;
        qn /= nn.norm();
        // std::cout << "qn: " << qn << "\n";
        // std::cout << "diffs before nn norm: " << diffs.transpose() << "\n";
        // diffs /= nn.norm();
        // std::cout << "qn: " << diffs.transpose() << "\n";
        double q = std::sqrt(
            std::abs((qn.transpose().eval() * (covarianceMatrix * qn).eval()).eval().value()));
        // std::cout << "q: " << q << "\n\n";
        Eigen::Vector3d n = nn / nn.norm();
        Eigen::Vector3d rotatedN = n;
        Vector rotatedNVector = Vector(rotatedN.x(), rotatedN.y(), rotatedN.z());

        Eigen::Vector3d rotatedNormal = rotate * normal;

        Eigen::Vector3d offsetP = p + q * rotatedNormal;
        Vector offsetPVector = Vector(offsetP.x(), offsetP.y(), offsetP.z());

        // if (secondHalf) rotatedNVector *= -1;
        constructor.addVertex(offsetPVector + centerVector, rotatedNVector, ColorRGB(1.0, 1.0, 1.0));
      }

      constructor.addIndicesToOffset(0, 1, 2);
      constructor.addIndicesToOffset(2, 1, 3);
    }
  }
  constructor.popIndicesNTimes(6);
}

void UncertaintyTensorOffsetSurfaceBuilder::precalculateDifftValues(
    DifftValues& vals, const MandelVector& t)
{
  MandelVector finiteDiff;
  finiteDiff.fill(0.0);

  for (auto i = 0; i < 6; ++i)
  {
    finiteDiff(i) = hHalf_;
    std::vector<double> dist(2);
    for (auto s = 0; s < 2; ++s)
    {
      auto index = 2 * i + s;
      const MandelVector& tMandel =
          (s == 0) ? MandelVector(t + finiteDiff) : MandelVector(t - finiteDiff);

      auto newT = symmetricTensorFromMandel(tMandel);
      double cl = newT.linearCertainty();
      double cp = newT.planarCertainty();
      vals.linear[index] = cl >= cp;
      auto AAndB = getAAndB(cl, cp, vals.linear[index], emphasis_);
      vals.A[index] = AAndB.first;
      vals.B[index] = AAndB.second;
      Eigen::Matrix3d rotateTranspose = newT.getEigenvectorsAsMatrix().transpose();
      Eigen::Matrix3d scaleInv = newT.getEigenvalues().asDiagonal().inverse();
      vals.undoScaleAndRotate[index] = scaleInv * rotateTranspose;
      vals.normEigvals[index] = newT.getNormalizedEigenvalues();
    }
    finiteDiff(i) = 0.0;
  }
}

MandelVector UncertaintyTensorOffsetSurfaceBuilder::getQn(
    const DifftValues& vals, const Eigen::Vector3d& p)
{
  MandelVector qn;
  for (int i = 0; i < 6; ++i)
  {
    std::vector<double> dist(2);
    for (int s = 0; s < 2; ++s)
    {
      auto index = 2 * i + s;
      Eigen::Vector3d newP = vals.undoScaleAndRotate[index] * p;
      dist[s] = evaluateSuperquadricImpl(vals.linear[index], newP, vals.A[index], vals.B[index]);
    }
    qn(i) = dist[0] - dist[1];
  }

  return qn;
}

double UncertaintyTensorOffsetSurfaceBuilder::evaluateSuperquadricImpl(
    bool linear, const Eigen::Vector3d& p, double A, double B)
{
  if (linear)
    return evaluateSuperquadricImplLinear(p, A, B);
  else
    return evaluateSuperquadricImplPlanar(p, A, B);
}

// Generate around x-axis
double UncertaintyTensorOffsetSurfaceBuilder::evaluateSuperquadricImplLinear(
    const Eigen::Vector3d& p, double A, double B)
{
  // auto orthoT = Dyadic3DTensor({Eigen::Vector3d(1,0,0),Eigen::Vector3d(0,1,0),Eigen::Vector3d(0,0,1)}, p);
  // auto cl = orthoT.linearCertainty();
  // auto cp = orthoT.planarCertainty();
  // auto AAndB = getAAndB(cl, cp, cl > cp, emphasis_);
  // A = AAndB.first;
  // B = AAndB.second;
  const double twoDivA = 2.0 / A;
  const double twoDivB = 2.0 / B;
  return GGU::spow(
             GGU::spow(std::abs(p.y()), twoDivA) + GGU::spow(std::abs(p.z()), twoDivA), A / B) +
         GGU::spow(std::abs(p.x()), twoDivB) - 1;
}

// Generate around z-axis
double UncertaintyTensorOffsetSurfaceBuilder::evaluateSuperquadricImplPlanar(
    const Eigen::Vector3d& p, double A, double B)
{
  const double twoDivA = 2.0 / A;
  const double twoDivB = 2.0 / B;
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

std::pair<double, double> TensorGlyphBuilder::getAAndB(
    double cl, double cp, bool linear, double emphasis)
{
  double pPower = GlyphGeomUtility::spow(1.0 - cp, emphasis);
  double lPower = GlyphGeomUtility::spow(1.0 - cl, emphasis);
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

  double cl = t_.linearCertainty();
  double cp = t_.planarCertainty();
  bool linear = cl >= cp;
  auto AAndB = getAAndB(cl, cp, linear, emphasis);
  SuperquadricPointParams params;
  params.A = AAndB.first;
  params.B = AAndB.second;

  SuperquadricPointParams normalParams;
  normalParams.A = 2.0 - params.A;
  normalParams.B = 2.0 - params.B;
  auto eigvecs = t_.getEigenvectors();
  auto cross = eigvecs[0].cross(eigvecs[1]);

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
        Point p = Point(evaluateSuperquadricPoint(linear, params));
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

Vector TensorGlyphBuilder::evaluateSuperquadricPoint(
    bool linear, const SuperquadricPointParams& params)
{
  if (linear)
    return evaluateSuperquadricPointLinear(params);
  else
    return evaluateSuperquadricPointPlanar(params);
}

// Generate around x-axis
Vector TensorGlyphBuilder::evaluateSuperquadricPointLinear(const SuperquadricPointParams& params)
{
  double x, y, z;
  x = GlyphGeomUtility::spow(params.cosPhi, params.B);
  y = -GlyphGeomUtility::spow(params.sinPhi, params.B) *
      GlyphGeomUtility::spow(params.sinTheta, params.A);
  z = GlyphGeomUtility::spow(params.sinPhi, params.B) *
      GlyphGeomUtility::spow(params.cosTheta, params.A);
  return Vector(x, y, z);
}

// Generate around z-axis
Vector TensorGlyphBuilder::evaluateSuperquadricPointPlanar(const SuperquadricPointParams& params)
{
  double x, y, z;
  x = GlyphGeomUtility::spow(params.sinPhi, params.B) *
      GlyphGeomUtility::spow(params.cosTheta, params.A);
  y = GlyphGeomUtility::spow(params.sinPhi, params.B) *
      GlyphGeomUtility::spow(params.sinTheta, params.A);
  z = GlyphGeomUtility::spow(params.cosPhi, params.B);
  return Vector(x, y, z);
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
    if (flatTensor_)
    {
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
