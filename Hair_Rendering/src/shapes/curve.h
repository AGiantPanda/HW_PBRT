
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_CURVE_H
#define PBRT_SHAPES_CURVE_H

// shapes/curve.h*
// panda:
// modeled with cubic Bezier splines, defined by four control points p0, p1, p2, p3
// intermediate points are: 
// p(u) = (1 − u)^3 * p0 + 3 * (1 − u)^2 * u * p1 + 3 * (1 − u) * u^2 * p2 + u^3 * p3
#include "shape.h"

struct CurveCommon;

// CurveType Declarations
// panda:
// Flat: oriented to face the ray being interesected, hair or fur
// Cylinder: a true cylinder shape such as spaghetti
// Ribbon: shapes that don't actually have a cylindrical cross section, blade of grass
enum CurveType { Flat, Cylinder, Ribbon };

// CurveCommon Declarations
// panda:
// stores the control points and other info that is shared across curve segments
struct CurveCommon {
    CurveCommon(const Point c[4], float w0, float w1, CurveType type,
                const Normal *norm);
    const CurveType type;
    Point cpObj[4];
    float width[2];
    Normal n[2];
    float normalAngle, invSinNormalAngle;
};

// Curve Declarations
class Curve : public Shape {
  public:
    // Curve Public Methods
    Curve(const Transform *ObjectToWorld, const Transform *WorldToObject,
          bool reverseOrientation, const CurveCommon* common,
          float uMin, float uMax)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
          common(common),
          uMin(uMin),
          uMax(uMax) {}
    BBox ObjectBound() const;
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
    float Area() const;
    Point Sample(float u1, float u2, Normal *Ns) const;

  private:
    // Curve Private Methods
    bool recursiveIntersect(const Ray &r, float *tHit,
                            float *rayEpsilon, DifferentialGeometry *dg, const Point cp[4],
                            const Transform &rayToObject, float u0, float u1,
                            int depth) const;

    // Curve Private Data
    const CurveCommon* common;
    const float uMin, uMax;
};

std::vector<Reference<Shape> > CreateCurveShape(const Transform *o2w,
                        const Transform *w2o,
                        bool reverseOrientation,
                        const ParamSet &params);


#endif  // PBRT_SHAPES_CURVE_H
