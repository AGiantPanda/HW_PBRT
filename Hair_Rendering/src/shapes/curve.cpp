
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

// shapes/curve.cpp*
#include "stdafx.h"
#define SHAPE_DEBUG
#ifdef SHAPE_DEBUG
#include "vdb-win/vdb.h"
#endif
#include "shapes/curve.h"
#include "paramset.h"
#include <cmath>
//#include "stats.h"

//STAT_MEMORY_COUNTER("Memory/Curves", curveBytes);
//STAT_PERCENT("Intersections/Ray-curve intersection tests", nHits, nTests);
//STAT_INT_DISTRIBUTION("Intersections/Curve refinement level", refinementLevel);
//STAT_COUNTER("Scene/Curves", nCurves);
//STAT_COUNTER("Scene/Split curves", nSplitCurves);

// Curve Utility Functions
static Point BlossomBezier(const Point p[4], float u0, float u1, float u2) {
    //Point a[3] = {LerpP(u0, p[0], p[1]), LerpP(u0, p[1], p[2]),
    //                LerpP(u0, p[2], p[3])};
    //Point b[2] = {LerpP(u1, a[0], a[1]), LerpP(u1, a[1], a[2])};
	Point a[3] = {
		(1 - u0) * p[0] + u0 * p[1],
		(1 - u0) * p[1] + u0 * p[2],
		(1 - u0) * p[2] + u0 * p[3]
	};
	Point b[2] = {
		(1 - u1) * a[0] + u1 * a[1],
		(1 - u1) * a[1] + u1 * a[2]
	};
	return (1 - u2) * b[0] + u2 * b[1];// LerpP(u2, b[0], b[1]);
}

inline void SubdivideBezier(const Point cp[4], Point cpSplit[7]) {
    cpSplit[0] = cp[0];
    cpSplit[1] = (cp[0] + cp[1]) / 2;
    cpSplit[2] = (cp[0] + 2 * cp[1] + cp[2]) / 4;
    cpSplit[3] = (cp[0] + 3 * cp[1] + 3 * cp[2] + cp[3]) / 8;
    cpSplit[4] = (cp[1] + 2 * cp[2] + cp[3]) / 4;
    cpSplit[5] = (cp[2] + cp[3]) / 2;
    cpSplit[6] = cp[3];
}

static Point EvalBezier(const Point cp[4], float u,
                          Vector *deriv = NULL) {
    //Point cp1[3] = {LerpP(u, cp[0], cp[1]), LerpP(u, cp[1], cp[2]),
    //                  LerpP(u, cp[2], cp[3])};
    //Point cp2[2] = {LerpP(u, cp1[0], cp1[1]), LerpP(u, cp1[1], cp1[2])};
	Point cp1[3] = {
		(1 - u) * cp[0] + u * cp[1],
		(1 - u) * cp[1] + u * cp[2],
		(1 - u) * cp[2] + u * cp[3]
	};
	Point cp2[2] = {
		(1 - u) * cp1[0] + u * cp1[1],
		(1 - u) * cp1[1] + u * cp1[2]
	};
    if (deriv) *deriv = (float)3 * (cp2[1] - cp2[0]);
	return (1 - u) * cp2[0] + u * cp2[1];// LerpP(u, cp2[0], cp2[1]);
}

// Curve Method Definitions
CurveCommon::CurveCommon(const Point c[4], float width0, float width1,
                         CurveType type, const Normal *norm)
    : type(type) {
    width[0] = width0;
    width[1] = width1;
    for (int i = 0; i < 4; ++i)
        cpObj[i] = c[i];
    if (norm) {
        n[0] = Normalize(norm[0]);
        n[1] = Normalize(norm[1]);
        normalAngle = std::acos(Clamp(Dot(n[0], n[1]), 0.0, 1.0));
        invSinNormalAngle = 1 / std::sin(normalAngle);
    }
    //++nCurves;
}

std::vector<Reference<Shape> > CreateCurve(
    const Transform *o2w, const Transform *w2o, bool reverseOrientation,
    const Point c[4], float w0, float w1, CurveType type,
    const Normal *norm, int splitDepth) {
    std::vector<Reference<Shape> > segments;
    CurveCommon* common =
        new CurveCommon(c, w0, w1, type, norm);
    const int nSegments = 1 << splitDepth;
    segments.reserve(nSegments);
    for (int i = 0; i < nSegments; ++i) {
        float uMin = i / (float)nSegments;
        float uMax = (i + 1) / (float)nSegments;
        segments.push_back(new Curve(o2w, w2o, reverseOrientation,
                                                   common, uMin, uMax));
        //++nSplitCurves;
    }
    //curveBytes += sizeof(CurveCommon) + nSegments * sizeof(Curve);
    return segments;
}

BBox Curve::ObjectBound() const {
    // Compute object-space control points for curve segment, _cpObj_
    Point cpObj[4];
    // panda:
    // the blossom p(u, u, u) gives the curve's value at position u
    cpObj[0] = BlossomBezier(common->cpObj, uMin, uMin, uMin);
    cpObj[1] = BlossomBezier(common->cpObj, uMin, uMin, uMax);
    cpObj[2] = BlossomBezier(common->cpObj, uMin, uMax, uMax);
    cpObj[3] = BlossomBezier(common->cpObj, uMax, uMax, uMax);
    BBox b =
        Union(BBox(cpObj[0], cpObj[1]), BBox(cpObj[2], cpObj[3]));
    float width[2] = {Lerp(uMin, common->width[0], common->width[1]),
                      Lerp(uMax, common->width[0], common->width[1])};
	float delta = std::max(width[0], width[1]) * 0.5f;
	b.Expand(delta);
	return b;
}

bool Curve::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
                      DifferentialGeometry *dg) const {
    //ProfilePhase p(isect ? Prof::CurveIntersect : Prof::CurveIntersectP);
    //++nTests;
    // Transform _Ray_ to object space
    Vector oErr, dErr;
	Ray ray;
	(*WorldToObject)(r, &ray);

    // Compute object-space control points for curve segment, _cpObj_
    Point cpObj[4];
    cpObj[0] = BlossomBezier(common->cpObj, uMin, uMin, uMin);
    cpObj[1] = BlossomBezier(common->cpObj, uMin, uMin, uMax);
    cpObj[2] = BlossomBezier(common->cpObj, uMin, uMax, uMax);
    cpObj[3] = BlossomBezier(common->cpObj, uMax, uMax, uMax);

    // Project curve control points to plane perpendicular to ray

    // Be careful to set the "up" direction passed to LookAt() to equal the
    // vector from the first to the last control points.  In turn, this
    // helps orient the curve to be roughly parallel to the x axis in the
    // ray coordinate system.
    //
    // In turn (especially for curves that are approaching stright lines),
    // we get curve bounds with minimal extent in y, which in turn lets us
    // early out more quickly in recursiveIntersect().
    Vector dx = Cross(ray.d, cpObj[3] - cpObj[0]);
    if (dx.LengthSquared() == 0) {
        // If the ray and the vector between the first and last control
        // points are parallel, dx will be zero.  Generate an arbitrary xy
        // orientation for the ray coordinate system so that intersection
        // tests can proceeed in this unusual case.
        Vector dy;
        CoordinateSystem(ray.d, &dx, &dy);
    }

    // panda:
    // create a ray-based coordinate system
    // reduce the number of operations
    Transform objectToRay = LookAt(ray.o, ray.o + ray.d, dx);
    Point cp[4] = {objectToRay(cpObj[0]), objectToRay(cpObj[1]),
                     objectToRay(cpObj[2]), objectToRay(cpObj[3])};

    // Before going any further, see if the ray's bounding box intersects
    // the curve's bounding box. We start with the y dimension, since the y
    // extent is generally the smallest (and is often tiny) due to our
    // careful orientation of the ray coordinate ysstem above.
    float maxWidth = std::max(Lerp(uMin, common->width[0], common->width[1]),
                              Lerp(uMax, common->width[0], common->width[1]));
    if (std::max(std::max(cp[0].y, cp[1].y), std::max(cp[2].y, cp[3].y)) +
            0.5f * maxWidth < 0 ||
        std::min(std::min(cp[0].y, cp[1].y), std::min(cp[2].y, cp[3].y)) -
            0.5f * maxWidth > 0)
        return false;

    // Check for non-overlap in x.
    if (std::max(std::max(cp[0].x, cp[1].x), std::max(cp[2].x, cp[3].x)) +
            0.5f * maxWidth < 0 ||
        std::min(std::min(cp[0].x, cp[1].x), std::min(cp[2].x, cp[3].x)) -
            0.5f * maxWidth > 0)
        return false;

    // Check for non-overlap in z.
    float rayLength = ray.d.Length();
    float zMax = rayLength * ray.maxt;
    if (std::max(std::max(cp[0].z, cp[1].z), std::max(cp[2].z, cp[3].z)) +
            0.5f * maxWidth < 0 ||
        std::min(std::min(cp[0].z, cp[1].z), std::min(cp[2].z, cp[3].z)) -
            0.5f * maxWidth > zMax)
        return false;

    // Compute refinement depth for curve, _maxDepth_
    float L0 = 0;
    for (int i = 0; i < 2; ++i)
        L0 = std::max(
            L0, std::max(
                    std::max(std::abs(cp[i].x - 2 * cp[i + 1].x + cp[i + 2].x),
                             std::abs(cp[i].y - 2 * cp[i + 1].y + cp[i + 2].y)),
                    std::abs(cp[i].z - 2 * cp[i + 1].z + cp[i + 2].z)));

    float eps =
        std::max(common->width[0], common->width[1]) * .05f;  // width / 20

    int r0;
    float v = 1.41421356237f * 6.f * L0 / (8.f * eps);
    if(v < 1)
        r0 = 0;
    else
    {
        uint32_t bits;
        memcpy(&bits, &v, sizeof(float));
        r0 = (bits >> 23) - 127 + (bits & (1 << 22) ? 1 : 0);
        r0 /= 2;
    }

    int maxDepth = Clamp(r0, 0, 10);
    //ReportValue(refinementLevel, maxDepth);

    return recursiveIntersect(ray, tHit, rayEpsilon, dg, cp, Inverse(objectToRay), uMin,
                              uMax, maxDepth);
}

// panda:
// test whether the given ray intersects the given curve segment over the given parametric range [u0, u1]
bool Curve::recursiveIntersect(const Ray &ray, float *tHit,
                               float *rayEpsilon, DifferentialGeometry *dg, const Point cp[4],
                               const Transform &rayToObject, float u0, float u1,
                               int depth) const {
    float rayLength = ray.d.Length();

    if (depth > 0) {
        // Split curve segment into sub-segments and test for intersection
        Point cpSplit[7];
        SubdivideBezier(cp, cpSplit);

        // For each of the two segments, see if the ray's bounding box
        // overlaps the segment before recursively checking for
        // intersection with it.
        bool hit = false;
        float u[3] = {u0, (u0 + u1) / 2.f, u1};
        // Pointer to the 4 control poitns for the current segment.
        const Point *cps = cpSplit;
        for (int seg = 0; seg < 2; ++seg, cps += 3) {
            float maxWidth =
                std::max(Lerp(u[seg], common->width[0], common->width[1]),
                         Lerp(u[seg + 1], common->width[0], common->width[1]));

            // As above, check y first, since it most commonly lets us exit
            // out early.
            if (std::max(std::max(cps[0].y, cps[1].y),
                         std::max(cps[2].y, cps[3].y)) +
                        0.5 * maxWidth < 0 ||
                std::min(std::min(cps[0].y, cps[1].y),
                         std::min(cps[2].y, cps[3].y)) -
                        0.5 * maxWidth > 0)
                continue;

            if (std::max(std::max(cps[0].x, cps[1].x),
                         std::max(cps[2].x, cps[3].x)) +
                        0.5 * maxWidth < 0 ||
                std::min(std::min(cps[0].x, cps[1].x),
                         std::min(cps[2].x, cps[3].x)) -
                        0.5 * maxWidth > 0)
                continue;

            float zMax = rayLength * ray.maxt;
            if (std::max(std::max(cps[0].z, cps[1].z),
                         std::max(cps[2].z, cps[3].z)) +
                        0.5 * maxWidth < 0 ||
                std::min(std::min(cps[0].z, cps[1].z),
                         std::min(cps[2].z, cps[3].z)) -
                        0.5 * maxWidth > zMax)
                continue;

            hit |= recursiveIntersect(ray, tHit, rayEpsilon, dg, cps, rayToObject,
                                      u[seg], u[seg + 1], depth - 1);
            // If we found an intersection and this is a shadow ray,
            // we can exit out immediately.
            if (hit && !tHit) return true;
        }
        return hit;
    } else {
        // Intersect ray with curve segment

        // Test ray against segment endpoint boundaries

        // Test sample point against tangent perpendicular at curve start
        float edge =
            (cp[1].y - cp[0].y) * -cp[0].y + cp[0].x * (cp[0].x - cp[1].x);
        if (edge < 0) return false;

        // Test sample point against tangent perpendicular at curve end
        edge = (cp[2].y - cp[3].y) * -cp[3].y + cp[3].x * (cp[3].x - cp[2].x);
        if (edge < 0) return false;

        // Compute line $w$ that gives minimum distance to sample point
        Vector2f segmentDirection = Point2f(cp[3]) - Point2f(cp[0]);
        float denom = segmentDirection.LengthSquared();
        if (denom == 0) return false;
        float w = Dot(-Vector2f(cp[0]), segmentDirection) / denom;

        // Compute $u$ coordinate of curve intersection point and _hitWidth_
        float u = Clamp(Lerp(w, u0, u1), u0, u1);
        float hitWidth = Lerp(u, common->width[0], common->width[1]);
        Normal nHit;
        if (common->type == Ribbon) {
            // Scale _hitWidth_ based on ribbon orientation
            float sin0 = std::sin((1 - u) * common->normalAngle) *
                         common->invSinNormalAngle;
            float sin1 =
                std::sin(u * common->normalAngle) * common->invSinNormalAngle;
            nHit = sin0 * common->n[0] + sin1 * common->n[1];
            hitWidth *= AbsDot(nHit, ray.d) / rayLength;
        }

        // Test intersection point against curve width
        Vector dpcdw;
        Point pc = EvalBezier(cp, Clamp(w, 0.0, 1.0), &dpcdw);
        float ptCurveDist2 = pc.x * pc.x + pc.y * pc.y;
        if (ptCurveDist2 > hitWidth * hitWidth * .25) return false;
        float zMax = rayLength * ray.maxt;
        if (pc.z < 0 || pc.z > zMax) return false;

        // Compute $v$ coordinate of curve intersection point
        float ptCurveDist = std::sqrt(ptCurveDist2);
        float edgeFunc = dpcdw.x * -pc.y + pc.x * dpcdw.y;
        float v = (edgeFunc > 0) ? 0.5f + ptCurveDist / hitWidth
                                 : 0.5f - ptCurveDist / hitWidth;

        // Compute hit _t_ and partial derivatives for curve intersection
        if (tHit != NULL) {
            // FIXME: this tHit isn't quite right for ribbons...
            *tHit = pc.z / rayLength;
            // Compute error bounds for curve intersection
            Vector pError(2 * hitWidth, 2 * hitWidth, 2 * hitWidth);

            // Compute $\dpdu$ and $\dpdv$ for curve intersection
            Vector dpdu, dpdv;
            EvalBezier(common->cpObj, u, &dpdu);
            if (common->type == Ribbon)
                dpdv = Normalize(Cross(nHit, dpdu)) * hitWidth;
            else {
                // Compute curve $\dpdv$ for flat and cylinder curves
                Vector dpduPlane = (Inverse(rayToObject))(dpdu);
                Vector dpdvPlane =
                    Normalize(Vector(-dpduPlane.y, dpduPlane.x, 0)) *
                    hitWidth;
                if (common->type == Cylinder) {
                    // Rotate _dpdvPlane_ to give cylindrical appearance
                    float theta = Lerp(v, -90., 90.);
                    Transform rot = Rotate(-theta, dpduPlane);
                    dpdvPlane = rot(dpdvPlane);
                }
                dpdv = rayToObject(dpdvPlane);
            }
            //*isect = (*ObjectToWorld)(SurfaceInteraction(
            //    ray(pc.z), pError, Point2f(u, v), -ray.d, dpdu, dpdv,
            //    Normal(0, 0, 0), Normal(0, 0, 0), ray.time, this));
            const Transform &o2w = *ObjectToWorld;
            Point pHit = ray(pc.z);
            *dg = DifferentialGeometry(o2w(pHit), o2w(dpdu), o2w(dpdv),
                                       o2w(Normal(0,0,0)), o2w(Normal(0,0,0)),
                                       u, v, this);
#ifdef SHAPE_DEBUG
			Point phit = o2w(pHit);
			vdb_color(1.f, 1.f, 1.f);
			vdb_point(phit[0], phit[1], phit[2]);
			vdb_color(0.f, 1.f, 0.f);
			vdb_line(phit[0], phit[1], phit[2], phit[0] + dg->nn[0], phit[1] + dg->nn[1], phit[2] + dg->nn[2]);
#endif
			*rayEpsilon = 1e-3f * *tHit;
        }
        //++nHits;
        return true;
    }
}

bool Curve::IntersectP(const Ray &r) const {
    float tHit;
    float rayEpsilon;
    DifferentialGeometry dg;
    return Intersect(r, NULL, NULL, NULL);
}

float Curve::Area() const {
    // Compute object-space control points for curve segment, _cpObj_
    Point cpObj[4];
    cpObj[0] = BlossomBezier(common->cpObj, uMin, uMin, uMin);
    cpObj[1] = BlossomBezier(common->cpObj, uMin, uMin, uMax);
    cpObj[2] = BlossomBezier(common->cpObj, uMin, uMax, uMax);
    cpObj[3] = BlossomBezier(common->cpObj, uMax, uMax, uMax);
    float width0 = Lerp(uMin, common->width[0], common->width[1]);
    float width1 = Lerp(uMax, common->width[0], common->width[1]);
    float avgWidth = (width0 + width1) * 0.5f;
    float approxLength = 0.f;
    for (int i = 0; i < 3; ++i)
        approxLength += Distance(cpObj[i], cpObj[i + 1]);
    return approxLength * avgWidth;
}

Point Curve::Sample(float u1, float u2, Normal *Ns) const {
    //LOG(FATAL) << "Curve::Sample not implemented.";
    return Point();
}

std::vector<Reference<Shape> > CreateCurveShape(const Transform *o2w,
                                                     const Transform *w2o,
                                                     bool reverseOrientation,
                                                     const ParamSet &params) {
    float width = params.FindOneFloat("width", 1.f);
    float width0 = params.FindOneFloat("width0", width);
    float width1 = params.FindOneFloat("width1", width);

    int ncp;
    const Point *cp = params.FindPoint("P", &ncp);
    if (ncp != 4) {
        Error(
            "Must provide 4 control points for \"curve\" primitive. "
            "(Provided %d).",
            ncp);
        return std::vector<Reference<Shape> >();
    }

    CurveType type;
    std::string curveType = params.FindOneString("type", "flat");
    if (curveType == "flat")
        type = Flat;
    else if (curveType == "ribbon")
        type = Ribbon;
    else if (curveType == "cylinder")
        type = Cylinder;
    else {
        Error("Unknown curve type \"%s\".  Using \"flat\".", curveType.c_str());
        type = Cylinder;
    }
    int nnorm;
    const Normal *n = params.FindNormal("N", &nnorm);
    if (n != NULL) {
        if (type != Ribbon) {
            Warning("Curve normals are only used with \"ribbon\" type curves.");
            n = NULL;
        } else if (nnorm != 2) {
            Error(
                "Must provide two normals with \"N\" parameter for ribbon "
                "curves. "
                "(Provided %d).",
                nnorm);
            return std::vector<Reference<Shape> >();
        }
    }

    int sd = params.FindOneFloat("splitdepth", 3);

    if (type == Ribbon && !n) {
        Error(
            "Must provide normals \"N\" at curve endpoints with ribbon "
            "curves.");
        return std::vector<Reference<Shape> >();
    } else
        return CreateCurve(o2w, w2o, reverseOrientation, cp, width0, width1,
                           type, n, sd);
}
