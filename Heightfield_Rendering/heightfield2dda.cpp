
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

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


// shapes/heightfield2.cpp*
#include "stdafx.h"
#include "shapes/heightfield2dda.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"

// Heightfield2dda Method Definitions
Heightfield2dda::Heightfield2dda(const Transform *o2w, const Transform *w2o,
        bool ro, int x, int y, const float *zs)
    : Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));

    // Initialize heightfield primitives for grid, maybe its better to use Triangle Class
	int ntris = 2 * (nx - 1) * (ny - 1);
	points = new Point[nx*ny];
	dpdus = new Vector[nx*ny];
	dpdvs = new Vector[nx*ny];
	// Compute heightfield2 vertex positions
	int pos = 0;
	for (int y = 0; y < ny; ++y) {
		for (int x = 0; x < nx; ++x) {
			points[pos].x = (float)x / (float)(nx - 1);
			points[pos].y = (float)y / (float)(ny - 1);
			points[pos].z = z[pos];

            // for differential geometry
            int minx = (x == 0 ? x : x - 1);
            int maxx = (x == nx - 1 ? x : x + 1);
            //dpdus.push_back(Vector(1.f, 0.f, (z[y * nx + maxx] - z[y * nx + minx]) / (maxx - minx) * nx));
			dpdus[pos].x = 1.f;
			dpdus[pos].y = 0.f;
			dpdus[pos].z = (z[y * nx + maxx] - z[y * nx + minx]) / (maxx - minx) * nx;

            int miny = (y == 0 ? y : y - 1);
            int maxy = (y == ny - 1 ? y : y + 1);
            //dpdvs.push_back(Vector(0.f, 1.f, (z[maxy * nx + x] - z[miny * nx + x]) / (maxy - miny) * ny));
			dpdvs[pos].x = 0.f;
			dpdvs[pos].y = 1.f;
			dpdvs[pos].z = (z[maxy * nx + x] - z[miny * nx + x]) / (maxy - miny) * ny;
			++pos;
		}
	}

    // Get the overall bounds
    bounds = Heightfield2dda::ObjectBound();
}


Heightfield2dda::~Heightfield2dda() {
    delete[] z;
	delete[] dpdus;
	delete[] dpdvs;
	delete[] points;
}


BBox Heightfield2dda::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return BBox(Point(0,0,minz), Point(1,1,maxz));
}


//for those can't do intersection test, use RefineFunc()
bool Heightfield2dda::CanIntersect() const {
    return true;
}


//ray & dg is in world space
bool Heightfield2dda::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const {
    //transform ray to object space
    Ray ray;
    (*WorldToObject)(r, &ray);
	*tHit = ray.maxt + 1;

    bool hitSomething = false;
	float rayT;

    // Check ray against overall grid bounds
	if (bounds.Inside(ray(ray.mint)))
		rayT = ray.mint;
	else if (!bounds.IntersectP(ray, &rayT))
	{
		return hitSomething;
	}
	Point gridIntersect = ray(rayT);

    // Set up 2D DDA for ray
	float NextCrossingT[2], DeltaT[2];
	int Step[2], Out[2], Pos[2];
	float width[2] = { 1.f / (nx - 1), 1.f / (ny - 1) };
	int invWidth[2] = { nx - 1, ny - 1 };

	for (int axis = 0; axis < 2; ++axis)
	{
		// Compute current voxel for axis
		Pos[axis] = Float2Int((gridIntersect[axis] - bounds.pMin[axis]) * invWidth[axis]);
		Pos[axis] = Clamp(Pos[axis], 0, invWidth[axis] - 1);

		if (ray.d[axis] >= 0)
		{
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT + 
				((Pos[axis] + 1) * width[axis] - (gridIntersect[axis] - bounds.pMin[axis])) / ray.d[axis];
			DeltaT[axis] = width[axis] / ray.d[axis];
			Step[axis] = 1;
			Out[axis] = invWidth[axis];
		}
		else
		{
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT + 
				((Pos[axis] + 1) * width[axis] - (gridIntersect[axis] - bounds.pMin[axis])) / ray.d[axis];
			DeltaT[axis] = -width[axis] / ray.d[axis];
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}

    // Walk ray through voxel grid
#define VERT(x,y) ((x)+(y)*nx)
	while (1)
	{
		// Check for intersection in current grid and advance to next
		// in every grid we have exactly 2 triangles
		Point p1, p2, p3;
		float t = ray.maxt + 1;
		int x = Pos[0], y = Pos[1];

		//the first triangle
		p1 = points[VERT(x, y)];
		p2 = points[VERT(x + 1, y)];
		p3 = points[VERT(x + 1, y + 1)];
		hitSomething |= TriangleIntersect(ray, p1, p2, p3, t);
		if (*tHit > t)
		{
			*tHit = t;
			// get the paramset
			Point phit = ray(t);
			*dg = DifferentialGeometry((*ObjectToWorld)(phit),
				(*ObjectToWorld)((p2 - p1) * invWidth[0]), (*ObjectToWorld)((p3 - p2) * invWidth[1]),
				Normal(0, 0, 0), Normal(0, 0, 0),
				phit.x, phit.y, this);
			*rayEpsilon = 1e-3f * *tHit;
		}
		
		//the second triangle
		p1 = points[VERT(x, y)];
		p2 = points[VERT(x + 1, y + 1)];
		p3 = points[VERT(x, y + 1)];
		hitSomething |= TriangleIntersect(ray, p1, p2, p3, t);
		if (*tHit > t) 
		{
			*tHit = t;
			// get the paramset
			Point phit = ray(t);
			*dg = DifferentialGeometry((*ObjectToWorld)(phit),
				(*ObjectToWorld)(p2 - p1) * invWidth[0], (*ObjectToWorld)(p3 - p2) * invWidth[1],
				Normal(0, 0, 0), Normal(0, 0, 0),
				phit.x, phit.y, this);
			*rayEpsilon = 1e-3f * *tHit;
		}

		// Specified Optimization for Heightfield
		// because heightfield is fully organized
		//if (hitSomething)
			//break;

		// Advance to next voxel
		// Find _stepAxis_ for stepping to next voxel
		int bits = (NextCrossingT[0] < NextCrossingT[1]) ? 0 : 1;
		if (ray.maxt < NextCrossingT[bits])
			break;
		Pos[bits] += Step[bits];
		if (Pos[bits] == Out[bits])
			break;
		NextCrossingT[bits] += DeltaT[bits];
	}
#undef VERT
    return hitSomething;
}


bool Heightfield2dda::IntersectP(const Ray &r) const {
	float tHit;
	float rayEpsilon;
	DifferentialGeometry dg;
	return Intersect(r, &tHit, &rayEpsilon, &dg);
}

bool Heightfield2dda::TriangleIntersect(const Ray &r, const Point &p1, const Point &p2, const Point &p3, float &tHit) const{
	Ray ray = r;

	Vector e1 = p2 - p1;
	Vector e2 = p3 - p1;
	Vector s1 = Cross(ray.d, e2);
	float divisor = Dot(s1, e1);
	
	if (divisor == 0.)
		return false;
	float invDivisor = 1.f / divisor;
	
	// Compute first barycentric coordinate
	Vector s = ray.o - p1;
	float b1 = Dot(s, s1) * invDivisor;
	if (b1 < 0. || b1>1.)
		return false;
		
	// Compute second barycentric coordinate
	Vector s2 = Cross(s, e1);
	float b2 = Dot(ray.d, s2) * invDivisor;
	if (b2 < 0. || b1 + b2 > 1.)
		return false;
		
	// Compute _t_ to intersection point
	float t = Dot(e2, s2) * invDivisor;
	if (t < ray.mint || t > ray.maxt)
		return false;

	tHit = t;
}

bool Heightfield2dda::TriangleIntersectP(const Ray &r, const Point &p1, const Point &p2, const Point &p3) const {
	return false;
}

//this compute the triangle mesh
void Heightfield2dda::Refine(vector<Reference<Shape> > &refined) const {
    int ntris = 2*(nx-1)*(ny-1);
    refined.reserve(ntris);
    int *verts = new int[3*ntris];
    Point *P = new Point[nx*ny];
    float *uvs = new float[2*nx*ny];
    int nverts = nx*ny;
    int x, y;
    // Compute heightfield2 vertex positions
    int pos = 0;
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
            P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
            P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
            P[pos].z = z[pos];
            ++pos;
        }
    }

    // Fill in heightfield2 vertex offset array
    int *vp = verts;
    for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {
#define VERT(x,y) ((x)+(y)*nx)
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y);
            *vp++ = VERT(x+1, y+1);
    
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y+1);
            *vp++ = VERT(x, y+1);
        }
#undef VERT
    }
    ParamSet paramSet;
    paramSet.AddInt("indices", verts, 3*ntris);
    paramSet.AddFloat("uv", uvs, 2 * nverts);
    paramSet.AddPoint("P", P, nverts);
    refined.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));
    delete[] P;
    delete[] uvs;
    delete[] verts;
}

//void Heightfield2dda::GetShadingGeometry(const Transform &o2w, const DifferentialGeometry &dg, DifferentialGeometry *dgShading) const {
//    float u = dg.u * (nx - 1), v = dg.v * (ny - 1);
//    int i = floor(u), j = floor(v);
//    float uu = u - i, vv = v - j;
//
//	Vector dpdu, dpdv;
//#define VERT(x,y) ((x)+(y)*nx)
//	// this interpolation is better
//    if (uu > vv){
//        dpdu = (1 - uu) * dpdus[VERT(i, j)] + (uu - vv) * dpdus[VERT(i + 1, j)] + vv * dpdus[VERT(i+1, j+1)];
//        dpdv = (1 - uu) * dpdvs[VERT(i, j)] + (uu - vv) * dpdvs[VERT(i + 1, j)] + vv * dpdvs[VERT(i+1, j+1)];
//    }
//    else
//    {
//        dpdu = (1 - vv) * dpdus[VERT(i, j)] + (vv - uu) * dpdus[VERT(i, j + 1)] + uu * dpdus[VERT(i+1, j+1)];
//        dpdv = (1 - vv) * dpdvs[VERT(i, j)] + (vv - uu) * dpdvs[VERT(i, j + 1)] + uu * dpdvs[VERT(i+1, j+1)];
//    }
//#undef VERT
//    *dgShading = DifferentialGeometry((dg.p),
//                        (o2w)(dpdu), (o2w)(dpdv),
//                        Normal(0, 0, 0), Normal(0, 0, 0),
//                        dg.u, dg.v, this);
//}

Heightfield2dda *CreateHeightfield2ddaShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield2dda(o2w, w2o, reverseOrientation, nu, nv, Pz);
}


