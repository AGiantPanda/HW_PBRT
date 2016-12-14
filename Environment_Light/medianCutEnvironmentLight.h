
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_LIGHTS_medianCutEnvironmentLight_H
#define PBRT_LIGHTS_medianCutEnvironmentLight_H

// lights/medianCutEnvironmentLight.h*
#include <vector>
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"

struct RegionLight {
	Point regionPos;		// the upper-left of the region
	int width;				// width of the region
	int height;				// height of the region
	float energy;			// total energy of the region
	Point lightPos;			// the light source position
	RGBSpectrum lightColor; // the light source color
	RegionLight(Point rP, int w, int h, float e = 0.f, Point lP = Point(0.f, 0.f, 0.f))
	{
		regionPos = rP;
		width = w;
		height = h;
		energy = e;
		lightPos = lP;
	}
};

// medianCutEnvironmentLight Declarations
class medianCutEnvironmentLight : public Light {
public:
    // medianCutEnvironmentLight Public Methods
	medianCutEnvironmentLight(const Transform &light2world, const Spectrum &power, int ns,
        const string &texmap);
    ~medianCutEnvironmentLight();
    Spectrum Power(const Scene *) const;
    bool IsDeltaLight() const { return true; }
    Spectrum Le(const RayDifferential &r) const;
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls,
        float time, Vector *wi, float *pdf, VisibilityTester *visibility) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
        float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
    void SHProject(const Point &p, float pEpsilon, int lmax, const Scene *scene,
        bool computeLightVis, float time, RNG &rng, Spectrum *coeffs) const;
private:
    // medianCutEnvironmentLight Private Data
    MIPMap<RGBSpectrum> *radianceMap;
    Distribution2D *distribution;
	int Width;
	int Height;
	float *energyMatrix;
	RGBSpectrum *colorMatrix;
	std::vector<RegionLight> regionLight;
	float constantPdf;
	float regionNums;

	// medianCutEnvironmentLight Private Methods
	float GetEnergy(int ul, int vl, int ur, int vr);
	RGBSpectrum GetColor(int ul, int vl, int ur, int vr);
	RegionLight DivideRegion(RegionLight &r1);
	void GetCentralPoint(RegionLight &rl);
	void GetCentroidPoint(RegionLight &rl);
};


medianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet);

#endif // PBRT_LIGHTS_medianCutEnvironmentLight_H
