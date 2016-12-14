
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


// lights/medianCutEnvironmentLight.cpp*
#include "stdafx.h"
#include "lights/medianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

// This is for debug usage, PanDeBug stands for "Panda's Debug", not "Pan Debug".
// Yes, Im boring
//#define PanDeBug

// medianCutEnvironmentLight Utility Classes
struct medianCutEnvironmentCube {
    // InfiniteAreaCube Public Methods
	medianCutEnvironmentCube(const medianCutEnvironmentLight *l, const Scene *s,
                     float t, bool cv, float pe)
        : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const medianCutEnvironmentLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};



// medianCutEnvironmentLight Method Definitions
medianCutEnvironmentLight::~medianCutEnvironmentLight() {
	delete[] energyMatrix;
	delete[] colorMatrix;
	delete distribution;
    delete radianceMap;
}


medianCutEnvironmentLight::medianCutEnvironmentLight(const Transform &light2world,
        const Spectrum &L, int ns, const string &texmap)
    : Light(light2world, ns) {
	constantPdf = 0.0f;
	regionNums = ns;

    int width = 0, height = 0;
    RGBSpectrum *texels = NULL;
    // Read texel data from _texmap_ into _texels_, the original probe image is stored here
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
	Width = width;
	Height = height;
	float normalize = (2 * M_PI) * (M_PI) / float(Height * Width);

// This shows the original texel
#ifdef PanDeBug
	float *debugimg = new float[width * height * 3];
	for (int v = 0; v < height; v++) {
		for (int u = 0; u < width; u++)
			texels[u + v*width].ToRGB(debugimg + (u + v*width) * 3);
	}
	WriteImage("texmaptest.png", debugimg, NULL, width,
		height, 0, 0, 0, 0);
	printf("width %d height %d\n", width, height);
	fflush(stdin);
#endif

	// Compute normalized texel
	for (int v = 0; v < height; v++) {
		float sinTheta = sinf(M_PI * float(v + .5f) / float(height));
		for (int u = 0; u < width; u++) {
			texels[u + v*width] = texels[u + v*width] * sinTheta * normalize;
		}
	}

    // Initialize sampling PDFs for infinite area light

    // Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    float *img = new float[width*height];
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            img[u+v*width] = radianceMap->Lookup(up, vp, filter).y();
            img[u+v*width] *= sinTheta * normalize;
        }
    }

    // Compute sampling distributions for rows and columns of image
    distribution = new Distribution2D(img, width, height);

#define VERT(u, v) ((u) + (v) * Width)
	// Compute energyMatrix and colorMatrix from environment map
	energyMatrix = new float[width * height];
	colorMatrix = new RGBSpectrum[width * height];
	energyMatrix[0] = img[0];
	colorMatrix[0] = texels[0];
	for (int v = 1; v < height; v++) {
		energyMatrix[VERT(0, v)] = energyMatrix[VERT(0, v - 1)] + img[VERT(0, v)];
		colorMatrix[VERT(0, v)] = colorMatrix[VERT(0, v - 1)] + texels[VERT(0, v)];
	}
	for (int u = 1; u < width; u++) {
		energyMatrix[VERT(u, 0)] = energyMatrix[VERT(u - 1, 0)] + img[VERT(u, 0)];
		colorMatrix[VERT(u, 0)] = colorMatrix[VERT(u - 1, 0)] + texels[VERT(u, 0)];
	}
	for (int v = 1; v < height; v++)
	{
		for (int u = 1; u < width; u++)
		{
			energyMatrix[VERT(u, v)] = energyMatrix[VERT(u, v - 1)] + energyMatrix[VERT(u - 1, v)]
					   				 - energyMatrix[VERT(u - 1, v - 1)]
									 + img[VERT(u, v)];
			colorMatrix[VERT(u, v)] = colorMatrix[VERT(u, v - 1)] + colorMatrix[VERT(u - 1, v)]
									- colorMatrix[VERT(u - 1, v - 1)]
									+ texels[VERT(u, v)];
		}
	}
#undef VERT

	// Subdivide the region until nSamples are enough
	regionLight.push_back(RegionLight(Point(0.f, 0.f, 0.f), width - 1, height - 1));
	for (int n = regionLight.size(); n < regionNums;)
	{
		for (int i = 0; i < n; i++)
		{
			regionLight.push_back(DivideRegion(regionLight[i]));
		}
		n = regionLight.size();
	}
	
// Shows the divide result
#ifdef PanDeBug
	for (int i = 0; i < regionLight.size(); i++) {
		int ul = regionLight[i].regionPos.x, vl = regionLight[i].regionPos.y;
		int ur = ul + regionLight[i].width, vr = vl + regionLight[i].height;
		for (int u = ul; u <= ur; u++) {
			debugimg[(u + vl*width) * 3] = 0.0f;
			debugimg[(u + vl*width) * 3+1] = 1.0f;
			debugimg[(u + vl*width) * 3+2] = 0.0f;

			debugimg[(u + vr*width) * 3] = 0.0f;
			debugimg[(u + vr*width) * 3+1] = 1.0f;
			debugimg[(u + vr*width) * 3+2] = 0.0f;
		}
		for (int v = vl; v <= vr; v++) {
			debugimg[(ul + v*width) * 3] = 0.0f;
			debugimg[(ul + v*width) * 3 + 1] = 1.0f;
			debugimg[(ul + v*width) * 3 + 2] = 0.0f;

			debugimg[(ur + v*width) * 3] = 0.0f;
			debugimg[(ur + v*width) * 3 + 1] = 1.0f;
			debugimg[(ur + v*width) * 3 + 2] = 0.0f;
		}
		printf("i: %d, energy = %f\n", i, regionLight[i].energy);
	}
	WriteImage("medianCuttest.png", debugimg, NULL, width,
		height, 0, 0, 0, 0);
	printf("width %d height %d\n", width, height);
	fflush(stdin);
#endif

	// Get Center/Centroid Point Attributes of each region
	// To compute the centroid point could time comsuming so ill put it here
	for (int i = 0; i < regionLight.size(); i++)
	{
		//GetCentralPoint(regionLight[i]);
		//GetCentroidPoint();
		RGBSpectrum color = RGBSpectrum(0.f);
		float centroidV = 0.f, centroidU = 0.f, totalE = 0.f;
#define VERT(u, v) ((u)+(v) * width)
		for (int v = regionLight[i].regionPos.y; v <= regionLight[i].regionPos.y + regionLight[i].height; v++) {
			for (int u = regionLight[i].regionPos.x; u <= regionLight[i].regionPos.x + regionLight[i].width; u++) {
				color += texels[VERT(u, v)];
				float f = img[VERT(u, v)] * img[VERT(u, v)];
				centroidV += v * f, centroidU += u * f, totalE += f;
			}
		}
#undef VERT
		regionLight[i].lightPos.x = centroidU / totalE / float(Width);
		regionLight[i].lightPos.y = centroidV / totalE / float(Height);
		regionLight[i].lightColor = color;
	}

// Shows the light source
#ifdef PanDeBug
	for (int i = 0; i < regionLight.size(); i++) {
		int u = regionLight[i].lightPos.x, v = regionLight[i].lightPos.y;
		debugimg[(u + v*width) * 3] = 1.0f;
		debugimg[(u + v*width) * 3 + 1] = 1.0f;
		debugimg[(u + v*width) * 3 + 2] = 1.0f;
	}
	WriteImage("lightSourcetest.png", debugimg, NULL, width,
		height, 0, 0, 0, 0);
	printf("width %d height %d\n", width, height);
	fflush(stdin);
#endif

	constantPdf = 1.0 / float(regionLight.size());
	//
	delete[] img;
#ifdef PanDeBug
	delete[] debugimg;
#endif
	delete[] texels;
}


Spectrum medianCutEnvironmentLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
        Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum medianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void medianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    // Project _InfiniteAreaLight_ to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _InfiniteAreaLight_ to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _InfiniteAreaLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                        (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project _InfiniteAreaLight_ to SH from cube map sampling
        SHProjectCube(medianCutEnvironmentCube(this, scene, time, computeLightVis,
                                       pEpsilon),
                      p, 200, lmax, coeffs);
    }
}


medianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new medianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum medianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // randomly choose one light source from RS
	RegionLight rl = regionLight[Floor2Int(ls.uComponent * regionLight.size())];

    // Convert infinite light sample point to direction
    float theta = rl.lightPos.y * M_PI, phi = rl.lightPos.x * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                              costheta));

    // Compute PDF for sampled infinite light direction
	*pdf = constantPdf;

    // Return radiance value for infinite light direction
    visibility->SetRay(p, pEpsilon, *wi, time);
    Spectrum Ls = Spectrum(rl.lightColor,
                           SPECTRUM_ILLUMINANT);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


float medianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
    PBRT_INFINITE_LIGHT_STARTED_PDF();
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
           (2.f * M_PI * M_PI * sintheta);
	PBRT_INFINITE_LIGHT_FINISHED_PDF();
    return p;
}



float medianCutEnvironmentLight::GetEnergy(int ul, int vl, int ur, int vr) {
#define VERT(u, v) ((u) + (v) * Width)
	float up, left, corner;
	up = (vl == 0) ? 0 : energyMatrix[VERT(ur, vl - 1)];
	left = (ul == 0) ? 0 : energyMatrix[VERT(ul - 1, vr)];
	corner = (ul == 0 || vl == 0) ? 0 : energyMatrix[VERT(ul - 1, vl - 1)];
	return energyMatrix[VERT(ur, vr)] - up - left + corner;
#undef VERT
}

RGBSpectrum medianCutEnvironmentLight::GetColor(int ul, int vl, int ur, int vr) {
#define VERT(u, v) ((u) + (v) * Width)
	RGBSpectrum up, left, corner;
	up = (vl == 0) ? 0 : colorMatrix[VERT(ur, vl - 1)];
	left = (ul == 0) ? 0 : colorMatrix[VERT(ul - 1, vr)];
	corner = (ul == 0 || vl == 0) ? 0 : colorMatrix[VERT(ul - 1, vl - 1)];
	return (colorMatrix[VERT(ur, vr)] - up - left + corner);
#undef VERT
}

RegionLight medianCutEnvironmentLight::DivideRegion(RegionLight &r1) {
	RegionLight r2(Point(0.f, 0.f, 0.f), 0, 0);
	int ul = r1.regionPos.x, vl = r1.regionPos.y;
	int ur = r1.regionPos.x + r1.width, vr = r1.regionPos.y + r1.height;
	float midWidth = vl + r1.height / 2;
	midWidth = r1.width * sinf(M_PI * float(midWidth + .5f) / float(Height));

	float es = GetEnergy(ul, vl, ur, vr) / 2;
	if (midWidth > r1.height) {
		// divide along the width
		int l = ul;
		int mid = ul + r1.width / 2;
		int r = ur;
		float leftEnergy = GetEnergy(ul, vl, mid, vr);
		if (ul == ur)
		{
			r1.energy = leftEnergy / 2.0f;
			r2.regionPos.x = ul;
			r2.regionPos.y = vl;
			r2.width = r1.width;
			r2.height = r1.height;
			r2.energy = r1.energy;
		}
		else
		{
			while (r - l > 1 && leftEnergy != es) {
				if (leftEnergy > es) {
					r = mid;
				}
				else {
					l = mid;
				}
				mid = (l + r) / 2;
				leftEnergy = GetEnergy(ul, vl, mid, vr);
			}
			r1.width = mid - ul;
			r1.energy = leftEnergy;
			r2.regionPos.x = mid + 1;
			r2.regionPos.y = vl;
			r2.width = ur - (mid + 1);
			r2.height = r1.height;
			r2.energy = GetEnergy(mid + 1, vl, ur, vr);
		}
	}
	else {
		// divide along the height
		int u = vl;
		int mid = vl + r1.height / 2;
		int b = vr;
		float upperEnergy = GetEnergy(ul, vl, ur, mid);
		if (vl == vr)
		{
			r1.energy = upperEnergy / 2.0f;
			r2.regionPos.x = ul;
			r2.regionPos.y = vl;
			r2.width = r1.width;
			r2.height = r1.height;
			r2.energy = r1.energy;
		}
		else
		{
			while (b - u > 1 && upperEnergy != es) {
				if (upperEnergy > es) {
					b = mid;
				}
				else {
					u = mid;
				}
				mid = (u + b) / 2;
				upperEnergy = GetEnergy(ul, vl, ur, mid);
			}
			r1.height = mid - vl;
			r1.energy = upperEnergy;
			r2.regionPos.x = ul;
			r2.regionPos.y = mid + 1;
			r2.width = r1.width;
			r2.height = vr - (mid + 1);
			r2.energy = GetEnergy(ul, mid + 1, ur, vr);
		}
	}
	return r2;
}

void medianCutEnvironmentLight::GetCentralPoint(RegionLight &rl) {
	rl.lightPos.x = float(rl.regionPos.x + rl.width / 2) / float(Width);
	rl.lightPos.y = float(rl.regionPos.y + rl.height / 2) / float(Height);
	rl.lightColor = GetColor(rl.regionPos.x, rl.regionPos.y, rl.regionPos.x + rl.width, rl.regionPos.y + rl.height);
}

void medianCutEnvironmentLight::GetCentroidPoint(RegionLight &rl) {

}

// Looks not used?? leave it alone first
Spectrum medianCutEnvironmentLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Compute direction for infinite light sample ray

    // Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) return Spectrum(0.f);

    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                    costheta));
    *Ns = (Normal)d;

    // Compute origin for infinite light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

    // Compute _InfiniteAreaLight_ ray PDF
    float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;
    if (sintheta == 0.f) *pdf = 0.f;
    Spectrum Ls = (radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}

#ifdef PanDeBug
#undef PanDeBug
#endif