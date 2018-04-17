
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

/*

See the writeup "The Implementation of a Hair Scattering Model" at
http://pbrt.org/hair.pdf for a description of the implementation here.

*/

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_HAIR_H
#define PBRT_MATERIALS_HAIR_H

// materials/hair.h*
#include "material.h"
#include "pbrt.h"
#include "reflection.h"
#include <array>


// HairMaterial Declarations
class HairMaterial : public Material {
  public:
    // HairMaterial Public Methods
    HairMaterial(Reference<Texture<Spectrum>> sa,
                 Reference<Texture<Spectrum>> c,
                 Reference<Texture<float>> eu,
                 Reference<Texture<float>> pheom,
                 Reference<Texture<float>> et,
                 Reference<Texture<float>> bm,
                 Reference<Texture<float>> bn,
                 Reference<Texture<float>> a){
          sigma_a = sa;
          color = c;
          eumelanin = eu;
          pheomelanin = pheom,
          eta = et;
          beta_m = bm;
          beta_n = bn;
          alpha = a; 
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const;


  private:
    // HairMaterial Private Data
    Reference<Texture<Spectrum>> sigma_a, color;
    Reference<Texture<float>> eumelanin, pheomelanin, eta;
    Reference<Texture<float>> beta_m, beta_n, alpha;
};


HairMaterial *CreateHairMaterial(const Transform &xform,const TextureParams &mp);

// HairBSDF Constants
static const int pMax = 3;
static const float SqrtPiOver8 = 0.626657069f;

// HairBSDF Declarations
class HairBSDF : public BxDF {
  public:
    // HairBSDF Public Methods
    HairBSDF(float h, float eta, const Spectrum &sigma_a, float beta_m,
             float beta_n, float alpha);
    Spectrum f(const Vector &wo, const Vector &wi) const;
    Spectrum Sample_f(const Vector &wo, Vector *wi, const Point2f &u,
                      float *pdf, BxDFType *sampledType) const;
    float Pdf(const Vector &wo, const Vector &wi) const;
    std::string ToString() const;
    static Spectrum SigmaAFromConcentration(float ce, float cp);
    static Spectrum SigmaAFromReflectance(const Spectrum &c, float beta_n);

  private:
    // HairBSDF Private Methods
    std::array<float, pMax + 1> ComputeApPdf(float cosThetaO) const;

    // HairBSDF Private Data
    const float h, gammaO, eta;
    const Spectrum sigma_a;
    const float beta_m, beta_n, alpha;
    float v[pMax + 1];
    float s;
    float sin2kAlpha[3], cos2kAlpha[3];
};

// General Utility Functions
inline float Sqr(float v) { return v * v; }
template <int n>
static float Pow(float v) {
    static_assert(n > 0, "Power can't be negative");
    float n2 = Pow<n / 2>(v);
    return n2 * n2 * Pow<n & 1>(v);
}

template <>
inline float Pow<1>(float v) {
    return v;
}
template <>
inline float Pow<0>(float v) {
    return 1;
}

inline float SafeASin(float x) {
    assert(x >= -1.0001 && x <= 1.0001);
    return std::asin(Clamp(x, -1.f, 1.f));
}

inline float SafeSqrt(float x) {
    //CHECK_GE(x, -1e-4);
    return std::sqrt(std::max(float(0), x));
}

// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
static uint32_t Compact1By1(uint32_t x) {
    // TODO: as of Haswell, the PEXT instruction could do all this in a
    // single instruction.
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x &= 0x55555555;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >> 1)) & 0x33333333;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >> 2)) & 0x0f0f0f0f;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >> 4)) & 0x00ff00ff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x >> 8)) & 0x0000ffff;
    return x;
}
//-----------------------------------------------------------
// Decomposes a sample ξ ∈ [0, 1) into a pair of samples
//-----------------------------------------------------------
static Point2f Demuxfloat(float f) {
    assert(f >= 0 && f < 1);
    uint64_t v = f * (1ull << 32);
    //CHECK_LT(v, 0x100000000);
    uint32_t bits[2] = {Compact1By1(v), Compact1By1(v >> 1)};
    return {bits[0] / float(1 << 16), bits[1] / float(1 << 16)};
}


#endif  // PBRT_MATERIALS_HAIR_H
