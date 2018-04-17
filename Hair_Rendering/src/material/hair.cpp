
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

// materials/hair.cpp*
#include "stdafx.h"
#include <array>
#include <numeric>
//#include "interaction.h"
#include "materials/hair.h"
#include "paramset.h"
#include "reflection.h"
//#include "sampling.h"
#include "spectrum.h"
#include "texture.h"
#include "textures/constant.h"
//[queenie]
#include <stdio.h>

// Hair Local Declarations
inline float I0(float x), LogI0(float x);

// Hair Local Functions
//-----------------------------------
// Longitudinal Scattering (Theta)
//-----------------------------------
static float Mp(float cosThetaI, float cosThetaO, float sinThetaI,
                float sinThetaO, float v) {
    float a = cosThetaI * cosThetaO / v;
    float b = sinThetaI * sinThetaO / v;
    float mp =
        (v <= .1)
            ? (std::exp(LogI0(a) - b - 1 / v + 0.6931f + std::log(1 / (2 * v))))
            : (std::exp(-b) * I0(a)) / (std::sinh(1 / v) * 2 * v);
    assert(!isinf(mp) && !isnan(mp));
    return mp;
}

inline float I0(float x) {
    float val = 0;
    float x2i = 1;
    int ifact = 1;
    int i4 = 1;
    // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
    for (int i = 0; i < 10; ++i) {
        if (i > 1) ifact *= i;
        val += x2i / (i4 * Sqr(ifact));
        x2i *= x * x;
        i4 *= 4;
    }
    return val;
}

inline float LogI0(float x) {
    if (x > 12)
        return x + 0.5 * (-std::log(2 * M_PI) + std::log(1 / x) + 1 / (8 * x));
    else
        return std::log(I0(x));
}

//-----------------------------------
// Attenuation(h = sinGamma)
//-----------------------------------
static std::array<Spectrum, pMax + 1> Ap(float cosThetaO, float eta, float h,
                                         const Spectrum &T) {
    std::array<Spectrum, pMax + 1> ap;
    // Compute $p=0$ attenuation at initial cylinder intersection
    float cosGammaO = SafeSqrt(1 - h * h);
    float cosTheta = cosThetaO * cosGammaO;
    // ------------------------------------------------------------
    // compute Fresnel reflection formula for dielectric materials 
    // input: cosThetaI, etaI, etaT
    // ------------------------------------------------------------
    float f = FrDielectric(cosTheta, 1.f, eta);
    
    ap[0] = f;

    // Compute $p=1$ attenuation term
    ap[1] = Sqr(1 - f) * T;

    // Compute attenuation terms up to $p=_pMax_$
    for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;

    // Compute attenuation term accounting for remaining orders of scattering
    ap[pMax] = ap[pMax - 1] * f * T / (Spectrum(1.f) - T * f);
    return ap;
}

//----------------------------------------------------------
// the angle between two ray(original and reflection/refraction)
//----------------------------------------------------------
inline float Phi(int p, float gammaO, float gammaT) {
    return 2 * p * gammaT - 2 * gammaO + p * M_PI;
}

//---------------------------------------------------------------------
// Use logistic distribution to represent the rouughness of surface
//---------------------------------------------------------------------
inline float Logistic(float x, float s) {
    x = std::abs(x);
    return std::exp(-x / s) / (s * Sqr(1 + std::exp(-x / s)));
}
//---------------------------------------------------------------------
// Integration of logistic distribution(CDF)
//---------------------------------------------------------------------
inline float LogisticCDF(float x, float s) {
    return 1 / (1 + std::exp(-x / s));
}
//---------------------------------------------------------------------
// Normalized logistic function over a range [a, b]
//---------------------------------------------------------------------
inline float TrimmedLogistic(float x, float s, float a, float b) {
    //CHECK_LT(a, b);
    return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}
//---------------------------------------------------------------------
// Computing the angular difference between φ and Φ(p, h) 
// and evaluating the azimuthal distribution with that angle
//---------------------------------------------------------------------
inline float Np(float phi, int p, float s, float gammaO, float gammaT) {
    float dphi = phi - Phi(p, gammaO, gammaT);
    // Remap _dphi_ to $[-\pi,\pi]$
    while (dphi > M_PI) dphi -= 2 * M_PI;
    while (dphi < -M_PI) dphi += 2 * M_PI;
    return TrimmedLogistic(dphi, s, -M_PI, M_PI);
}

static float SampleTrimmedLogistic(float u, float s, float a, float b) {
    //CHECK_LT(a, b);
    float k = LogisticCDF(b, s) - LogisticCDF(a, s);
    float x = -s * std::log(1 / (u * k + LogisticCDF(a, s)) - 1);
    assert(!isnan(x));
    return Clamp(x, a, b);
}

// HairMaterial Method Definitions
// Change it into GetBSDF
//---------------------------------------------------------------------------
// h: sin_gamma(入射角) = h 
// eta: the refraction index of the interior of hair
// sigma_a : the absorption coefficient of the hair interior
// beta_m : longitudinal roughness of the hair
// beta_n : azimuthal roughness
// alpha: the angle that the small scales on the surface of hair are offset from the base cylinder(degrees)
//---------------------------------------------------------------------------
BSDF *HairMaterial::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    DifferentialGeometry dgs;
    dgs = dgShading;
    float e = eta->Evaluate(dgs);
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn, e);

    float bm = beta_m->Evaluate(dgs);
    float bn = beta_n->Evaluate(dgs);
    float a = Radians(alpha->Evaluate(dgs));

    Spectrum sig_a;
    if (sigma_a)
        sig_a = sigma_a->Evaluate(dgs).Clamp();
    else if (color) {
        Spectrum c = color->Evaluate(dgs).Clamp();
        sig_a = HairBSDF::SigmaAFromReflectance(c, bn);
    } else {
        assert(eumelanin || pheomelanin);
        sig_a = HairBSDF::SigmaAFromConcentration(
            std::max(float(0), eumelanin ? eumelanin->Evaluate(dgs) : 0),
            std::max(float(0), pheomelanin ? pheomelanin->Evaluate(dgs) : 0));
    }

    // Offset along width
    float h = -1 + 2 * dgs.v;
    bsdf->Add(BSDF_ALLOC(arena, HairBSDF)(h, e, sig_a, bm, bn, a));
    return bsdf;
}

HairMaterial *CreateHairMaterial(const Transform &xform, const TextureParams &mp) {
    Reference<Texture<Spectrum>> sigma_a =
        mp.GetSpectrumTextureOrNull("sigma_a");
    Reference<Texture<Spectrum>> color =
        mp.GetSpectrumTextureOrNull("color");
    Reference<Texture<float>> eumelanin =
        mp.GetFloatTextureOrNull("eumelanin");
    Reference<Texture<float>> pheomelanin =
        mp.GetFloatTextureOrNull("pheomelanin");
    if (sigma_a) {
        if (color)
            Warning(
                "Ignoring \"color\" parameter since \"sigma_a\" was provided.");
        if (eumelanin)
            Warning(
                "Ignoring \"eumelanin\" parameter since \"sigma_a\" was "
                "provided.");
        if (pheomelanin)
            Warning(
                "Ignoring \"pheomelanin\" parameter since \"sigma_a\" was "
                "provided.");
    } else if (color) {
        if (sigma_a)
            Warning(
                "Ignoring \"sigma_a\" parameter since \"color\" was provided.");
        if (eumelanin)
            Warning(
                "Ignoring \"eumelanin\" parameter since \"color\" was "
                "provided.");
        if (pheomelanin)
            Warning(
                "Ignoring \"pheomelanin\" parameter since \"color\" was "
                "provided.");
    } else if (eumelanin || pheomelanin) {
        if (sigma_a)
            Warning(
                "Ignoring \"sigma_a\" parameter since "
                "\"eumelanin\"/\"pheomelanin\" was provided.");
        if (color)
            Warning(
                "Ignoring \"color\" parameter since "
                "\"eumelanin\"/\"pheomelanin\" was provided.");
    } else {
        // Default: brown-ish hair.
        sigma_a = new ConstantTexture<Spectrum>(HairBSDF::SigmaAFromConcentration(1.3, 0.));
    }

    Reference<Texture<float>> eta = mp.GetFloatTexture("eta", 1.55f);
    Reference<Texture<float>> beta_m = mp.GetFloatTexture("beta_m", 0.3f);
    Reference<Texture<float>> beta_n = mp.GetFloatTexture("beta_n", 0.3f);
    Reference<Texture<float>> alpha = mp.GetFloatTexture("alpha", 2.f);

    return new HairMaterial(sigma_a, color, eumelanin, pheomelanin, eta, beta_m,
                            beta_n, alpha);
}

// HairBSDF Method Definitions
HairBSDF::HairBSDF(float h, float eta, const Spectrum &sigma_a, float beta_m,
                   float beta_n, float alpha)
    : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)),
      h(h),
      gammaO(SafeASin(h)),
      eta(eta),
      sigma_a(sigma_a),
      beta_m(beta_m),
      beta_n(beta_n),
      alpha(alpha) {
    assert(h >= -1 && h <= 1);
    assert(beta_m >= 0 && beta_m <= 1);
    assert(beta_n >= 0 && beta_n <= 1);
    //------------------------------------------------
    // Compute longitudinal variance from $\beta_m$]
    //------------------------------------------------
    static_assert(
        pMax >= 3,
        "Longitudinal variance code must be updated to handle low pMax");
    v[0] = Sqr(0.726f * beta_m + 0.812f * Sqr(beta_m) + 3.7f * Pow<20>(beta_m));
    v[1] = .25 * v[0];
    v[2] = 4 * v[0];
    for (int p = 3; p <= pMax; ++p)
        // TODO: is there anything better here?
        v[p] = v[2];
    //-------------------------------------------------------
    // Compute azimuthal logistic scale factor from $\beta_n$
    //-------------------------------------------------------
    s = SqrtPiOver8 *
        (0.265f * beta_n + 1.194f * Sqr(beta_n) + 5.372f * Pow<22>(beta_n));
    assert(!isnan(s));

    // Compute $\alpha$ terms for hair scales
    sin2kAlpha[0] = std::sin(alpha);
    cos2kAlpha[0] = SafeSqrt(1 - Sqr(sin2kAlpha[0]));
    //compute sin(2^k α) and cos(2^k α)
    for (int i = 1; i < 3; ++i) {
        sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
        cos2kAlpha[i] = Sqr(cos2kAlpha[i - 1]) - Sqr(sin2kAlpha[i - 1]);
    }
}

Spectrum HairBSDF::f(const Vector &wo, const Vector &wi) const {
    // Compute hair coordinate system terms related to _wo_
    float sinThetaO = wo.x;
    float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
    float phiO = std::atan2(wo.z, wo.y);

    // Compute hair coordinate system terms related to _wi_
    float sinThetaI = wi.x;
    float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
    float phiI = std::atan2(wi.z, wi.y);

    // Compute $\cos \thetat$ for refracted ray
    float sinThetaT = sinThetaO / eta;
    float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

    // Compute $\gammat$ for refracted ray
    float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;//modified eta
    float sinGammaT = h / etap;
    float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
    float gammaT = SafeASin(sinGammaT);

    // Compute the transmittance _T_ of a single path through the cylinder
    Spectrum T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));

    // Evaluate hair BSDF
    float phi = phiI - phiO;
    std::array<Spectrum, pMax + 1> ap = Ap(cosThetaO, eta, h, T);
    Spectrum fsum(0.);
    for (int p = 0; p < pMax; ++p) {
        //--------------------------------------------------------------
        // Take care of the rotation angle caused by alpha tilt in hair
        //--------------------------------------------------------------
        // Compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
        float sinThetaIp, cosThetaIp;
        if (p == 0) {
            //p = 0, where θi is rotated by 2α
            sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];
            cosThetaIp = cosThetaI * cos2kAlpha[1] - sinThetaI * sin2kAlpha[1];
        }

        // Handle remainder of $p$ values for hair scale tilt
        else if (p == 1) {
            //p = 1, the rotation is by −α 
            sinThetaIp = sinThetaI * cos2kAlpha[0] - cosThetaI * sin2kAlpha[0];
            cosThetaIp = cosThetaI * cos2kAlpha[0] + sinThetaI * sin2kAlpha[0];
        } else if (p == 2) {
            //p = 2, −4α
            sinThetaIp = sinThetaI * cos2kAlpha[2] - cosThetaI * sin2kAlpha[2];
            cosThetaIp = cosThetaI * cos2kAlpha[2] + sinThetaI * sin2kAlpha[2];
        } else {
            sinThetaIp = sinThetaI;
            cosThetaIp = cosThetaI;
        }

        // Handle out-of-range $\cos \thetai$ from scale adjustment
        cosThetaIp = std::abs(cosThetaIp);
        fsum += Mp(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, v[p]) * ap[p] *
                Np(phi, p, s, gammaO, gammaT);
    }

    // Compute contribution of remaining terms after _pMax_ 
    // -------------------------------------------------------------------
    // Use uniform distribution N(φ)=1/(2π) for the azimuthal distribution
    // -------------------------------------------------------------------
    fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * ap[pMax] /
            (2.f * M_PI);
    if (AbsCosTheta(wi) > 0) fsum /= AbsCosTheta(wi);
    assert(!isinf(fsum.y()) && !isnan(fsum.y()));
    return fsum;
}

std::array<float, pMax + 1> HairBSDF::ComputeApPdf(float cosThetaO) const {
    // Compute array of $A_p$ values for _cosThetaO_
    float sinThetaO = SafeSqrt(1 - cosThetaO * cosThetaO);

    // Compute $\cos \thetat$ for refracted ray
    float sinThetaT = sinThetaO / eta;
    float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

    // Compute $\gammat$ for refracted ray
    float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
    float sinGammaT = h / etap;
    float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
    float gammaT = SafeASin(sinGammaT);

    // Compute the transmittance _T_ of a single path through the cylinder
    Spectrum T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));
    std::array<Spectrum, pMax + 1> ap = Ap(cosThetaO, eta, h, T);

    // Compute $A_p$ PDF from individual $A_p$ terms
    std::array<float, pMax + 1> apPdf;
    float sumY =
        std::accumulate(ap.begin(), ap.end(), float(0),
                        [](float s, const Spectrum &ap) { return s + ap.y(); });
    for (int i = 0; i <= pMax; ++i) apPdf[i] = ap[i].y() / sumY;
    return apPdf;
}

Spectrum HairBSDF::Sample_f(const Vector &wo, Vector *wi, const Point2f &u2,
                            float *pdf, BxDFType *sampledType) const {
    // Compute hair coordinate system terms related to _wo_
    float sinThetaO = wo.x;
    float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
    float phiO = std::atan2(wo.z, wo.y);

    // Derive four random samples from _u2_
    Point2f u[2] = {Demuxfloat(u2[0]), Demuxfloat(u2[1])};

    // Determine which term $p$ to sample for hair scattering
    std::array<float, pMax + 1> apPdf = ComputeApPdf(cosThetaO);
    int p;
    for (p = 0; p < pMax; ++p) {
        if (u[0][0] < apPdf[p]) break;
        u[0][0] -= apPdf[p];
    }

    // Sample $M_p$ to compute $\thetai$
    u[1][0] = std::max(u[1][0], float(1e-5));
    float cosTheta =
        1 + v[p] * std::log(u[1][0] + (1 - u[1][0]) * std::exp(-2 / v[p]));
    float sinTheta = SafeSqrt(1 - Sqr(cosTheta));
    float cosPhi = std::cos(2 * M_PI * u[1][1]);
    float sinThetaI = -cosTheta * sinThetaO + sinTheta * cosPhi * cosThetaO;
    float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));

    // Update sampled $\sin \thetai$ and $\cos \thetai$ to account for scales
    float sinThetaIp = sinThetaI, cosThetaIp = cosThetaI;
    if (p == 0) {
        sinThetaIp = sinThetaI * cos2kAlpha[1] - cosThetaI * sin2kAlpha[1];
        cosThetaIp = cosThetaI * cos2kAlpha[1] + sinThetaI * sin2kAlpha[1];
    } else if (p == 1) {
        sinThetaIp = sinThetaI * cos2kAlpha[0] + cosThetaI * sin2kAlpha[0];
        cosThetaIp = cosThetaI * cos2kAlpha[0] - sinThetaI * sin2kAlpha[0];
    } else if (p == 2) {
        sinThetaIp = sinThetaI * cos2kAlpha[2] + cosThetaI * sin2kAlpha[2];
        cosThetaIp = cosThetaI * cos2kAlpha[2] - sinThetaI * sin2kAlpha[2];
    }
    sinThetaI = sinThetaIp;
    cosThetaI = cosThetaIp;

    // Sample $N_p$ to compute $\Delta\phi$

    // Compute $\gammat$ for refracted ray
    float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
    float sinGammaT = h / etap;
    float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
    float gammaT = SafeASin(sinGammaT);
    float dphi;
    if (p < pMax)
        dphi =
            Phi(p, gammaO, gammaT) + SampleTrimmedLogistic(u[0][1], s, -M_PI, M_PI);
    else
        dphi = 2 * M_PI * u[0][1];

    // Compute _wi_ from sampled hair scattering angles
    float phiI = phiO + dphi;
    *wi = Vector(sinThetaI, cosThetaI * std::cos(phiI),
                   cosThetaI * std::sin(phiI));

    // Compute PDF for sampled hair scattering direction _wi_
    *pdf = 0;
    for (int p = 0; p < pMax; ++p) {
        // Compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
        float sinThetaIp, cosThetaIp;
        if (p == 0) {
            sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];
            cosThetaIp = cosThetaI * cos2kAlpha[1] - sinThetaI * sin2kAlpha[1];
        }

        // Handle remainder of $p$ values for hair scale tilt
        else if (p == 1) {
            sinThetaIp = sinThetaI * cos2kAlpha[0] - cosThetaI * sin2kAlpha[0];
            cosThetaIp = cosThetaI * cos2kAlpha[0] + sinThetaI * sin2kAlpha[0];
        } else if (p == 2) {
            sinThetaIp = sinThetaI * cos2kAlpha[2] - cosThetaI * sin2kAlpha[2];
            cosThetaIp = cosThetaI * cos2kAlpha[2] + sinThetaI * sin2kAlpha[2];
        } else {
            sinThetaIp = sinThetaI;
            cosThetaIp = cosThetaI;
        }

        // Handle out-of-range $\cos \thetai$ from scale adjustment
        cosThetaIp = std::abs(cosThetaIp);
        *pdf += Mp(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, v[p]) *
                apPdf[p] * Np(dphi, p, s, gammaO, gammaT);
    }
    *pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) *
            apPdf[pMax] * (1 / (2 * M_PI));
    // if (std::abs(wi->x) < .9999) CHECK_NEAR(*pdf, Pdf(wo, *wi), .01);
    return f(wo, *wi);
}

float HairBSDF::Pdf(const Vector &wo, const Vector &wi) const {
    // Compute hair coordinate system terms related to _wo_
    float sinThetaO = wo.x;
    float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
    float phiO = std::atan2(wo.z, wo.y);

    // Compute hair coordinate system terms related to _wi_
    float sinThetaI = wi.x;
    float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
    float phiI = std::atan2(wi.z, wi.y);

    // Compute $\cos \thetat$ for refracted ray
    float sinThetaT = sinThetaO / eta;
    float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

    // Compute $\gammat$ for refracted ray
    float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
    float sinGammaT = h / etap;
    float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
    float gammaT = SafeASin(sinGammaT);

    // Compute PDF for $A_p$ terms
    std::array<float, pMax + 1> apPdf = ComputeApPdf(cosThetaO);

    // Compute PDF sum for hair scattering events
    float phi = phiI - phiO;
    float pdf = 0;
    for (int p = 0; p < pMax; ++p) {
        // Compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
        float sinThetaIp, cosThetaIp;
        if (p == 0) {
            sinThetaIp = sinThetaI * cos2kAlpha[1] + cosThetaI * sin2kAlpha[1];
            cosThetaIp = cosThetaI * cos2kAlpha[1] - sinThetaI * sin2kAlpha[1];
        }

        // Handle remainder of $p$ values for hair scale tilt
        else if (p == 1) {
            sinThetaIp = sinThetaI * cos2kAlpha[0] - cosThetaI * sin2kAlpha[0];
            cosThetaIp = cosThetaI * cos2kAlpha[0] + sinThetaI * sin2kAlpha[0];
        } else if (p == 2) {
            sinThetaIp = sinThetaI * cos2kAlpha[2] - cosThetaI * sin2kAlpha[2];
            cosThetaIp = cosThetaI * cos2kAlpha[2] + sinThetaI * sin2kAlpha[2];
        } else {
            sinThetaIp = sinThetaI;
            cosThetaIp = cosThetaI;
        }

        // Handle out-of-range $\cos \thetai$ from scale adjustment
        cosThetaIp = std::abs(cosThetaIp);
        pdf += Mp(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, v[p]) *
               apPdf[p] * Np(phi, p, s, gammaO, gammaT);
    }
    pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) *
           apPdf[pMax] * (1 / (2 * M_PI));
    return pdf;
}

std::string HairBSDF::ToString() const {
    char tmp[256];
    sprintf(tmp,"[ Hair h: %f gammaO: %f eta: %f beta_m: %f beta_n: %f alpha: %f "
        "v[0]: %f s: %f sigma_a: ", h, gammaO, eta, beta_m, beta_n, alpha,
        v[0], s);
    return std::string(tmp)+
        sigma_a.ToString() +
        std::string("  ]");
}

Spectrum HairBSDF::SigmaAFromConcentration(float ce, float cp) {
    float sigma_a[3];
    float eumelaninSigmaA[3] = {0.419f, 0.697f, 1.37f};
    float pheomelaninSigmaA[3] = {0.187f, 0.4f, 1.05f};
    for (int i = 0; i < 3; ++i)
        sigma_a[i] = (ce * eumelaninSigmaA[i] + cp * pheomelaninSigmaA[i]);
    return Spectrum::FromRGB(sigma_a);
}

Spectrum HairBSDF::SigmaAFromReflectance(const Spectrum &c, float beta_n) {
    Spectrum sigma_a;
    for (int i = 0; i < Spectrum::nSamples_tmp; ++i)
        sigma_a[i] = Sqr(std::log(c[i]) /
                         (5.969f - 0.215f * beta_n + 2.532f * Sqr(beta_n) -
                          10.73f * Pow<3>(beta_n) + 5.574f * Pow<4>(beta_n) +
                          0.245f * Pow<5>(beta_n)));
    return sigma_a;
}