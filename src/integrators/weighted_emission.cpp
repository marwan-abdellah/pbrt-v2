
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


// integrators/weighted_emission.cpp*
#include "stdafx.h"
#include "integrators/weighted_emission.h"
#include "volumes/volumegrid.h"
#include "paramset.h"
#include "iostream"
#include <typeinfo>

static int PRINT;

// WeightedEmissionIntegrator Method Definitions
void WeightedEmissionIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
                                        const Scene *scene) {
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}


Spectrum WeightedEmissionIntegrator::Transmittance(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
    if (!scene->volumeRegion) return Spectrum(1.f);
    float step, offset;
    if (sample) {
        step = stepSize;
        offset = sample->oneD[tauSampleOffset][0];
    }
    else {
        step = 4.f * stepSize;
        offset = rng.RandomFloat();
    }
    Spectrum tau = scene->volumeRegion->tau(ray, step, offset);
    return Exp(-tau);
}


Spectrum WeightedEmissionIntegrator::Li(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray,
        const Sample *sample, RNG &rng, Spectrum *T,
        MemoryArena &arena) const {

    PRINT = 0;
    VolumeGridDensity *vr = static_cast<VolumeGridDensity*>(scene->volumeRegion);

#ifdef DISPLAY_BB
    // Data Bounding Box
    BBox BoundingBox= vr->WorldBound();
    printf("BB-MIN::[%f, %f, %f] \n", BoundingBox.pMin.x, BoundingBox.pMin.y, BoundingBox.pMin.z);
    printf("BB-MAX::[%f, %f, %f] \n", BoundingBox.pMax.x, BoundingBox.pMax.y, BoundingBox.pMax.z);
#endif

    Assert(sample != NULL);
    float t0, t1;
    if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f) {
        *T = Spectrum(1.f);
        return 0.f;
    }
    // Do "weighted emission"-only volume integration in _vr_
    Spectrum Lv(0.);
    Spectrum ShiftSpectrum(0.);

    // Prepare for volume integration stepping
    int nSamples = Ceil2Int((t1-t0) / stepSize);
    float step = (t1 - t0) / nSamples;
    Spectrum Tr(1.f);
    Point p = ray(t0), pPrev;
    Vector w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step;
    for (int i = 0; i < nSamples; ++i, t0 += step) {

        // Advance to sample at _t0_ and update _T_
        pPrev = p;
        p = ray(t0);

        // std::cout << "p::x " << p.x << std::endl;
        // std::cout << "p::y " << p.y << std::endl;
        // std::cout << "p::z " << p.z << std::endl;

        Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);
        Spectrum stepTau = vr->tau(tauRay,
                                   .5f * stepSize, rng.RandomFloat());

        Tr *= Exp(-stepTau);

        // Possibly terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) {
                Tr = 0.f;
                break;
            }
            Tr /= continueProb;
        }

        // Compute emission-only source term at _p_
        Lv += Tr * vr->Lve(p, w, ray.time);

#ifdef RGBSpectrum
        // Spectrum color components RGB
        // printf("Spectrum RGB[%f, %f, %f] \n", Lv.getR(), Lv.getG() , Lv.getB());
#endif
    }
    *T = Tr;



#ifdef RGBSpectrum
    ShiftSpectrum = Lv;
    ShiftSpectrum.setG(Lv.getB());
    ShiftSpectrum.setB(Lv.getG());
    return ShiftSpectrum * step;
#endif


    // Printing the spectrum
    // for (int iSample = 0; iSample < 30; ++iSample)
    // {   if (Lv.getC()[iSample] != 0)
    //        printf("The SPD[%d] = %f \n", iSample, Lv.getC()[iSample]);
    // }


    return Lv * step;
}


WeightedEmissionIntegrator *CreateWeightedEmissionVolumeIntegrator(const ParamSet &params) {
    float stepSize  = params.FindOneFloat("stepsize", 1.f);
    std::cout << std::endl;
    std::cout << "WeightedEmissionIntegrator ................." << std::endl;
    return new WeightedEmissionIntegrator(stepSize);
}


