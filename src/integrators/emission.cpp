
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


// integrators/emission.cpp*
#include "stdafx.h"
#include "integrators/emission.h"
#include "paramset.h"

// EmissionIntegrator Method Definitions
/**
 *
 */
void EmissionIntegrator::RequestSamples(Sampler *sampler,
                                        Sample *sample,
                                        const Scene *scene) {
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}


Spectrum EmissionIntegrator::Transmittance(const Scene *scene,
                                           const Renderer *renderer,
                                           const RayDifferential &ray,
                                           const Sample *sample,
                                           RNG &rng,
                                           MemoryArena &arena) const {
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

// this methid is responsible for calculating the second term in equation 16.2
// it doesn't take into account the first term that is resulting from the
// intersection of the rays with any surface in the scene.
Spectrum EmissionIntegrator::Li(const Scene *scene, // scene
                                const Renderer *renderer, // renderer
                                const RayDifferential &ray, // for each ray
                                const Sample *sample, // sample along the ray
                                RNG &rng,
                                Spectrum *T, // Transmittance along the ray
                                MemoryArena &arena) const {

    // get the volume from the scene
    VolumeRegion *vr = scene->volumeRegion;

    // check that there some samples along the ray,
    // otherwise exit
    Assert(sample != NULL);

    // t0 -> ray::mint
    // t1 -> tay::maxt
    float t0, t1;

    // calculate the intersection points between the ray and the volume
    // retrun 0 if no volume or no intersection between the volume and the ray
    // this means that we don't have any addition from that volume to the ray
    if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f) {
        *T = Spectrum(1.f);
        return 0.f;
    }

    // do emission-only volume integration in _vr_
    // Lv is the calculated value at the ray origin due to the emissio  n
    Spectrum Lv(0.);

    // prepare for volume integration stepping
    // for approximation
    int nSamples = Ceil2Int((t1-t0) / stepSize);

    // calculating the real step size
    float step = (t1 - t0) / nSamples;

    // used for summating up the transmiattance value accross the ray
    Spectrum Tr(1.f);

    // points along the ray
    // current and previous
    Point p = ray(t0), pPrev;

    // ray direction
    Vector w = -ray.d;

    //TODO !!!
    t0 += sample->oneD[scatterSampleOffset][0] * step;

///////////////////////////////////////////////////////////////////////////////////////////
    // loop over the samples of the ray and calculate the radinace added by each sample
    // to eventually calculate the radiance at the origin of the ray "pixel"
    for (int i = 0; i < nSamples; ++i, t0 += step) {
        // Advance to sample at _t0_ and update _T_
        pPrev = p;
        p = ray(t0);

        // a ray called tauRay to calculate the transmittace between 2 successive points
        Ray tauRay(pPrev,       // start
                   p - pPrev,   // direction vector
                   0.f,         // start
                   1.f,         // end
                   ray.time,    // time
                   ray.depth);  // depth

        // calcaulate transmittance between two points having the ray and the volume
        Spectrum stepTau = vr->tau(tauRay,              // ray in the volume
                                   .5f * stepSize,      // TODO "0.5" step !!!
                                   rng.RandomFloat());  // offset

        // T(A->C) = T(A->B)T(B->C)
        // that's why it is initialized by "1"
        Tr *= Exp(-stepTau);

///////////////////////////////////////////////////////////////////////////////////////////
        // Possibly terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) {
                Tr = 0.f;
                break;
            }
            Tr /= continueProb;
        }
///////////////////////////////////////////////////////////////////////////////////////////

        // Compute emission-only source term at _p_
        // we have the transmittance from the volume uptil the point p
        // this gives us flexibelity to calculate the source term
        // Lve * Tr would give us the
        Lv += Tr * vr->Lve(p, w, ray.time);
    }
///////////////////////////////////////////////////////////////////////////////////////////

    // return the transmittance along the whole volume to be used later
    *T = Tr;

    // return the Li along the ray in the volume
    // TODO multiplication by the step
    return Lv * step;
}


EmissionIntegrator *CreateEmissionVolumeIntegrator(const ParamSet &params) {
    float stepSize  = params.FindOneFloat("stepsize", 1.f);
    return new EmissionIntegrator(stepSize);
}


