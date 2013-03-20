
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

// this is a testing integrator for the fluorescence effects
// integrators/abdellah.cpp*
#include "stdafx.h"
#include "integrators/abdellah.h"
#include "scene.h"
#include "paramset.h"
#include "montecarlo.h"

// Abdellah Method Definitions
// Get the samples along the ray
void Abdellah::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}

//Calculate the transmittance
Spectrum Abdellah::Transmittance(const Scene *scene,
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






void Abdellah::Test(){


    float lambdaInitial = 400;
    float lambdaFinal = 700;
    int lambdaSamples = 30;
    float lambdaStep = (lambdaFinal - lambdaInitial) / (lambdaSamples);

    int lambdaCtr;



    for (lambdaCtr = 0; lambdaCtr < lambdaSamples; lambdaCtr++)
    {

        // Power value
        float powerValue;




    }


}








// RADIANCE CALCULATION
Spectrum Abdellah::Li(const Scene *scene,
                      const Renderer *renderer,
                      const RayDifferential &ray,
                      const Sample *sample, RNG &rng,
                      Spectrum *T,
                      MemoryArena &arena) const {

    // This is our volume region to be integrated in the scene...
    VolumeRegion *vr = scene->volumeRegion;


    // Intersection points between the ray and the volume
    float t0, t1;

    // if NO VOLUME
    // if VOLUME DOESN'T INTERSECT THE RAY
    // if NO INTERSECTION
    if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f) {
        *T = 1.f;
        return 0.f;
    }

    // Do single scattering volume integration in _vr_
    Spectrum Lv(0.);

    // Prepare for volume integration stepping
    int nSamples = Ceil2Int((t1-t0) / stepSize);
    float step = (t1 - t0) / nSamples;
    Spectrum Tr(1.f);
    Point p = ray(t0), pPrev;
    Vector w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step;

    // Compute sample patterns for single scattering samples
    float *lightNum = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightNum, rng);
    float *lightComp = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightComp, rng);
    float *lightPos = arena.Alloc<float>(2*nSamples);
    LDShuffleScrambled2D(1, nSamples, lightPos, rng);
    uint32_t sampOffset = 0;



    // Tracing ray
    for (int i = 0; i < nSamples; ++i, t0 += step) {

        // Advance to sample at _t0_ and update _T_
        pPrev = p;

        p = ray(t0);

        // Calculate tau along the ray s
        Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);

        // Calculate tau in the volume
        Spectrum stepTau = vr->tau(tauRay,
                                   .5f * stepSize, rng.RandomFloat());

        // Calculate the trnasmittance
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

        // Compute single-scattering source term at _p_
        // This is basically due to the emissio as calculated before
        Lv += Tr * vr->Lve(p, w, ray.time);

        // Get the scattering coeffecient sigma_s
        // Defined in the scene descriptor
        Spectrum ss = vr->sigma_s(p, w, ray.time);

        // Checking if the scattering coeff is != 0 (BLACK BODY)
        // Neither the there are no lights in the scene
        if (!ss.IsBlack() && scene->lights.size() > 0) {

            // Number of lights in the scene
            int nLights = scene->lights.size();

            // ???
            int ln = min(Floor2Int(lightNum[sampOffset] * nLights),
                         nLights-1);

            // Getting light in the scene
            Light *light = scene->lights[ln];

            // Add contribution of _light_ due to scattering at _p_
            float pdf;

            // Checking the obstructions between the light and the point
            VisibilityTester vis;
            Vector wo;

            // light sample
            LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
                           lightPos[2*sampOffset+1]);

            // Getting the contribution of the light
            Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
            
            // Light has an effet "check"
            if (!L.IsBlack() && pdf > 0.f && vis.Unoccluded(scene)) {

                // Get the contribution due to the light source for the elsatic scattering
                Spectrum Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);

                // Elastic Scattering
                if(1)
                {
                    Lv += Tr * ss * vr->p(p, w, -wo, ray.time) * Ld * float(nLights) /
                            pdf;
                }



                // Inelastic scattering
                if (1)
                {
                    // Get the contribution from the in elastic scattering and use the
                    // isotropic phase function

                    // Radiance due to the light
                    // In that part, I should account for the spectral shift with some spectral
                    // shifting function and also multiply with the quantum yield
                    float quantumYield = 0.79;
                    Spectrum Lfluro(0.);

                    Lfluro = quantumYield * Tr * ss * (1 / (4.f * M_PI)) * Ld * float(nLights) /
                            pdf;

                    // Do the spectral shift
                    float* sampleValues = Lfluro.getSapectrumSamples();
                    float energyCount = 0;
                    for (int i = 0; i < nSamples; i++)
                    {
                        energyCount += sampleValues[i];
                        sampleValues[i] = 0;
                    }

                    // Zero Spectrum
                    Lfluro.zeroSpectrum();

                    // Shift it to the selected wavelength region

                    sampleValues[9] = energyCount/3;
                    sampleValues[10] = energyCount/3;
                    sampleValues[11] = energyCount/3;

                    Lfluro.setSpectrumValues(sampleValues);
                    Lv += Lfluro;

                    // Conservation of enegry
                    Lv /= 2;

                    free(sampleValues);
                }
            }
        }

        ++sampOffset;
    }
    *T = Tr;

//    float* sampleValues = Lv.getSapectrumSamples();

//    // Get all the energy of the spectrum
//    // and add it to the green components only
//    float energyCount = 0;
//    for (int i = 0; i < nSamples; i++)
//        energyCount += sampleValues[i];

//    // Zero Spectrum
//    Lv.zeroSpectrum();


//    // Green wavelengths : 530 [13] - 600 [20]
//    const int lambda_start = 18;
//    const int lambda_end = 19;
//    for (int i = lambda_start; i < lambda_end; i++)
//    {
//        int N = lambda_end - lambda_start;
//        sampleValues[i] = energyCount / N;
//    }

//    Lv.setSpectrumValues(sampleValues);

//    free(sampleValues);

    return Lv * step;
}




Abdellah *CreateAbdellahIntegrator(const ParamSet &params) {
    float stepSize  = params.FindOneFloat("stepsize", 1.f);
    return new Abdellah(stepSize);
}


