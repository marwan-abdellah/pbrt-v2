
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
#include "volumes/volumegrid.h"

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

// RADIANCE CALCULATION
Spectrum Abdellah::Li(const Scene *scene,
                      const Renderer *renderer,
                      const RayDifferential &ray,
                      const Sample *sample, RNG &rng,
                      Spectrum *T,
                      MemoryArena &arena) const {

    // This is our volume region to be integrated in the scene...
    VolumeRegion *vr = scene->volumeRegion;
    // DensityRegion *vr = dynamic_cast<DensityRegion*> (scene->volumeRegion);


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

        // The new point "p" which is refered by "p'" in the equation.
        p = ray(t0);

        // This is the ray that will advance in the volume
        Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);

        // Calculate tau in the volume
        // Step size is reduced by half ???
        // May be for increasing quality
        // Offset is randomized ???
        // The first sample is placed randomly in the first segment "p.881"
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

        // Compute single-scattering source term at _p_
        // This is basically due to the emission as calculated before
        // First term in the transfer equation (16.4)
        Lv += Tr * vr->Lve(p, w, ray.time);

        // Calculate the second part of equation in section 16.4
        // The direct contribution from the light source RADIANCE
        // Get sigma scatering (ss)
        // Make some checkes for reducing the rendering time
        // * Check if the volume is not scattering (black body which only absorbs)
        // * Check if the scene doesn't have any light source (DARK scene)
        // *** Consider only a single light source in the scene ***
        // * Check if light is occluded or not

        // Get the scattering coeffecient sigma_s
        // Defined in the scene descriptor
        Spectrum ss = vr->sigma_s(p, w, ray.time);

        // Checking if the scattering coeff is != 0 (BLACK BODY)
        // Neither the there are no lights in the scene
        if (!ss.IsBlack() && scene->lights.size() > 0) {

            // Number of lights in the scene
            int nLights = scene->lights.size();

            // ???
            // Light index
            // Getting a single light from the light sources in the scene ???
            // Selecting a light to sample
            // A single light is selected and sampled at each point along the ray
            // The light contribution is scaled by the number of lights
            // Sample one light strategy
            int ln = min(Floor2Int(lightNum[sampOffset] * nLights),
                         nLights-1);

            // Getting light in the scene
            Light *light = scene->lights[ln];

            // Add contribution of _light_ due to scattering at _p_
            float pdf;

            // Checking the obstructions between the light and the point
            VisibilityTester vis;

            // Observer direction
            Vector wo;

            // Light sample
            LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
                           lightPos[2*sampOffset+1]);

            // The second term of single scattering equation in section 16.4
            // Taking a light sample on the "selected" light source
            // Calculate the radiance at point "p" (should be p')
            // It returns the distribution of light pdf at (p, wo)
            // Light only contributes to the radiance if pdf > 0
            Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
            
            // Check if the spectrum is black
            // Check if the light has contribution along the ray
            // Check if the light is unoccluded
            // If yes, then add that contribution along the ray
            if (!L.IsBlack() && pdf > 0.f && vis.Unoccluded(scene)) {

                // Ld at the point along the ray "p'"
                Spectrum Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);

                int ELASTIC_SCATTER = 0;
                int INELASTIC_SCATTER = 1;

                // Elastic Scattering
                if(ELASTIC_SCATTER)
                {
                    // The multiplication with the number of light is a bit confusing
                    // The division by the pdf is a bit confusing
                    // This according to the "sample one light strategy"
                    Lv += Tr * ss * vr->p(p, w, -wo, ray.time) * Ld * float(nLights) /
                            pdf;
                }

                // Inelastic scattering
                // Spectral shift
                if (INELASTIC_SCATTER)
                {
                    // Get the contribution from the in elastic scattering and use the
                    // isotropic phase function

                    // Consider the radiance due to the light "Ld"
                    // In that part, I should account for the spectral shift with some spectral
                    // shifting function and also multiply with the quantum yield
                    // Qunatum yield is obtained from Wikipedia for GFP
                    float quantumYield = 0.79;

                    // Flu. radiance
                    Spectrum Lfluro(0.);

                    // Radiance due to fluorescence along the ray at point "p'"
                    // Tr is less than 1
                    // Tr is a function of distance
                    Lfluro = quantumYield * Tr * ss * (1 / (4.f * M_PI)) * Ld * float(nLights) /
                            pdf;

                    // Do the spectral shift
                    float* sampleValues = Lfluro.getSapectrumSamples();
                    float energyCount = 0;
                    for (int i = 0; i < nSpectralSamples; i++)
                    {
                        energyCount += sampleValues[i];
                        sampleValues[i] = 0;
                    }

                    // Zero Spectrum
                    Lfluro.zeroSpectrum();

                    // Shift it to the selected wavelength region
                    // 11 corresponds to 510 nm
                    sampleValues[11] = energyCount/1;

                    Lfluro.setSpectrumValues(sampleValues);
                    Lv += Lfluro;

                    // Conservation of enegry is needed when we account for both cases
                    if (ELASTIC_SCATTER)
                        Lv /= 2;

                    free(sampleValues);
                }
            }
        }

        ++sampOffset;
    }
    *T = Tr;

    return Lv * step;
}




Abdellah *CreateAbdellahIntegrator(const ParamSet &params) {
    float stepSize  = params.FindOneFloat("stepsize", 1.f);
    return new Abdellah(stepSize);
}


