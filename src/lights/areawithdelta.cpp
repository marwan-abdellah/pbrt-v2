
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


// lights/areawithdelta.cpp*
#include "stdafx.h"
#include "lights/areawithdelta.h"
#include "paramset.h"
#include "montecarlo.h"

// DiffuseAreaLight Method Definitions
AreaWithDelta::~AreaWithDelta() {
    delete shapeSet;
}


AreaWithDelta::AreaWithDelta(const Transform &light2world,
        const Spectrum &le, int ns, const Reference<Shape> &s)
    : AreaLight(light2world, ns) {
    Lemit = le;
    shapeSet = new ShapeSet(s);
    area = shapeSet->Area();
}


Spectrum AreaWithDelta::Power(const Scene *) const {
    return Lemit * area * M_PI;
}


AreaLight *CreateDiffuseAreaLightWithDelta(const Transform &light2world, const ParamSet &paramSet,
        const Reference<Shape> &shape) {

    printf("createddddddddddddddddd");
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new AreaWithDelta(light2world, L * sc, nSamples, shape);
}


Spectrum AreaWithDelta::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    PBRT_AREA_LIGHT_STARTED_SAMPLE();
    Normal ns;

    // For a given point p, get a light sample ls and the normal on the
    // plane
    Point ps = shapeSet->Sample(p, ls, &ns);

    // Gets the direction between the point on the light source and the
    // point on the ray
    *wi = Normalize(ps - p);
    // *wi =  Vector(0, 1, 0); // This restricts the direction to be normal only

    // Calculate the pdf
    *pdf = shapeSet->Pdf(p, *wi);

    // Check if there is no blocker between the sample point on the ray
    // and the sample point on the light source
    visibility->SetSegment(p, pEpsilon, ps, 1e-3f, time);

    // Retrun the radiance from that sample
    // Check for the right directionality of the light
    Spectrum Ls = L(ps, ns, -*wi);

    PBRT_AREA_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


float AreaWithDelta::Pdf(const Point &p, const Vector &wi) const {
    return shapeSet->Pdf(p, wi);
}


Spectrum AreaWithDelta::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_AREA_LIGHT_STARTED_SAMPLE();
    Point org = shapeSet->Sample(ls, Ns);
    Vector dir;
    dir.x = Ns->x;
    dir.y = Ns->y;
    dir.z = Ns->z;

    UniformSampleSphere(u1, u2);
    if (Dot(dir, *Ns) < 0.) dir *= -1.f;
    *ray = Ray(org, dir, 1e-3f, INFINITY, time);
    *pdf = shapeSet->Pdf(org) * INV_TWOPI;
    Spectrum Ls = L(org, *Ns, dir);
    PBRT_AREA_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


