
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


// shapes/rectangle.cpp*
#include "stdafx.h"
#include "shapes/rectangle.h"
#include "paramset.h"
#include "montecarlo.h"

// Rectangle Method Definitions
Rectangle::Rectangle(const Transform *o2w, const Transform *w2o, bool ro,
           float x_, float y_, float height_)
    : Shape(o2w, w2o, ro) {

    // X * Y
    x = x_;
    y = y_;

    // WHERE THE PLAN OF THE RECTANGLE INTERSECTS
    height = height_;
}


BBox Rectangle::ObjectBound() const {

    // BOUNDING BOX WILL BE IN THE LOCAL COORDINATES
    return BBox(Point(-x/2, -y/2, height),
                Point( x/2,  y/2, height));
}


bool Rectangle::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
                     DifferentialGeometry *dg) const {

    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);

    // Compute plane intersection with the rectangle
    // RETURN FALSE IF THE RAY IS PARALLEL TO THE PLANE
    // I.E. THE Z COMPONENT IS ALMOST EQUAL TO ZERO
    if (fabsf(ray.d.z) < 1e-7)
        return false;

    // CHECK IF THE Z-COMPONENT OF THE RAY IS EQUAL TO THE HEIGHT
    float thit = (height - ray.o.z) / ray.d.z;
    if (thit < ray.mint || thit > ray.maxt)
        return false;

    // NOW WE WANT TO CALCULATE phit WHERE THE RAY INTERSECTS THE PLANE
    Point phit = ray(thit);

    // IF phit IS OUTSIDE THE RECTANGLE THEN RETURN FALSE
    if (!((phit.x < x/2 && phit.y < y/2) ||
            (phit.x > -x/2 && phit.y < y/2) ||
            (phit.x > -x/2 && phit.y > y/2) ||
            (phit.x < x/2 && phit.y > -y/2)))
        return false;

    // The normal is the same over the rectangle
    Normal dndu(0,0,0), dndv(0,0,0);

    // Initialize _DifferentialGeometry_ from parametric information
    const Transform &o2w = *ObjectToWorld;
    // *dg = DifferentialGeometry(o2w(phit), NULL, NULL,
       //                        o2w(dndu), o2w(dndv), 0, 0, this);

    // Update _tHit_ for quadric intersection
    *tHit = thit;

    // Compute _rayEpsilon_ for quadric intersection
    *rayEpsilon = 5e-4f * *tHit;

    return true;
}


bool Rectangle::IntersectP(const Ray &r) const {

    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);

    // Compute plane intersection for disk
    if (fabsf(ray.d.z) < 1e-7) return false;
    float thit = (height - ray.o.z) / ray.d.z;
    if (thit < ray.mint || thit > ray.maxt)
        return false;

    // NOW WE WANT TO CALCULATE phit WHERE THE RAY INTERSECTS THE PLANE
    Point phit = ray(thit);

    // IF phit IS OUTSIDE THE RECTANGLE THEN RETURN FALSE
    if (!((phit.x < x/2 && phit.y < y/2) ||
            (phit.x > -x/2 && phit.y < y/2) ||
            (phit.x > -x/2 && phit.y > y/2) ||
            (phit.x < x/2 && phit.y > -y/2)))
        return false;

    return true;
}


float Rectangle::Area() const {
    return x * y;
}


Rectangle *CreateRectangleShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    float x = params.FindOneFloat("x", 1.);
    float y = params.FindOneFloat("y", 1.);
    float height = params.FindOneFloat("height", 0.);
    return new Rectangle(o2w, w2o, reverseOrientation,
                         x, y, height);
}


Point Rectangle::Sample(float u1, float u2, Normal *Ns) const {

    // SAMPLING THE RECTANGLE TO RETURN A POINT ON IT
    Point p;

    ConcentricSampleDisk(u1, u2, &p.x, &p.y);

    p.x *= x/2;
    p.y *= y/2;
    p.z = height;

    // NORMAL ON THE PLANE
    *Ns = Normalize((*ObjectToWorld)(Normal(0,0,1)));
    if (ReverseOrientation) *Ns *= -1.f;
    return (*ObjectToWorld)(p);
}


