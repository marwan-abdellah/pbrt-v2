
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
//#include "3rdparty/ilmbase-1.0.2/"

#include "3rdparty/ilmbase-1.0.2/ImathMatrix.h"
#include "3rdparty/ilmbase-1.0.2/ImathMatrixAlgo.h"

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

        // Compute plane intersection for disk
        // Checks if the plane is parallel to the ray or not
        // We can get the direction of the ray
        // If the Z component of the direction of the ray is zero
        // then, the ray is parallel to the plane and in such case
        // there is no intersection point between the ray and the plane.
        if (fabsf(ray.d.z) < 1e-7)
            return false;

        // Now, the direction of the ray is not parallel to the plane
        // We have to check if the intersection happens or not
        // We have to compute the parametric t where the ray intersects the plane
        // We want to find t such that the z-component of the ray intersects the plane
        // The ray "line" equation is l = l0 + (l1 - l0) * t
        // l1 - l0 will give us the distance between the two points on the plane
        // Then t is the ratio and in such case it should be between 0 and 1
        // Considering that the rectangle completely lies in the z plane
        /// distance = l1 - l0
        /// thit = (l - l0) / distance

        // But since we assume that the plane is located at height
        // Then, the point l is at height on the plane
        /// l = height
        float thit = (height - ray.o.z) / ray.d.z;

        // Then we check if the thit is between the ratio of 0 and 1 that is mapped
        // between ray.mint and ray.maxt, if not retrun false
        if (thit < ray.mint || thit > ray.maxt)
            return false;

        // Then we see if the point lies inside the disk or not
        // Substitute the thit in the ray equation to get hit point on the ray
        Point phit = ray(thit);

        // We have to make sure that the interesction lies inside the plane
        if (!(phit.x < x/2 && phit.x > -x/2 && phit.y < y/2 && phit.y > -y/2))
            return false;

        // Assuming that the plane is formed from the following 4 points
        // P0, P1, P2, P3
        //
        // p0 *---------------* p1
        //    |               |
        //    |               |
        //    |               |
        //    |       O       |
        //    |               |
        //    |               |
        //    |               |
        // p2 *---------------* p3  -> X
        //
        // P0 @ (-x/2, y/2)
        // P1 @ (x/2, y/2)
        // P2 @ (-x/2, -y/2)
        // P3 @ (x/2, -y/2)
        Point P0(-x/2, y/2, height), P1(x/2, y/2, height);
        Point P2(-x/2, -y/2, height), P3(x/2, -y/2, height);

        /// Now, we have to find the parametric form of the plane in terms of (u,v)
        /// Plane equation can be formed by at least 3 points P0, P1, P2
        /// P0 -> P1 (vector 1)
        /// P0 -> p2 (vector 2)
        /// An arbitrary point on the plane p is found in the following parametric form
        /// P = P0 + (P1 - P0) u + (P2 - P0) v
        /// Now we need to express two explicit equations of u and v
        /// So, we have to construct the system of equation that solves for u and v
        ///
        /// Since we have found the intersection point between the plane and the line
        /// we have to use it to formalize the system of equations that will be used
        /// to find the parametric form of the plane
        /// Plane equation is : P = P0 + (P1 - P0) u + (P2 - P0) v
        /// Ray equation is : l = l0 + (l1 - l0) * thit
        /// But l = P, then
        /// l0 + (l1 - l0) * thit = P0 + (P1 - P0) * u + (P2 - P0) * v
        /// l0 - P0 = (l0 - l1) * thit +  (P1 - P0) * u + (P2 - P0) * v
        /// MAPPING : l0 = ray.o
        /// [l0.x - P0.x] = [l0.x - l1.x P1.x - P0.x P2.x - P0.x] [t]
        /// [l0.y - P0.y] = [l0.y - l1.y P1.y - P0.y P2.y - P0.y] [u]
        /// [l0.z - P0.z] = [l0.z - l1.z P1.z - P0.z P2.z - P0.z] [v]
        ///
        /// Then, we should find the inverse of the matrix in order to
        /// solve for u,v and t for check

        // System AX = B
        float a11 = ray.o.x - 0;
        float a12 = P1.x - P0.x;
        float a13 = P2.x - P0.x;
        float a21 = ray.o.y - 0;
        float a22 = P1.y - P0.y;
        float a23 = P2.y - P0.y;
        float a31 = ray.o.y - height;
        float a32 = P1.z - P0.z;
        float a33 = P2.z - P0.z;

        float b1 = -7;
        float b2 = -2;
        float b3 = 14;

        float x1 = 0;
        float x2 = 0;
        float x3 = 0;

        Imath::M33f A(a11,a12,a13,a21,a22,a23,a31,a32,a33), AInverted;
        Imath::V3f X(x1, x2, x3);
        Imath::V3f B(b1,b2, b3);

        // This operation has been checked and working for getting
        // the correct inverse of the matrix A
        AInverted = A.invert(false);

        x1 =  AInverted[0][0] * B[0] +
                AInverted[0][1] * B[1] +
                AInverted[0][2] * B[2];

        x2 =  AInverted[1][0] * B[0] +
                AInverted[1][1] * B[1] +
                AInverted[1][2] * B[2];

        x3 =  AInverted[2][0] * B[0] +
                AInverted[2][1] * B[1] +
                AInverted[2][2] * B[2];

        /// Then we have u = something, and v = something
        ///
        /// Then we come for the derivatives, so we have to find the derivatives
        /// from the parametric forms defined above for the plane equations
        /// dpdu = (P1 - P0)
        /// dpdv = (P2 - P0)
        ///
        /// For the normal we have the always fixed direction in y
        /// So the derivative for the normal is zero
        /// dndu = (0, 0, 0) dndv = (0, 0, 0)
        ///
        /// Then we can construct the DifferentilGeometry and go ahead

        // Find parametric representation of disk hit
        float u = x2;
        float v = x3;

        Vector dpdu(P1.x - P0.x, P1.y - P0.y, P1.z - P0.z);
        Vector dpdv(P2.x - P0.x, P2.y - P0.y, P2.z - P0.z);
        Normal dndu(0,0,0), dndv(0,0,0);


        // Initialize _DifferentialGeometry_ from parametric information
        const Transform &o2w = *ObjectToWorld;
        *dg = DifferentialGeometry(o2w(phit), o2w(dpdu), o2w(dpdv),
                                   o2w(dndu), o2w(dndv), u, v, this);

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


Point Rectangle::Sample(float u, float v, Normal *Ns) const {

    // u and v are rendom samples on the surface of the light source
    Point P0(-x/2, y/2, height), P1(x/2, y/2, height);
    Point P2(-x/2, -y/2, height), P3(x/2, -y/2, height);

    Point p = P0 + (P1 - P0) * u + (P2 - P0) * v;

    Normal n = Normal(Cross(P2-P0, P1-P0));
    *Ns = Normalize(n);


    // NORMAL ON THE PLANE
    *Ns = Normalize((*ObjectToWorld)(n));
    if (ReverseOrientation) *Ns *= -1.f;
    return (*ObjectToWorld)(p);
}


float Rectangle::Pdf(const Point &p, const Vector &wi) const
{
    // Intersect sample ray with area light geometry
    DifferentialGeometry dgLight;
    Ray ray(p, wi, 1e-3f);
    ray.depth = -1; // temporary hack to ignore alpha mask
    float thit, rayEpsilon;
    if (!Intersect(ray, &thit, &rayEpsilon, &dgLight)) return 0.;

    float pdf;
    Point pObject = (*WorldToObject)(p);

    // TODO: Check with the normal (n.p)
    if (!(pObject.x < x/2 && pObject.x > -x/2
            && pObject.y < y/2 && pObject.y > -y/2))
        return 0;

    // N & wi

    Point P0(-x/2, y/2, height);
    Point P1(x/2, y/2, height);
    Point P2(-x/2, -y/2, height);

    Point P0_W = (*WorldToObject)(P0);
    Point P1_W = (*WorldToObject)(P1);
    Point P2_W = (*WorldToObject)(P2);


    Normal n = Normal(Cross(P2_W-P0_W, P1_W-P0_W));
    Normal Ns = Normalize(n);

    Normal Nwi(wi.x, wi.y, wi.z);

    if (Dot(Ns,  Nwi) < 0.0010)
        pdf = 1.0;

    if (isinf(pdf)) pdf = 0.f;
        return pdf;
}
