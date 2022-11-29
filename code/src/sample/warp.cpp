/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Vector3f Warp::sampleUniformHemisphere(Sampler *sampler, const Normal3f &pole) {
    // Naive implementation using rejection sampling
    Vector3f v;
    do {
        v.x() = 1.f - 2.f * sampler->next1D();
        v.y() = 1.f - 2.f * sampler->next1D();
        v.z() = 1.f - 2.f * sampler->next1D();
    } while (v.squaredNorm() > 1.f);

    if (v.dot(pole) < 0.f)
        v = -v;
    v /= v.norm();

    return v;
}

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float r = sqrt(sample[0]), phi = 2 * M_PI * sample[1];
    return Point2f(r * cos(phi), r * sin(phi));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
   return float(p.norm() <= 1) * INV_PI;
}

Vector3f Warp::squareToUniformCylinder(const Point2f &sample){
    float z = 2 * sample[0] - 1, phi = 2 * M_PI * sample[1];
    return Vector3f(cos(phi), sin(phi), z);
}

Vector3f Warp::squareToUniformSphereCap(const Point2f &sample, float cosThetaMax) {
    return squareToUniformSphere(Point2f((sample[0] * (1 - cosThetaMax) + 1 + cosThetaMax)/2.f, sample[1]));
}

float Warp::squareToUniformSphereCapPdf(const Vector3f &v, float cosThetaMax) {
    return float(abs(v.norm() - 1) < Epsilon && v[2] - cosThetaMax > Epsilon) * 0.5 * INV_PI / (1 - cosThetaMax);
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    Vector3f cylinder_sample = squareToUniformCylinder(sample);
    float r = sqrt(1 - cylinder_sample[2] *  cylinder_sample[2]);
    return Vector3f(r * cylinder_sample[0], r * cylinder_sample[1], cylinder_sample[2]);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return float(abs(v.norm() - 1) < Epsilon) * 0.25 * INV_PI;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    return squareToUniformSphere(Point2f((sample[0] + 1)/2.f, sample[1]));
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return float(abs(v.norm() - 1) < Epsilon && v[2] > Epsilon) * 0.5 * INV_PI;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Point2f disk_sample = squareToUniformDisk(sample);
    return Vector3f(disk_sample[0], disk_sample[1], sqrt(1 - disk_sample[0] * disk_sample[0] - disk_sample[1] * disk_sample[1]));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    return float(abs(v.norm() - 1) < Epsilon && v[2] > Epsilon) * v[2] * INV_PI;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float theta = atan(alpha * sqrt(-log(1 - sample[0]))), phi = 2 * M_PI * sample[1];
    return Vector3f(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    float theta = acos(m[2]), alpha_square = alpha * alpha;
    return float(abs(m.norm() - 1) < Epsilon && m[2] > Epsilon) * exp( -pow(tan(theta), 2) / alpha_square ) / (M_PI * alpha_square * pow(cos(theta), 3));
}

Vector3f Warp::squareToUniformTriangle(const Point2f &sample) {
    float su1 = sqrtf(sample.x());
    float u = 1.f - su1, v = sample.y() * su1;
    return Vector3f(u,v,1.f-u-v);
}

Vector3f Warp::squareToHenyeyGreenstein(const Point2f &sample, float g) {

    //use the inverse method
    if(g == 0.0){
        return Warp::squareToUniformSphere(sample);
    }
    float fraction = (1.0f - g*g) / (1.0f - g + 2.0f * g * sample[0]);
    float theta =  acos(1.0f / (2.0f * g) * (1.0f +  g * g - fraction * fraction));
    float phi = 2.0 * M_PI * sample[1];
    return Vector3f(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
}

float Warp::squareToHenyeyGreensteinPdf(const Vector3f &m, float g){
    float cos_theta = m[2];
    float frac = 1.0f + g * g - 2.0f * g * cos_theta;
    if(frac <= 0.0f) return 0.0f;
    float pdf =  ((1.0f - g * g) / (4.0f * M_PI * pow(frac, 1.5f)));
    return pdf;
}

NORI_NAMESPACE_END
