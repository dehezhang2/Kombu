/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Romain Prévost

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

#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN

class AreaEmitter : public Emitter {
public:
    AreaEmitter(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
    }

    virtual std::string toString() const override {
        return tfm::format(
                "AreaLight[\n"
                "  radiance = %s,\n"
                "]",
                m_radiance.toString());
    }

    virtual Color3f eval(const EmitterQueryRecord & lRec) const override {
        if(!m_shape)
            throw NoriException("There is no shape attached to this Area light!");
        return (lRec.n.dot(-lRec.wi) > 0 ? m_radiance : Color3f(0.f));
    }

    virtual Color3f sample(EmitterQueryRecord & lRec, const Point2f & sample) const override {
        if(!m_shape)
            throw NoriException("There is no shape attached to this Area light!");
        ShapeQueryRecord sRec(lRec.ref);
        m_shape->sampleSurface(sRec, sample);
        // light source
        lRec.p = sRec.p;
        // normal vector of the light source
        lRec.n = sRec.n;
        // direction from intersection to the light source
        lRec.wi = (lRec.p - lRec.ref).normalized();
        lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, (lRec.p - lRec.ref).norm() - Epsilon);
        lRec.pdf = pdf(lRec);
        if(lRec.pdf > 0.f ) return eval(lRec) / lRec.pdf;
        return Color3f(0.f);
    }

    virtual float pdf(const EmitterQueryRecord &lRec) const override {
        if(!m_shape)
            throw NoriException("There is no shape attached to this Area light!");
        // https://www.pbr-book.org/3ed-2018/Light_Transport_I_Surface_Reflection/Sampling_Light_Sources#sec:sampling-lights
        // Lecture note direct illumination I, p22
        float cos_theta = lRec.n.dot(-lRec.wi);
        if(cos_theta <= 0) return 0.f;
        float pdf_A = m_shape->pdfSurface({lRec.ref, lRec.p});
        float ray_norm = (lRec.ref - lRec.p).squaredNorm();
        if(ray_norm == 0) return 0.f;
        float jacobian = abs(cos_theta) / ray_norm;
        return pdf_A / jacobian;
    }

    virtual Color3f samplePhoton(Ray3f &ray, const Point2f &sample1, const Point2f &sample2) const override {
        if(!m_shape)
            throw NoriException("There is no shape attached to this Area light!");
        ShapeQueryRecord sRec;
        m_shape -> sampleSurface(sRec, sample1);
        Vector3f direction = Warp::squareToCosineHemisphere(sample2);
        direction = Frame(sRec.n).toWorld(direction);
        ray = Ray3f(sRec.p, direction);
        EmitterQueryRecord lRec(sRec.p + direction, sRec.p, sRec.n);
        return eval(lRec) * M_PI / m_shape->pdfSurface(sRec);
    }

     bool onSurface() const override{return true;}

protected:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaEmitter, "area")
NORI_NAMESPACE_END