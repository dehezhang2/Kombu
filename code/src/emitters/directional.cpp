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

class DirectionalEmitter : public Emitter {
public:
    DirectionalEmitter(const PropertyList &props) {
        m_irradiance = props.getColor("radiance");
        m_localToWorld = props.getTransform("toWorld", Transform());
    }

    virtual std::string toString() const override {
        return tfm::format(
                "Directionlight[\n"
                "  radiance = %s,\n"
                "  transform = %s,\n"
                "]",
                m_irradiance.toString(),
                m_localToWorld.toString()
                );
    }

    virtual Color3f eval(const EmitterQueryRecord & lRec) const override {
        return Color3f(0.f);
    }

    virtual Color3f sample(EmitterQueryRecord & lRec, const Point2f & sample) const override {
        float radius = lRec.bSphere_radius * 1.1f;
        Vector3f d = m_localToWorld * Normal3f(0.f, 0.f, 1.f);
        Point3f disk_center = lRec.bSphere_center - d * radius;
        float distance = (lRec.ref - disk_center).dot(d);
        if(distance < 0) return 0.f;
        
        // light source
        lRec.p = lRec.ref - distance * d;
        // normal vector of the light source
        lRec.n = Normal3f(d);
        // direction from intersection to the light source
        lRec.wi = -d.normalized();
        lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, (lRec.p - lRec.ref).norm() - Epsilon);
        lRec.pdf = pdf(lRec);
        lRec.isDelta = true;
        return m_irradiance;
    }

    float pdf(const EmitterQueryRecord &lRec) const override {
        return 1.0f;
    }

    bool isDirectional() const override{
        return true;
    }

protected:
    Color3f m_irradiance;
    Transform m_localToWorld;
};

NORI_REGISTER_CLASS(DirectionalEmitter, "directional")
NORI_NAMESPACE_END