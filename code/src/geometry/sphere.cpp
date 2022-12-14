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

#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Sphere : public Shape {
public:
    Sphere(const PropertyList & propList) {
        m_position = propList.getPoint3("center", Point3f());
        m_radius = propList.getFloat("radius", 1.f);

        m_bbox.expandBy(m_position - Vector3f(m_radius));
        m_bbox.expandBy(m_position + Vector3f(m_radius));
    }

    virtual BoundingBox3f getBoundingBox(uint32_t index) const override { return m_bbox; }

    virtual Point3f getCentroid(uint32_t index) const override { return m_position; }

    virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const override {
	    /* to be implemented */
        float a = ray.d.squaredNorm();
        float b = 2 * (ray.o - m_position).dot(ray.d);
        float c = (ray.o - m_position).squaredNorm() - m_radius * m_radius;
        float delta = b * b - 4 * a * c;
        if(delta < 0) return false;
        float t0 = (-b - sqrt(delta)) /(2 * a);
        float t1 = (-b + sqrt(delta)) /(2 * a);
        t = ((t0 >= ray.mint && t0 <= ray.maxt) ? t0 : t1);
        return (t >= ray.mint && t <= ray.maxt);
    }

    virtual void setHitInformation(uint32_t index, const Ray3f &ray, Intersection & its) const override {
        /* to be implemented */
        its.p = ray.o + its.t * ray.d;
        Vector3f normalized_coor = (its.p - m_position).normalized();
        its.uv[0] = atan2(normalized_coor[1], normalized_coor[0]) / (2 * M_PI) + 0.5;
        its.uv[1] = 1-acos(normalized_coor[2]) / M_PI;
        its.geoFrame = Frame(normalized_coor);
		its.shFrame = Frame(normalized_coor);

        Vector3f local = its.p - m_position;
        // its.dpdu = m_objectToWorld(Vector(-local.y, local.x, 0) * (2*M_PI));
        its.dpdu = Vector3f(-local.y(), local.x(), 0) * (2*M_PI);
        float zrad = std::sqrt(local.x()*local.x() + local.y()*local.y());
        float theta = safe_acos(local.z()/m_radius);
        if (zrad > 0) {
            float invZRad = 1.0f / zrad,
                  cosPhi = local.x() * invZRad,
                  sinPhi = local.y() * invZRad;
            // its.dpdv = m_objectToWorld(Vector(local.z * cosPhi, local.z * sinPhi,
            //         -std::sin(theta)*m_radius) * M_PI);
            its.dpdv = Vector3f(local.z() * cosPhi, local.z() * sinPhi,
                    -std::sin(theta)*m_radius) * M_PI;
        } else {
            // avoid a singularity
            const float cosPhi = 0, sinPhi = 1;
            // its.dpdv = m_objectToWorld(Vector(local.z * cosPhi, local.z * sinPhi,
            //         -std::sin(theta)*m_radius) * M_PI);
            its.dpdv = Vector3f(local.z() * cosPhi, local.z() * sinPhi,
                    -std::sin(theta)*m_radius) * M_PI;
        }
    }

    virtual void sampleSurface(ShapeQueryRecord & sRec, const Point2f & sample) const override {
        Vector3f q = Warp::squareToUniformSphere(sample);
        sRec.p = m_position + m_radius * q;
        sRec.n = q;
        sRec.pdf = std::pow(1.f/m_radius,2) * Warp::squareToUniformSpherePdf(Vector3f(0.0f,0.0f,1.0f));
    }
    virtual float pdfSurface(const ShapeQueryRecord & sRec) const override {
        return std::pow(1.f/m_radius,2) * Warp::squareToUniformSpherePdf(Vector3f(0.0f,0.0f,1.0f));
    }


    virtual std::string toString() const override {
        return tfm::format(
                "Sphere[\n"
                "  center = %s,\n"
                "  radius = %f,\n"
                "  bsdf = %s,\n"
                "  emitter = %s\n"
                "]",
                m_position.toString(),
                m_radius,
                m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
                m_emitter ? indent(m_emitter->toString()) : std::string("null"));
    }

protected:
    Point3f m_position;
    float m_radius;
};

NORI_REGISTER_CLASS(Sphere, "sphere");
NORI_NAMESPACE_END
