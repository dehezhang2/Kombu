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

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    virtual Color3f eval(const BSDFQueryRecord &) const override {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    virtual float pdf(const BSDFQueryRecord &) const override {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const override {
        float cos_theta_i = Frame::cosTheta(bRec.wi);
        
        bRec.measure = EDiscrete;
        float F = fresnel(cos_theta_i, m_extIOR, m_intIOR);
        if(sample[0] < F){
            bRec.eta = 1.0f;
            bRec.wo = Vector3f(
                -bRec.wi.x(),
                -bRec.wi.y(),
                bRec.wi.z()
            );
        } else {
            bRec.eta = m_intIOR / m_extIOR;
            float eta1_over_eta2 = 1.0f / bRec.eta;
            Normal3f n(0.0f, 0.0f, 1.0f);
            // check the direction of wi
            if(cos_theta_i < 0){
                eta1_over_eta2 = bRec.eta;
                n[2] = -n[2];
                cos_theta_i = -cos_theta_i;
            }
            bRec.wo = - eta1_over_eta2 * (bRec.wi - cos_theta_i * n)
                      - sqrt(1.0f - eta1_over_eta2 * eta1_over_eta2 * (1 - cos_theta_i * cos_theta_i)) * n;
            bRec.wo = bRec.wo.normalized();
            return eta1_over_eta2 * eta1_over_eta2;
        }
        return 1.0f;
    }

    virtual std::string toString() const override {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
