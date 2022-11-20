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
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the microfacet normal distribution D
    float evalBeckmann(const Normal3f &m) const {
        float temp = Frame::tanTheta(m) / m_alpha,
              ct = Frame::cosTheta(m), ct2 = ct*ct;

        return std::exp(-temp*temp) 
            / (M_PI * m_alpha * m_alpha * ct2 * ct2);
    }

    /// Evaluate Smith's shadowing-masking function G1 
    float smithBeckmannG1(const Vector3f &v, const Normal3f &m) const {
        float tanTheta = Frame::tanTheta(v);

        /* Perpendicular incidence -- no shadowing/masking */
        if (tanTheta == 0.0f)
            return 1.0f;

        /* Can't see the back side from the front and vice versa */
        if (m.dot(v) * Frame::cosTheta(v) <= 0)
            return 0.0f;

        float a = 1.0f / (m_alpha * tanTheta);
        if (a >= 1.6f)
            return 1.0f;
        float a2 = a * a;

        /* Use a fast and accurate (<0.35% rel. error) rational
           approximation to the shadowing-masking function */
        return (3.535f * a + 2.181f * a2) 
             / (1.0f + 2.276f * a + 2.577f * a2);
    }

    /// Evaluate the BRDF for the given pair of directions
    virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);
        Normal3f wh = (bRec.wi + bRec.wo).normalized();
        float D = evalBeckmann(wh);
        float F = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
    	float G = smithBeckmannG1(bRec.wi, wh)* smithBeckmannG1(bRec.wo, wh);
        float jacobian = 1.f / (4.f * Frame::cosTheta(bRec.wi)* Frame::cosTheta(bRec.wo));
        return m_kd * INV_PI + m_ks * D * F * G * jacobian;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    virtual float pdf(const BSDFQueryRecord &bRec) const override {
    	// throw NoriException("MicrofacetBRDF::pdf(): not implemented!");
        if (Frame::cosTheta(bRec.wo) <= 0) return 0.0f;
        Normal3f wh = (bRec.wi + bRec.wo).normalized();
        float cos_theta_o = Frame::cosTheta(bRec.wo);
        float cos_theta_h = Frame::cosTheta(wh);
        float cos_theta_h_i = bRec.wo.dot(wh);
        float D = evalBeckmann(wh);
        float jacobian = 1.f / (4.f * abs(cos_theta_h_i));
        return m_ks * D * cos_theta_h * jacobian + (1 - m_ks) * cos_theta_o * INV_PI;
    }

    /// Sample the BRDF
    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const override {
    	// throw NoriException("MicrofacetBRDF::sample(): not implemented!");
        float cos_theta_i = Frame::cosTheta(bRec.wi);
        if (cos_theta_i <= 0) return Color3f(0.0f);
        bRec.measure = ESolidAngle;
        // Warp from beckman or cosine semisphere according to ks
        if(_sample[0] < m_ks){
            Point2f sample(_sample[0] / m_ks, _sample[1]);
            Vector3f wh = Warp::squareToBeckmann(sample, m_alpha);
            float cos_theta_whi = bRec.wi.dot(wh);
            bRec.wo = (2.f * wh * cos_theta_whi - bRec.wi).normalized();
        } else {
            Point2f sample((_sample[0] - m_ks) / (1.f- m_ks), _sample[1]);
            bRec.wo = Warp::squareToCosineHemisphere(sample);
        }
        // Notice that we reverse wo and wi during sampling
        float cos_theta_o = Frame::cosTheta(bRec.wo);
    	if (cos_theta_o <= 0.f) return 0.f;
        return eval(bRec) * cos_theta_o / pdf(bRec);
    }

    virtual std::string toString() const override {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
