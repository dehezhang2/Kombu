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
#include <kombu/disney_utils.h>
#include <nori/common.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN


/// Ideal dielectric BSDF
class DisneyPrincipledGlass : public BSDF {
public:
    DisneyPrincipledGlass(const PropertyList &propList) {
        m_base_color = propList.getColor("base_color",0.5f);
        m_roughness = propList.getFloat("roughness", 0.5f);
        m_anisotropic = propList.getFloat("anisotropic", 0.0f);
        m_eta = propList.getFloat("eta", 1.f);
    }

    virtual Color3f eval(const BSDFQueryRecord & bRec) const override {
        float dotNWo = Frame::cosTheta(bRec.wo);
        float dotNWi = Frame::cosTheta(bRec.wi);

        bool reflect = dotNWo*dotNWi > 0;

        // Flip the shading frame if it is inconsistent with the geometry normal

        Color3f baseColor = m_base_color;
        float roughness = m_roughness;
        float eta = dotNWi > 0 ? m_eta : 1.f / m_eta;
        Vector3f wh = (bRec.wo+bRec.wi).normalized();
        if (!reflect) {
            // "Generalized half-vector" from Walter et al.
            // See "Microfacet Models for Refraction through Rough Surfaces"
            wh = (bRec.wi + bRec.wo * eta).normalized();
        }

        if(wh.dot(Vector3f(0.f,0.f,1.f))<0){
            wh = -wh;
        }

        float dotHWi = wh.dot(bRec.wi);
        float dotHWo = wh.dot(bRec.wo);

        float Fg = fresnel_dielectric(dotHWi,eta);

        float  D = GTR2(wh.z(),roughness);
        Frame frame(wh);
        float  G = smith_masking_gtr2(frame.toLocal(bRec.wi),roughness) *
                smith_masking_gtr2(frame.toLocal(bRec.wo),roughness);
        // Reflect case
        Color3f ret = 0.f;
        if(reflect){
            ret = baseColor * Fg * D * G  / (4 * abs(wh.dot(bRec.wi)));
        }
        // Refract Case
        else {
            // float eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
            float eta_factor = 1.f;
            float h_dot_out = wh.dot(bRec.wo);
            float h_dot_in = wh.dot(bRec.wi);
            float deom = h_dot_in+eta * h_dot_out;
            float F = Fg;
            Color3f res = baseColor * (eta_factor * (1 - F) * D * G  * fabs(h_dot_out * h_dot_in)) /
            (fabs(h_dot_in) * deom * deom );
            ret = sqrt(baseColor) * (1-Fg) * D * G * abs(h_dot_out * h_dot_in) /
                    ( fabs(h_dot_in) * abs(deom * deom) ) ;
        }
        return ret;
    }

    virtual float pdf(const BSDFQueryRecord & bRec) const override {
        float dotNWo = Frame::cosTheta(bRec.wo);
        float dotNWi = Frame::cosTheta(bRec.wi);

        bool reflect = dotNWo*dotNWi > 0;

        // Homework 1: implement this!
        float eta = dotNWi > 0 ? m_eta : 1.f / m_eta;
        // float eta = bRec.eta;
        assert(eta > 0);

        Vector3f wh = (bRec.wo+bRec.wi).normalized();
        if (!reflect) {
            // "Generalized half-vector" from Walter et al.
            // See "Microfacet Models for Refraction through Rough Surfaces"
            wh = (bRec.wi + bRec.wo * eta).normalized();
        }
        if(wh.dot(Vector3f(0.f,0.f,1.f))<0){
            wh = -wh;
        }

        float roughness = m_roughness;
        // Clamp roughness to avoid numerical issues.
        roughness = clamp(roughness, Epsilon, 1.f);

        // We sample the visible normals, also we use F to determine
        // whether to sample reflection or refraction
        // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
        float h_dot_in = wh.dot(bRec.wi);

        float F = fresnel_dielectric(h_dot_in, eta);
        float D = GTR2(wh.z(), roughness);
        float G_in = smith_masking_gtr2(bRec.wi, roughness);
        float prob = 0;
        if (reflect) {
            prob = (F * D * G_in) / (4 * fabs(Frame::cosTheta(bRec.wi)));
        } else {
            float h_dot_out = wh.dot(bRec.wo);
            float sqrt_denom = h_dot_in + eta * h_dot_out;
            float dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
            prob = (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / Frame::cosTheta(bRec.wi));
            // prob = std::min(1.0f,prob);
        }
        prob = std::max(Epsilon,prob);
        return prob;

    }

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const override {
    
        float dotNWi = Frame::cosTheta(bRec.wi);
        float eta = dotNWi > 0 ? m_eta : 1.f / m_eta;
            
        float roughness = m_roughness;
        // Clamp roughness to avoid numerical issues.
        roughness = clamp(roughness, Epsilon, 1.f);
        // Sample a micro normal and transform it to world space -- this is our half-vector.
        float alpha = roughness * roughness;
        Vector3f local_dir_in = bRec.wi;
        Vector3f local_micro_normal =
                sample_visible_normals(local_dir_in, alpha, sample);

        Vector3f wh = local_micro_normal;
        // Flip half-vector if it's below surface
        if (wh.z() < 0) {
            wh = -wh;
        }

        // Now we need to decide whether to reflect or refract.
        // We do this using the Fresnel term.
        float h_dot_in = wh.dot(bRec.wi);

        float F = fresnel_dielectric(h_dot_in, eta);

        float const rnd_param_w = m_sampler.nextFloat();
        if (rnd_param_w <= F) {
            // Reflection
            bRec.wo =  (-bRec.wi + 2 * bRec.wi.dot(wh) * wh).normalized();
            float cos_theta_o = Frame::cosTheta(bRec.wo);
            bRec.eta = 0.f;
            // bRec.measure = ESolidAngle;

            // set eta to 0 since we are not transmitting
            // return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness,false};
        } else {
            // Refraction
            // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
            // (note that our eta is eta2 / eta1, and l = -dir_in)
            float h_dot_out_sq = 1.f - (1.f - h_dot_in * h_dot_in) / (eta * eta);
            assert(h_dot_out_sq>0);
            // flip half_vector if needed
            if (h_dot_in < 0) {
                wh = -wh;
            }
            float h_dot_out= sqrt(h_dot_out_sq);
            bRec.wo = -bRec.wi / eta + (fabs(h_dot_in) / eta - h_dot_out) * wh;
            bRec.eta = eta;
            // return BSDFSampleRecord{refracted, eta, roughness,true};
        }
        Color3f ret = eval(bRec)* abs(Frame::cosTheta(bRec.wo))/pdf(bRec);

        return ret;
    }

    virtual std::string toString() const override {
        return tfm::format(
            "DisneyPrincipledGlass[\n"
            "  base_color = %f,%f,%f\n"
            "  roughness = %f\n"
            "  anisotropic = %f\n"
            "  eta = %f\n"
            "]",
            m_base_color[0],m_base_color[1],m_base_color[2],
            m_roughness,
            m_anisotropic,m_eta);
    }
private:
    Color3f m_base_color;
    float m_roughness;
    float m_anisotropic;
    float m_eta; 
    mutable pcg32 m_sampler;
};

NORI_REGISTER_CLASS(DisneyPrincipledGlass, "principle_glass");
NORI_NAMESPACE_END
