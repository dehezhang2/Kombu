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
#include <kombu/principledhelpers.h>
#include <kombu/microfacet.h>
#include <nori/warp.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN



/// Ideal dielectric BSDF
class DisneyPrincipled : public BSDF {
public:
    DisneyPrincipled(const PropertyList &propList) {
        m_base_color = propList.getColor("base_color",0.5f);
        m_roughness = propList.getFloat("roughness", 0.5f);
        m_has_anisotropic = propList.has("anisotropic");
        m_anisotropic = propList.getFloat("anisotropic", 0.0f);
        m_has_spec_trans = propList.has("spec_trans");
        m_spec_trans = propList.getFloat("spec_trans", 0.0f);
        m_has_sheen = propList.has("sheen");
        m_sheen = propList.getFloat("sheen", 0.0f);
        m_has_sheen_tint = propList.has("sheen_tint");
        m_sheen_tint = propList.getFloat("sheen_tint", 0.0f);
        m_has_flatness = propList.has("flatness");
        m_flatness = propList.getFloat("flatness", 0.0f);
        m_has_spec_tint = propList.has("spec_tint");
        m_spec_tint = propList.getFloat("spec_tint", 0.0f);
        m_has_metallic = propList.has("metallic");
        m_metallic = propList.getFloat("metallic", 0.0f);
        m_has_clearcoat = propList.has("clearcoat");
        m_clearcoat = propList.getFloat("clearcoat", 0.0f);
        m_clearcoat_gloss = propList.getFloat("clearcoat_gloss", 0.0f);
        m_spec_srate = propList.getFloat("main_specular_sampling_rate", 1.0f);
        m_clearcoat_srate = propList.getFloat("clearcoat_sampling_rate", 1.0f);
        m_diff_refl_srate = propList.getFloat("diffuse_reflectance_sampling_rate", 1.0f);

        m_specular = propList.getFloat("specular", 0.5f);
        m_eta_specular = false;
        if (m_has_spec_trans && m_specular==0.f){
            m_specular = 1e-3f;
        }
        m_eta = 2.0f / (1.0f - sqrt(0.08f * m_specular)) - 1.0f;

    }

    virtual Color3f eval(const BSDFQueryRecord & bRec) const override {
        bool active = true;
        float cos_theta_i = Frame::cosTheta(bRec.wi);
        // Ignore perfectly grazing configurations
        active &= (cos_theta_i!=0.0f);

        if (!active)
            return 0.0f;

        // Store the weights.
        float anisotropic = m_has_anisotropic ? m_anisotropic : 0.0f,
              roughness = m_roughness,
              flatness = m_has_flatness ? m_flatness: 0.0f,
              spec_trans = m_has_spec_trans ? m_spec_trans: 0.0f,
              metallic = m_has_metallic ? m_metallic : 0.0f,
              clearcoat = m_has_clearcoat ? m_clearcoat: 0.0f,
              sheen = m_has_sheen ? m_sheen : 0.0f;
        Color3f base_color = m_base_color;

        // Weights for BRDF and BSDF major lobes.
        float brdf = (1.0f - metallic) * (1.0f - spec_trans),
              bsdf = (1.0f - metallic) * spec_trans;

        float cos_theta_o = Frame::cosTheta(bRec.wo);

        // Reflection and refraction masks.
        bool reflect = cos_theta_i * cos_theta_o > 0.0f;
        bool refract = cos_theta_i * cos_theta_o < 0.0f;

        // Masks for the side of the incident ray (wi.z<0)
        bool front_side = cos_theta_i > 0.0f;
        float inv_eta   = 1.0f/m_eta;

        // Eta value w.r.t. ray instead of the object.
        float eta_path     = front_side? m_eta: inv_eta;
        float inv_eta_path = front_side? inv_eta: m_eta;

        // Main specular reflection and transmission lobe
        auto axay = calc_dist_params(anisotropic, roughness,m_has_anisotropic);
        float ax = axay.first,
              ay = axay.second;
        MicrofacetDistribution spec_dist(MicrofacetType::GGX, ax, ay);

        // Halfway vector
        Vector3f wh =
                (bRec.wi + bRec.wo * (reflect? 1.0f: eta_path)).normalized();

        // Make sure that the halfway vector points outwards the object
        if (Frame::cosTheta(wh)<0.f){
            wh = -wh;
        }

        // Dielectric Fresnel
        auto fresnel_ret = disney_fresnel(bRec.wi.dot(wh), m_eta);
        float F_spec_dielectric = std::get<0>(fresnel_ret);
        float cos_theta_t = std::get<1>(fresnel_ret);
        float eta_it = std::get<2>(fresnel_ret);
        float eta_ti = std::get<3>(fresnel_ret);


        bool reflection_compatibilty =
                mac_mic_compatibility(wh, bRec.wi, bRec.wo, cos_theta_i, true);
        bool refraction_compatibilty =
                mac_mic_compatibility(wh, bRec.wi, bRec.wo, cos_theta_i, false);
        // Masks for evaluating the lobes.
        // Specular reflection mask
        bool spec_reflect_active = active && reflect &&
                reflection_compatibilty &&
                (F_spec_dielectric > 0.0f);

        // Clearcoat mask
        bool clearcoat_active = m_has_clearcoat && active &&
                (clearcoat > 0.0f) && reflect &&
                reflection_compatibilty && front_side;

        // Specular transmission mask
        bool spec_trans_active = m_has_spec_trans && active && (bsdf > 0.0f) &&
                refract && refraction_compatibilty &&
                (F_spec_dielectric < 1.0f);

        // Diffuse, retro and fake subsurface mask
        bool diffuse_active = active && (brdf > 0.0f) && reflect && front_side;

        // Sheen mask
        bool sheen_active = m_has_sheen && active && (sheen > 0.0f) &&
                reflect && (1.0f - metallic > 0.0f) && front_side;

        // Evaluate the microfacet normal distribution
        float D = spec_dist.eval(wh);

        // Smith's shadowing-masking function
        float G = spec_dist.G(bRec.wi, bRec.wo, wh);

        // Initialize the final BSDF value.
        Color3f value(0.0f);

        // Main specular reflection evaluation
        if (spec_reflect_active) {
            // No need to calculate luminance if there is no color tint.
            float lum = m_has_spec_tint
                    ? luminance(base_color)
                    : 1.0f;
            float spec_tint =
                    m_has_spec_tint ? m_spec_tint: 0.0f;

            // Fresnel term
            Color3f F_principled = principled_fresnel(
                    F_spec_dielectric, metallic, spec_tint, base_color, lum,
                    bRec.wi.dot(wh), front_side, bsdf,m_eta,m_has_metallic,
                    m_has_spec_tint);

            // Adding the specular reflection component
            if(spec_reflect_active){
                value+= F_principled * D * G / (4.0f * abs(cos_theta_i));
            }
        }

        // Main specular transmission evaluation
        if (m_has_spec_trans && spec_trans_active) {

            /* Account for the solid angle compression when tracing
               radiance. This is necessary for bidirectional methods. */
            float scale = 1.0f;

            // Adding the specular transmission component
            if(spec_trans_active){
                value += sqrt(base_color) * bsdf *
                    abs((scale * (1.0f - F_spec_dielectric) * D * G * eta_path *
                    eta_path * bRec.wi.dot(wh) * bRec.wo.dot(wh)) /
                    (cos_theta_i * sqr(bRec.wi.dot(wh) +
                    eta_path * bRec.wo.dot(wh))));
            }
        }

        // Secondary isotropic specular reflection.
        if (m_has_clearcoat && clearcoat_active) {
            float clearcoat_gloss = m_clearcoat_gloss;

            // Clearcoat lobe uses the schlick approximation for Fresnel
            // term.
            float Fcc = calc_schlick(0.04f, bRec.wi.dot(wh),m_eta);

            /* Clearcoat lobe uses GTR1 distribution. Roughness is mapped
             * between 0.1 and 0.001. */
            GTR1Isotropic mfacet_dist(disney_lerp(0.1f, 0.001f, clearcoat_gloss));
            float Dcc = mfacet_dist.eval(wh);

            // Shadowing shadowing-masking term
            float G_cc = clearcoat_G(bRec.wi, bRec.wo, wh, 0.25f);

            // Adding the clearcoat component.
            if(clearcoat_active){
                value+=(clearcoat * 0.25f) * Fcc * Dcc * G_cc * abs(cos_theta_o);
            }
        }

        // Evaluation of diffuse, retro reflection, fake subsurface and
        // sheen.
        if (diffuse_active) {
            float Fo = schlick_weight(abs(cos_theta_o)),
            Fi = schlick_weight(abs(cos_theta_i));

            // Diffuse
            float f_diff = (1.0f - 0.5f * Fi) * (1.0f - 0.5f * Fo);

            float cos_theta_d = wh.dot(bRec.wo);
            float Rr          = 2.0f * roughness * sqr(cos_theta_d);

            // Retro reflection
            float f_retro = Rr * (Fo + Fi + Fo * Fi * (Rr - 1.0f));

            if (m_has_flatness) {
                /* Fake subsurface implementation based on Hanrahan Krueger
                   Fss90 used to "flatten" retro reflection based on
                   roughness.*/
                float Fss90 = Rr / 2.0f;
                float Fss =
                        disney_lerp(1.0f, Fss90, Fo) * disney_lerp(1.0f, Fss90, Fi);

                float f_ss = 1.25f * (Fss * (1.0f / (abs(cos_theta_o) +
                        abs(cos_theta_i)) -
                                0.5f) +
                                        0.5f);

                // Adding diffuse, retro and fake subsurface evaluation.
                if(diffuse_active){
                    value += brdf * abs(cos_theta_o) * base_color *
                             M_PI *
                             (disney_lerp(f_diff + f_retro, f_ss, flatness));
                }

            } else {
                // Adding diffuse, retro evaluation. (no fake ss.)
                if(diffuse_active){
                    value+=brdf * abs(cos_theta_o) * base_color *
                        M_PI * (f_diff + f_retro);
                }
            }
            // Sheen evaluation
            if (m_has_sheen && sheen_active) {
                float Fd = schlick_weight(abs(cos_theta_d));

                // Tint the sheen evaluation towards the base color.
                if (m_has_sheen_tint) {
                    float sheen_tint = m_sheen_tint;

                    // Luminance evaluation
                    float lum = luminance(base_color);

                    // Normalize color with luminance and tint the result.
                    Color3f c_tint = base_color / lum;
                    if(lum<=0.0f){
                        c_tint = 1.0f;
                    }
                    Color3f c_sheen = disney_lerp(1.0f, c_tint, sheen_tint);

                    // Adding sheen evaluation with tint.
                    if(sheen_active){
                        value += sheen * (1.0f - metallic) * Fd * c_sheen *
                                 abs(cos_theta_o);
                    }
                } else {
                    // Adding sheen evaluation without tint.
                    if(sheen_active){
                        value += sheen * (1.0f - metallic) * Fd * abs(cos_theta_o);
                    }
                }
            }
        }
        if (!active){
            return 0.f;
        }
        return value;
    }

    virtual float pdf(const BSDFQueryRecord &bRec) const override {

        float cos_theta_i = Frame::cosTheta(bRec.wi);
        // Ignore perfectly grazing configurations.
        bool active = true;
        active &= cos_theta_i != 0.0f;

        if (!active)
            return 0.0f;

        // Store the weights.
        float anisotropic =
                m_has_anisotropic ? m_anisotropic: 0.0f,
                roughness = m_roughness,
                spec_trans =
                        m_has_spec_trans ? m_spec_trans: 0.0f;
        float metallic = m_has_metallic ? m_metallic: 0.0f,
        clearcoat =
                m_has_clearcoat ? m_clearcoat: 0.0f;

        // BRDF and BSDF major lobe weights
        float brdf = (1.0f - metallic) * (1.0f - spec_trans),
        bsdf = (1.0f - metallic) * spec_trans;

        // Masks if incident direction is inside (wi.z<0)
        bool front_side = cos_theta_i > 0.0f;

        // Eta w.r.t. light path.
        float eta_path    = (front_side? m_eta: 1.f/(m_eta));
        float cos_theta_o = Frame::cosTheta(bRec.wo);

        bool reflect = cos_theta_i * cos_theta_o > 0.0f;
        bool refract = cos_theta_i * cos_theta_o < 0.0f;

        // Halfway vector calculation
        Vector3f wh = (bRec.wi + bRec.wo * (reflect? 1.0f: eta_path)).normalized();

        // Make sure that the halfway vector points outwards the object
        wh = mulsign(wh, Frame::cosTheta(wh));

        // Main specular distribution for reflection and transmission.
        auto axay = calc_dist_params(anisotropic, roughness,m_has_anisotropic);
        float ax = axay.first;
        float ay = axay.second;
        MicrofacetDistribution spec_distr(MicrofacetType::GGX, ax, ay);

        // Dielectric Fresnel calculation
        auto disney_fresnel_out =  disney_fresnel(bRec.wi.dot(wh), m_eta);
        float F_spec_dielectric = std::get<0>(disney_fresnel_out);
        float cos_theta_t = std::get<1>(disney_fresnel_out);
        float eta_it = std::get<2>(disney_fresnel_out);
        float eta_ti = std::get<3>(disney_fresnel_out);

        // Defining the probabilities
        float prob_spec_reflect = (
                front_side?
                m_spec_srate * (1.0f - bsdf * (1.0f - F_spec_dielectric)):
                F_spec_dielectric);
        float prob_spec_trans =
                m_has_spec_trans
                ? (front_side?
                    (m_spec_srate * bsdf * (1.0f - F_spec_dielectric)):
                    (1.0f - F_spec_dielectric))
                : 0.0f;
        float prob_clearcoat =
                m_has_clearcoat
                ? (front_side? 0.25f * clearcoat * m_clearcoat_srate:0.0f)
                : 0.0f;
        float prob_diffuse =
                (front_side? brdf * m_diff_refl_srate: 0.f);

        // Normalizing the probabilities.
        float rcp_tot_prob = 1.f/(prob_spec_reflect + prob_spec_trans +
                prob_clearcoat + prob_diffuse);
        prob_spec_reflect *= rcp_tot_prob;
        prob_spec_trans *= rcp_tot_prob;
        prob_clearcoat *= rcp_tot_prob;
        prob_diffuse *= rcp_tot_prob;

        /* Calculation of dwh/dwo term. Different for reflection and
         transmission. */
        float dwh_dwo_abs;
        if (m_has_spec_trans) {
            float dot_wi_h = bRec.wi.dot(wh);
            float dot_wo_h = bRec.wo.dot(wh);
            dwh_dwo_abs    = abs((reflect? 
                                  1.f/(4.0f * dot_wo_h):
                                  (sqr(eta_path) * dot_wo_h) /
                                   sqr(dot_wi_h + eta_path * dot_wo_h)));
        } else {
            dwh_dwo_abs = abs(1.f/(4.0f * bRec.wo.dot(wh)));
        }

        // Initializing the final pdf value.
        float pdf(0.0f);

        // Macro-micro surface compatibility mask for reflection.
        bool mfacet_reflect_macmic =
                mac_mic_compatibility(wh, bRec.wi, bRec.wo, cos_theta_i, true) && reflect;

        // Adding main specular reflection pdf
        if(mfacet_reflect_macmic){
            pdf+=prob_spec_reflect *
                spec_distr.pdf(mulsign(bRec.wi, cos_theta_i), wh) * dwh_dwo_abs;
        }
        // Adding cosine hemisphere reflection pdf
        if(reflect){
            pdf+=prob_diffuse * Warp::squareToCosineHemispherePdf(bRec.wo);
        }

        // Main specular transmission
        if (m_has_spec_trans) {
            // Macro-micro surface mask for transmission.
            bool mfacet_trans_macmic =
                    mac_mic_compatibility(wh, bRec.wi, bRec.wo, cos_theta_i, false) &&
                    refract;

            // Adding main specular transmission pdf
            if(mfacet_trans_macmic){
                pdf += prob_spec_trans *
                    spec_distr.pdf(mulsign(bRec.wi, cos_theta_i), wh) *
                    dwh_dwo_abs;
            }
        }
        // Adding the secondary specular reflection pdf.(clearcoat)
        if (m_has_clearcoat) {
            float clearcoat_gloss = m_clearcoat_gloss;
            GTR1Isotropic cc_dist(disney_lerp(0.1f, 0.001f, clearcoat_gloss));
            if(mfacet_reflect_macmic){
                pdf+=prob_clearcoat * cc_dist.pdf(wh) * dwh_dwo_abs;
            }
        }
        return pdf;
    }

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const override {
        float sample1 = m_sampler.nextFloat();

        float cos_theta_i = Frame::cosTheta(bRec.wi);
        bool active = true;
        // Ignoring perfectly grazing incoming rays
        active &= (cos_theta_i!=0.0f);

        if (!active)
            return 0.0f;

        // Store the weights.
        float anisotropic = m_has_anisotropic ? m_anisotropic: 0.0f,
        roughness = m_roughness,
        spec_trans = m_has_spec_trans ? m_spec_trans : 0.0f,
        metallic = m_has_metallic ? m_metallic : 0.0f,
        clearcoat = m_has_clearcoat ? m_clearcoat : 0.0f;

        // Weights of BSDF and BRDF major lobes
        float brdf = (1.0f - metallic) * (1.0f - spec_trans),
        bsdf = m_has_spec_trans ? (1.0f - metallic) * spec_trans : 0.0f;

        // Mask for incident side. (wi.z<0)
        bool front_side = cos_theta_i > 0.0f;

        // Defining main specular reflection distribution
        auto axay = calc_dist_params(anisotropic, roughness,m_has_anisotropic);
        float ax = axay.first;
        float ay = axay.second;

        MicrofacetDistribution spec_distr(MicrofacetType::GGX, ax, ay);
        Normal3f m_spec = std::get<0>(
                spec_distr.sample(mulsign(bRec.wi, cos_theta_i), sample));

        // Fresnel coefficient for the main specular.
        auto disney_fresnel_out = disney_fresnel(bRec.wi.dot(m_spec), m_eta);
        float F_spec_dielectric = std::get<0>(disney_fresnel_out);
        float cos_theta_t = std::get<1>(disney_fresnel_out);
        float eta_it = std::get<2>(disney_fresnel_out);
        float eta_ti = std::get<3>(disney_fresnel_out);

        // If BSDF major lobe is turned off, we do not sample the inside
        // case.
        active &= (front_side || (bsdf > 0.0f));

        // Probability definitions
        /* Inside  the material, just microfacet Reflection and
           microfacet Transmission is sampled. */
        float prob_spec_reflect = (
                front_side?
                m_spec_srate * (1.0f - bsdf * (1.0f - F_spec_dielectric)):
                F_spec_dielectric);
        float prob_spec_trans =
                m_has_spec_trans
                ? (front_side?
                    m_spec_srate * bsdf * (1.0f - F_spec_dielectric):
                    (1.0f - F_spec_dielectric))
                : 0.0f;
        // Clearcoat has 1/4 of the main specular reflection energy.
        float prob_clearcoat =
                m_has_clearcoat
                ? (front_side? 0.25f * clearcoat * m_clearcoat_srate: 0.0f)
                : 0.0f;
        float prob_diffuse = (front_side? brdf * m_diff_refl_srate: 0.0f);

        // Normalizing the probabilities.
        float rcp_tot_prob = 1.f/(prob_spec_reflect + prob_spec_trans +
                prob_clearcoat + prob_diffuse);
        prob_spec_trans *= rcp_tot_prob;
        prob_clearcoat *= rcp_tot_prob;
        prob_diffuse *= rcp_tot_prob;

        // Sampling mask definitions
        float curr_prob(0.0f);
        bool sample_diffuse = active && (sample1 < prob_diffuse);
        curr_prob += prob_diffuse;
        bool sample_clearcoat = m_has_clearcoat && active &&
                (sample1 >= curr_prob) &&
                (sample1 < curr_prob + prob_clearcoat);
        curr_prob += prob_clearcoat;
        bool sample_spec_trans = m_has_spec_trans && active &&
                (sample1 >= curr_prob) &&
                (sample1 < curr_prob + prob_spec_trans);
        curr_prob += prob_spec_trans;
        bool sample_spec_reflect = active && (sample1 >= curr_prob);

        // Eta will be changed in transmission.
        bRec.eta = 1.0f;

        // Main specular reflection sampling
        if (sample_spec_reflect) {
            Vector3f wo = reflect(bRec.wi, m_spec);
            if(sample_spec_reflect){
                bRec.wo = wo;
            }
            // dr::masked(bs.sampled_component, sample_spec_reflect) = 3;
            // dr::masked(bs.sampled_type, sample_spec_reflect) =
            //         +BSDFFlags::GlossyReflection;

            /* Filter the cases where macro and micro surfaces do not agree
             on the same side and reflection is not successful*/
            bool reflect = cos_theta_i * Frame::cosTheta(wo) > 0.0f;
            active &=
                    (!sample_spec_reflect ||
                    (mac_mic_compatibility(Vector3f(m_spec),
                                           bRec.wi, wo, cos_theta_i, true) &&
                    reflect));
        }
        // The main specular transmission sampling
        if (m_has_spec_trans && sample_spec_trans) {
            Vector3f wo = refract(bRec.wi, m_spec, cos_theta_t, eta_ti);
            if(sample_spec_trans){
                bRec.wo = wo;
                bRec.eta = eta_it;
            }
            // dr::masked(bs.sampled_component, sample_spec_trans) = 2;
            // dr::masked(bs.sampled_type, sample_spec_trans) =
            //         +BSDFFlags::GlossyTransmission;
            // dr::masked(bs.eta, sample_spec_trans) = eta_it;

            /* Filter the cases where macro and micro surfaces do not agree
             on the same side and refraction is successful. */
            bool refract = cos_theta_i * Frame::cosTheta(wo) < 0.0f;
            active &= (!sample_spec_trans ||
                    (mac_mic_compatibility(Vector3f(m_spec),
                                           bRec.wi, wo, cos_theta_i,
                                           false) &&
                                           refract));
        }
        // The secondary specular reflection sampling (clearcoat)
        if (m_has_clearcoat && sample_clearcoat) {
            float clearcoat_gloss = m_clearcoat_gloss;

            // Clearcoat roughness is mapped between 0.1 and 0.001.
            GTR1Isotropic cc_dist(disney_lerp(0.1f, 0.001f, clearcoat_gloss));
            Normal3f m_cc                = cc_dist.sample(sample);
            Vector3f wo                         = reflect(bRec.wi, m_cc);
            if(sample_clearcoat){
                bRec.wo = wo;
            }
            // dr::masked(bs.sampled_component, sample_clearcoat) = 1;
            // dr::masked(bs.sampled_type, sample_clearcoat) =
                    // +BSDFFlags::GlossyReflection;

            /* Filter the cases where macro and microfacets do not agree on
             the same side and reflection is not successful. */
            bool reflect = cos_theta_i * Frame::cosTheta(wo) > 0.0f;
            active &= (!sample_clearcoat ||
                    (mac_mic_compatibility(Vector3f(m_cc),
                                           bRec.wi, wo,
                                           cos_theta_i, true) &&
                                           reflect));
        }
        // Cosine hemisphere reflection sampling
        if (sample_diffuse) {
            Vector3f wo = Warp::squareToCosineHemisphere(sample);
            if(sample_diffuse){
                bRec.wo = wo;
            }
            // dr::masked(bs.sampled_component, sample_diffuse) = 0;
            // dr::masked(bs.sampled_type, sample_diffuse) =
            //         +BSDFFlags::DiffuseReflection;
            bool reflect = cos_theta_i * Frame::cosTheta(wo) > 0.0f;
            active &= (!sample_diffuse || reflect);
        }

        float res_pdf = pdf(bRec);
        active &= res_pdf > 0.0f;
        Color3f result = eval(bRec);
        if (!active){
            return 0.f;
        }
        return result / res_pdf;
    }

    virtual std::string toString() const override {
        return tfm::format(
            "DisneyPrincipled[\n"
            "  base_color = %f,%f,%f\n"
            "  specular_transmission = %f\n"
            "  metallic = %f\n"
            "  flatness = %f\n"
            "  specular = %f\n"
            "  roughness = %f\n"
            "  specular_tint = %f\n"
            "  anisotropic = %f\n"
            "  sheen = %f\n"
            "  sheen_tint = %f\n"
            "  clearcoat = %f\n"
            "  clearcoat_gloss = %f\n"
            "]",
            m_base_color[0],m_base_color[1],m_base_color[2],
            m_spec_trans,m_metallic, m_flatness,
            m_specular,m_roughness,m_spec_tint,
            m_anisotropic,m_sheen,m_sheen_tint,
            m_clearcoat,m_clearcoat_gloss);
    }
private:
    Color3f m_base_color;
    float m_roughness;
    float m_anisotropic;
    float m_sheen;
    float m_sheen_tint;
    float m_spec_trans;
    float m_flatness;
    float m_spec_tint;
    float m_clearcoat;
    float m_clearcoat_gloss;
    float m_metallic;
    float m_specular;
    float m_eta;
    float m_diff_refl_srate;
    float m_spec_srate;
    float m_clearcoat_srate;

    bool m_has_clearcoat;
    bool m_has_sheen;
    bool m_has_spec_trans;
    bool m_has_metallic;
    bool m_has_spec_tint;
    bool m_has_sheen_tint;
    bool m_has_anisotropic;
    bool m_has_flatness;

    bool m_eta_specular;

    mutable pcg32 m_sampler;
};

NORI_REGISTER_CLASS(DisneyPrincipled, "principled");
NORI_NAMESPACE_END
