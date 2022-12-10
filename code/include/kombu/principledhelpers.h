#if !defined(__NORI_PRINCIPLEDHELPERS_H)
#define __NORI_PRINCIPLEDHELPERS_H

#include <nori/object.h>
#include <nori/frame.h>


NORI_NAMESPACE_BEGIN

inline float rsqrt(float a){
    return 1.0f/(sqrt(a));
}
inline float fmadd(float a, float b, float c){
    return a*b+c;
} 
inline float sqr(float a){
    return a*a;
}
inline float safe_sqrt(float a){
    return sqrt(std::max(a,0.f));
}
inline float disney_lerp(float a, float b, float t){
    return fmadd(b, t, fmadd(-a, t, a));
}
inline Color3f disney_lerp(float a, float b, Color3f t){
    Color3f c = -a*t + a;
    return b*t +c;
}
inline Color3f disney_lerp(float a, Color3f b, float t){
    Color3f c = -a*t + a;
    return b*t +c;
}

inline float fnmadd(float a, float b, float c){
    return -a*b+c;
}
template <typename T>
inline T mulsign(T v1, float v2){
    if(v2>=0){
        return v1;
    }
    else{
        return -v1;
    }
}

template <typename T>
inline T mulsign_neg(T v1, float v2){
    if(v2>=0){
        return -v1;
    }
    else{
        return v1;
    }
}

Vector3f reflect(const Vector3f &wi, const Normal3f &m) {
    return Vector3f(m)* 2.f * wi.dot(m) - wi;
}

Vector3f refract(const Vector3f &wi, const Normal3f &m, float cos_theta_t,
                         float eta_ti) {
    float temp = fmadd(wi.dot(m),eta_ti, cos_theta_t);
    return m*temp - wi * eta_ti;
}

std::pair<float, float> calc_dist_params(float anisotropic,
                                         float roughness,
                                         bool has_anisotropic){
    float roughness_2 = roughness*roughness;
    if (!has_anisotropic) {
        float a = std::max(0.001f, roughness_2);
        return { a, a };
    }
    float aspect = sqrt(1.0f - 0.9f * anisotropic);
    return { std::max(0.001f, roughness_2 / aspect),
             std::max(0.001f, roughness_2 * aspect) };
}


std::tuple<float, float, float, float> disney_fresnel(float cos_theta_i, float eta) {
    auto outside_mask = cos_theta_i >= 0.f;

    float rcp_eta = 1.f/(eta),
          eta_it = (outside_mask? eta: rcp_eta),
          eta_ti = (outside_mask? rcp_eta: eta);

    /* Using Snell's law, calculate the squared sine of the
       angle between the surface normal and the transmitted ray */
    float cos_theta_t_sqr =
        fnmadd(fnmadd(cos_theta_i, cos_theta_i, 1.f), eta_ti * eta_ti, 1.f);

    /* Find the absolute cosines of the incident/transmitted rays */
    float cos_theta_i_abs = abs(cos_theta_i);
    float cos_theta_t_abs = safe_sqrt(cos_theta_t_sqr);

    auto index_matched = (eta==1.f),
         special_case  = index_matched || (cos_theta_i_abs== 0.f);

    float r_sc = (index_matched? float(0.f): float(1.f));

    /* Amplitudes of reflected waves */
    float a_s = fnmadd(eta_it, cos_theta_t_abs, cos_theta_i_abs) /
                fmadd(eta_it, cos_theta_t_abs, cos_theta_i_abs);

    float a_p = fnmadd(eta_it, cos_theta_i_abs, cos_theta_t_abs) /
                fmadd(eta_it, cos_theta_i_abs, cos_theta_t_abs);

    float r = 0.5f * (sqr(a_s) + sqr(a_p));

    if(special_case){
        r = r_sc;
    }

    /* Adjust the sign of the transmitted direction */
    float cos_theta_t = cos_theta_t_abs;
    if(cos_theta_i<0){
        cos_theta_t = -cos_theta_t_abs;
    }

    return { r, cos_theta_t, eta_it, eta_ti };
}


float schlick_weight(float cos_i) {
    float m = clamp(1.0f - cos_i, 0.0f, 1.0f);
    return sqr(sqr(m)) * m;
}

Color3f calc_schlick(Color3f R0, float cos_theta_i,float eta){
    bool outside_mask = cos_theta_i >= 0.0f;
    float rcp_eta     = 1.0f/(eta),
    eta_it      = (outside_mask? eta: rcp_eta),
    eta_ti      = (outside_mask? rcp_eta: eta);

    float cos_theta_t_sqr = fnmadd(
            fnmadd(cos_theta_i, cos_theta_i, 1.0f), sqr(eta_ti), 1.0f);
    float cos_theta_t = safe_sqrt(cos_theta_t_sqr);
    return (eta_it > 1.0f?
            disney_lerp(schlick_weight(abs(cos_theta_i)), 1.0f, R0):
            disney_lerp(schlick_weight(cos_theta_t), 1.0f, R0));
}

float calc_schlick(float R0, float cos_theta_i,float eta){
    bool outside_mask = cos_theta_i >= 0.0f;
    float rcp_eta     = 1.0f/(eta),
    eta_it      = (outside_mask? eta: rcp_eta),
    eta_ti      = (outside_mask? rcp_eta: eta);

    float cos_theta_t_sqr = fnmadd(
            fnmadd(cos_theta_i, cos_theta_i, 1.0f), sqr(eta_ti), 1.0f);
    float cos_theta_t = safe_sqrt(cos_theta_t_sqr);
    return (eta_it > 1.0f?
            disney_lerp(schlick_weight(abs(cos_theta_i)), 1.0f, R0):
            disney_lerp(schlick_weight(cos_theta_t), 1.0f, R0));
}

float schlick_R0_eta(float eta){
    return sqr((eta - 1.0f) / (eta + 1.0f));
}

Color3f principled_fresnel(const float &F_dielectric, const float &metallic,
                     const float &spec_tint,
                     const Color3f &base_color,
                     const float &lum, const float &cos_theta_i,
                     const bool front_side,
                     const float &bsdf, const float &eta,
                     bool has_metallic, bool has_spec_tint) {
    // Outside mask based on micro surface
    bool outside_mask = cos_theta_i >= 0.0f;
    float rcp_eta = 1.f/(eta);
    float eta_it  = (outside_mask? eta: rcp_eta);
    Color3f F_schlick = (0.0f);

    // Metallic component based on Schlick.
    if (has_metallic) {
        F_schlick += metallic * calc_schlick(
                base_color, cos_theta_i,eta);
    }

    // Tinted dielectric component based on Schlick.
    if (has_spec_tint) {
        Color3f c_tint = base_color / lum;
        if (lum<=0.0f){
            c_tint = 1.0f;
        }
        Color3f F0_spec_tint =
                c_tint * schlick_R0_eta(eta_it);
        F_schlick +=
                (1.0f - metallic) * spec_tint *
                calc_schlick(F0_spec_tint, cos_theta_i,eta);
    }

    // Front side fresnel.
    Color3f F_front =
            (1.0f - metallic) * (1.0f - spec_tint) * F_dielectric + F_schlick;
    /* For back side there is no tint or metallic, just true dielectric
       fresnel.*/
    return (front_side? F_front: bsdf * F_dielectric);
}


bool mac_mic_compatibility(const Vector3f &m,
                            const Vector3f &wi,
                            const Vector3f &wo,
                            const float &cos_theta_i,
                            bool reflection) {
    if (reflection) {
        return (wi.dot(mulsign(m, cos_theta_i)) > 0.0f) &&
        (wo.dot(mulsign(m, cos_theta_i)) > 0.0f);
    } else {
        return (wi.dot(mulsign(m, cos_theta_i)) > 0.0f) &&
        (wo.dot(mulsign_neg(m, cos_theta_i)) > 0.0f);
    }
}

float luminance(const Color3f &c) {
    return c[0] * 0.212671f + c[1] * 0.715160f + c[2] * 0.072169f;
}


class GTR1Isotropic {
public:
    GTR1Isotropic(float alpha) : m_alpha(alpha){};

    float eval(const Vector3f &m) const {
        float cos_theta  = Frame::cosTheta(m),
        cos_theta2 = sqr(cos_theta), alpha2 = sqr(m_alpha);

        float result = (alpha2 - 1.f) / (M_PI * log(alpha2) *
                (1.f + (alpha2 - 1.f) * cos_theta2));

        return (result * cos_theta > 1e-20f? result: 0.f);
    }

    float pdf(const Vector3f &m) const {
        return (m.z() < 0.f? 0.f: Frame::cosTheta(m) * eval(m));
    }

    Normal3f sample(const Point2f &sample) const {
        float sin_phi = sin((2.f * M_PI) *sample.x());
        float cos_phi = cos((2.f * M_PI) *sample.x());
        float alpha2            = sqr(m_alpha);

        float cos_theta2 =
                (1.f - pow(alpha2, 1.f - sample.y())) / (1.f - alpha2);

        float sin_theta = sqrt(std::max(0.f, 1.f - cos_theta2)),
        cos_theta = sqrt(std::max(0.f, cos_theta2));
        return Normal3f(cos_phi * sin_theta, sin_phi * sin_theta, cos_theta);
    }

private:
    float m_alpha;
};


float smith_ggx1(const Vector3f &v, const Vector3f &wh,
                 const float &alpha) {
    float alpha_2     = sqr(alpha),
    cos_theta   = abs(Frame::cosTheta(v)),
    cos_theta_2 = sqr(cos_theta),
    tan_theta_2 = (1.0f - cos_theta_2) / cos_theta_2;

    float result =
            2.0f * 1.0f/(1.0f + sqrt(1.0f + alpha_2 * tan_theta_2));

    // Perpendicular incidence -- no shadowing/masking
    if (v.z()==1.f){
        result = 1.f;
    }
    /* Ensure consistent orientation (can't see the back
       of the microfacet from the front and vice versa) */
    if(v.dot(wh) * Frame::cosTheta(v) <= 0.f){
        result = 0.f;
    }
    return result;
}

float clearcoat_G(const Vector3f &wi, const Vector3f &wo,
                  const Vector3f&wh, const float &alpha) {
    return smith_ggx1(wi, wh, alpha) * smith_ggx1(wo, wh, alpha);
}



NORI_NAMESPACE_END
#endif