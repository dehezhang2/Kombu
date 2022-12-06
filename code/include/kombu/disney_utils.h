#if !defined(__NORI_DISNEYTUILS_H)
#define __NORI_DISNEYTUILS_H

#include <nori/object.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN


template <typename T>
inline T schlick_fresnel(const T &F0, float cos_theta) {
    return F0 + (1.f - F0) *
        pow(std::max(1.f - cos_theta, 0.f), 5.f);
}

/// Fresnel equation of a dielectric interface.
/// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
/// n_dot_i: abs(cos(incident angle))
/// n_dot_t: abs(cos(transmission angle))
/// eta: eta_transmission / eta_incident
inline float fresnel_dielectric(float n_dot_i, float n_dot_t, float eta) {
    assert(n_dot_i >= 0 && n_dot_t >= 0 && eta > 0);
    float rs = (n_dot_i - eta * n_dot_t) / (n_dot_i + eta * n_dot_t);
    float rp = (eta * n_dot_i - n_dot_t) / (eta * n_dot_i + n_dot_t);
    float F = (rs * rs + rp * rp) / 2;
    return F;
}

/// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
/// This is a specialized version for the code above, only using the incident angle.
/// The transmission angle is derived from 
/// n_dot_i: cos(incident angle) (can be negative)
/// eta: eta_transmission / eta_incident
inline float fresnel_dielectric(float n_dot_i, float eta) {
    assert(eta > 0);
    float n_dot_t_sq = 1 - (1 - n_dot_i * n_dot_i) / (eta * eta);
    if (n_dot_t_sq <= 0) {
        // total internal reflection
        return 1.f;
    }
    float n_dot_t = sqrt(n_dot_t_sq);
    return fresnel_dielectric(fabs(n_dot_i), n_dot_t, eta);
}

inline float GTR2(float n_dot_h, float roughness) {
    float alpha = roughness * roughness;
    float a2 = alpha * alpha;
    float t = 1 + (a2 - 1) * n_dot_h * n_dot_h;
    return a2 / (M_PI * t*t);
}

inline float GGX(float n_dot_h, float roughness) {
    return GTR2(n_dot_h, roughness);
}

/// The masking term models the occlusion between the small mirrors of the microfacet models.
/// See Eric Heitz's paper "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"
/// for a great explanation.
/// https://jcgt.org/published/0003/02/03/paper.pdf
/// The derivation is based on Smith's paper "Geometrical shadowing of a random rough surface".
/// Note that different microfacet distributions have different masking terms.
inline float smith_masking_gtr2(const Vector3f &v_local, float roughness) {
    float alpha = roughness * roughness;
    float a2 = alpha * alpha;
    Vector3f v2(v_local.x()*v_local.x(),v_local.y()*v_local.y(),v_local.z()*v_local.z());
    float Lambda = (-1.f + sqrt(1.f + (v2.x() * a2 + v2.y() * a2) / v2.z())) / 2.f;
    return 1.f / (1.f + Lambda);
}

inline float smith_masking_gtr2_anisotropy(const Vector3f &v_local, float alphaX,float alphaY) {
    float ax2 =  alphaX * alphaX;
    float ay2 = alphaY * alphaY;
    Vector3f v2(v_local.x()*v_local.x(),v_local.y()*v_local.y(),v_local.z()*v_local.z());
    float Lambda = (-1.f + sqrt(1.f + (v2.x() * ax2 + v2.y() * ay2) / v2.z())) / 2.f;
    return 1.f / (1.f + Lambda);
}

/// See "Sampling the GGX Distribution of Visible Normals", Heitz, 2018.
/// https://jcgt.org/published/0007/04/01/
inline Vector3f sample_visible_normals(const Vector3f &local_dir_in, float alpha, const Point2f &rnd_param) {
    // The incoming direction is in the "ellipsodial configuration" in Heitz's paper
    if (local_dir_in.z() < 0) {
        // Ensure the input is on top of the surface.
        return -sample_visible_normals(-local_dir_in, alpha, rnd_param);
    }

    // Transform the incoming direction to the "hemisphere configuration".
    Vector3f hemi_dir_in = 
        Vector3f(alpha * local_dir_in.x(), alpha * local_dir_in.y(), local_dir_in.z()).normalized();

    // Parameterization of the projected area of a hemisphere.
    // First, sample a disk.
    float r = sqrt(rnd_param.x());
    float phi = 2 * M_PI * rnd_param.y();
    float t1 = r * cos(phi);
    float t2 = r * sin(phi);
    // Vertically scale the position of a sample to account for the projection.
    float s = (1.f + hemi_dir_in.z()) / 2.f;
    t2 = (1.f - s) * sqrt(1.f - t1 * t1) + s * t2;
    // Point in the disk space
    Vector3f disk_N(t1, t2, sqrt(std::max(0.f, 1.f - t1*t1 - t2*t2)));

    // Reprojection onto hemisphere -- we get our sampled normal in hemisphere space.
    Frame hemi_frame(hemi_dir_in);
    Vector3f hemi_N = hemi_frame.toWorld(disk_N);

    // Transforming the normal back to the ellipsoid configuration
    return Vector3f(alpha * hemi_N.x(), alpha * hemi_N.y(), std::max(0.f, hemi_N.z())).normalized();
}



NORI_NAMESPACE_END
#endif