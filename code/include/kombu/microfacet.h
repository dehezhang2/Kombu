#if !defined(__KOMBU_MICROFACET_H)
#define __KOMBU_MICROFACET_H

#include <nori/object.h>
#include <nori/frame.h>
#include <kombu/principledhelpers.h>

NORI_NAMESPACE_BEGIN

enum class MicrofacetType : uint32_t {
    /// Beckmann distribution derived from Gaussian random surfaces
    Beckmann = 0, // not implemented

    /// GGX: Long-tailed distribution for very rough surfaces (aka. Trowbridge-Reitz distr.)
    GGX = 1
};



inline Point2f square_to_uniform_disk_concentric(const Point2f&sample) {

    float x = fmadd(2.f, sample.x(), -1.f),
          y = fmadd(2.f, sample.y(), -1.f);

    /* Modified concentric map code with less branching (by Dave Cline), see
       http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html

      Original non-vectorized version:

        Value phi, r;
        if (x == 0 && y == 0) {
            r = phi = 0;
        } else if (x * x > y * y) {
            r = x;
            phi = (dr::Pi / 4.f) * (y / x);
        } else {
            r = y;
            phi = (dr::Pi / 2.f) - (x / y) * (dr::Pi / 4.f);
        }
    */

    bool is_zero         =  x==0.f&&
                            y==0.f,
         quadrant_1_or_3 = abs(x) < abs(y);

    float r  = (quadrant_1_or_3? y: x),
          rp = (quadrant_1_or_3? x: y);

    float phi = 0.25f * M_PI* rp / r;
    if (quadrant_1_or_3){
        phi = 0.5f * M_PI- phi;
    }
    if (is_zero){
        phi = 0.f;
    }
    float s = sin(phi);
    float c = cos(phi);
    return { r * c, r * s };
}

// template <typename float, typename Spectrum>
class MicrofacetDistribution {
public:

    /**
     * Create an isotropic microfacet distribution of the specified type
     *
     * \param type
     *     The desired type of microfacet distribution
     * \param alpha
     *     The surface roughness
     */
    MicrofacetDistribution(MicrofacetType type, float alpha, bool sample_visible = true)
        : m_type(type), m_alpha_u(alpha), m_alpha_v(alpha),
          m_sample_visible(sample_visible) {
        configure();
    }

    /**
     * Create an anisotropic microfacet distribution of the specified type
     *
     * \param type
     *     The desired type of microfacet distribution
     * \param alpha_u
     *     The surface roughness in the tangent direction
     * \param alpha_v
     *     The surface roughness in the bitangent direction
     */
    MicrofacetDistribution(MicrofacetType type, float alpha_u, float alpha_v,
                           bool sample_visible = true)
        : m_type(type), m_alpha_u(alpha_u), m_alpha_v(alpha_v),
          m_sample_visible(sample_visible) {
        configure();
    }


public:
    /// Return the distribution type
    MicrofacetType type() const { return m_type; }

    /// Return the roughness (isotropic case)
    float alpha() const { return m_alpha_u; }

    /// Return the roughness along the tangent direction
    float alpha_u() const { return m_alpha_u; }

    /// Return the roughness along the bitangent direction
    float alpha_v() const { return m_alpha_v; }

    /// Return whether or not only visible normals are sampled?
    bool sample_visible() const { return m_sample_visible; }

    /// Is this an isotropic microfacet distribution?
    bool is_isotropic() const {
        return m_alpha_u == m_alpha_v;
    }

    /// Is this an anisotropic microfacet distribution?
    bool is_anisotropic() const { return m_alpha_u != m_alpha_v; }

    /// Scale the roughness values by some constant
    void scale_alpha(float value) {
        m_alpha_u *= value;
        m_alpha_v *= value;
    }

    /**
     * \brief Evaluate the microfacet distribution function
     *
     * \param m
     *     The microfacet normal
     */
    float eval(const Vector3f &m) const {
        float alpha_uv = m_alpha_u * m_alpha_v,
              cos_theta         = Frame::cosTheta(m),
              cos_theta_2       = cos_theta*cos_theta,
              result;

        if (m_type == MicrofacetType::Beckmann) {
            // Beckmann distribution function for Gaussian random surfaces
            result = exp(-((m.x() / m_alpha_u)*(m.x() / m_alpha_u) +
                               (m.y() / m_alpha_v)*(m.y() / m_alpha_v)) /
                             cos_theta_2) /
                     (M_PI * alpha_uv * cos_theta_2*cos_theta_2);
        } else {
            // GGX / Trowbridge-Reitz distribution function
            result =
                1.0f/(M_PI * alpha_uv *
                        pow(pow(m.x() / m_alpha_u,2) +
                            pow(m.y() / m_alpha_v,2) + pow(m.z(),2),2) );
        }

        // Prevent potential numerical issues in other stages of the model
        return (result * cos_theta > 1e-20f)? result: 0.f;
    }

    /**
     * \brief Returns the density function associated with
     * the \ref sample() function.
     *
     * \param wi
     *     The incident direction (only relevant if visible normal sampling is used)
     *
     * \param m
     *     The microfacet normal
     */
    float pdf(const Vector3f &wi, const Vector3f &m) const {
        float result = eval(m);

        if (m_sample_visible)
            result *= smith_g1(wi, m) * abs(wi.dot(m)) / Frame::cosTheta(wi);
        else
            result *= Frame::cosTheta(m);

        return result;
    }

    /**
     * \brief Draw a sample from the microfacet normal distribution
     *  and return the associated probability density
     *
     * \param wi
     *    The incident direction. Only used if
     *    visible normal sampling is enabled.
     *
     * \param sample
     *    A uniformly distributed 2D sample
     *
     * \return A tuple consisting of the sampled microfacet normal
     *         and the associated solid angle density
     */
    std::pair<Normal3f, float> sample(const Vector3f &wi,
                                      const Point2f &sample) const {
        if (!m_sample_visible) {
            float sin_phi, cos_phi, cos_theta, cos_theta_2, alpha_2, pdf;

            // Sample azimuth component (identical for Beckmann & GGX)
            if (is_isotropic()) {
                sin_phi = sin((2.f * M_PI) * sample.y());
                cos_phi = cos((2.f * M_PI) * sample.y());
                alpha_2 = m_alpha_u * m_alpha_u;
            } else {
                float ratio  = m_alpha_v / m_alpha_u,
                      tmp    = ratio * tan((2.f * M_PI) * sample.y());

                cos_phi = rsqrt(fmadd(tmp, tmp, 1.f));

                if (abs(sample.y() - .5f) - .25f >=0){
                    cos_phi = -cos_phi;
                }

                sin_phi = cos_phi * tmp;

                alpha_2 = 1.0f/((cos_phi / m_alpha_u)*(cos_phi / m_alpha_u) +
                                (sin_phi / m_alpha_v)*(sin_phi / m_alpha_v));
            }

            // Sample elevation component
            if (m_type == MicrofacetType::Beckmann) {
                // Beckmann distribution function for Gaussian random surfaces
                cos_theta = rsqrt(fmadd(-alpha_2, log(1.f - sample.x()), 1.f));
                cos_theta_2 = cos_theta*cos_theta;

                // Compute probability density of the sampled position
                float cos_theta_3 = std::max(cos_theta_2 * cos_theta, 1e-20f);
                pdf = (1.f - sample.x()) / (M_PI * m_alpha_u * m_alpha_v * cos_theta_3);
            } else {
                // GGX / Trowbridge-Reitz distribution function
                float tan_theta_m_2 = alpha_2 * sample.x() / (1.f - sample.x());
                cos_theta = rsqrt(1.f + tan_theta_m_2);
                cos_theta_2 = cos_theta*cos_theta;

                // Compute probability density of the sampled position
                float temp = 1.f + tan_theta_m_2 / alpha_2,
                      cos_theta_3 = std::max(cos_theta_2 * cos_theta, 1e-20f);
                pdf = 1.0f/(M_PI * m_alpha_u * m_alpha_v * cos_theta_3 * temp*temp);
            }

            float sin_theta = sqrt(1.f - cos_theta_2);

            return {
                Normal3f(cos_phi * sin_theta,
                         sin_phi * sin_theta,
                         cos_theta),
                pdf
            };
        } else {
            // Visible normal sampling.
            float sin_phi, cos_phi, cos_theta;

            // Step 1: stretch wi
            Vector3f wi_p = Vector3f(
                m_alpha_u * wi.x(),
                m_alpha_v * wi.y(),
                wi.z()
            ).normalized();

            sin_phi = Frame::sinPhi(wi_p);
            cos_phi = Frame::cosPhi(wi_p);
            cos_theta = Frame::cosTheta(wi_p);

            // Step 2: simulate P22_{wi}(slope.x, slope.y, 1, 1)
            Vector2f slope = sample_visible_11(cos_theta, sample);

            // Step 3: rotate & unstretch
            slope = Vector2f(
                fmadd(cos_phi, slope.x(), -sin_phi * slope.y()) * m_alpha_u,
                fmadd(sin_phi, slope.x(), cos_phi * slope.y()) * m_alpha_v);

            // Step 4: compute normal & PDF
            Normal3f m = Vector3f(-slope.x(), -slope.y(), 1).normalized();

            float pdf = eval(m) * smith_g1(wi, m) * abs(wi.dot(m)) /
                        Frame::cosTheta(wi);

            return { m, pdf };
        }
    }

    /// Smith's separable shadowing-masking approximation
    float G(const Vector3f &wi, const Vector3f &wo, const Vector3f &m) const {
        return smith_g1(wi, m) * smith_g1(wo, m);
    }

    /**
     * \brief Smith's shadowing-masking function for a single direction
     *
     * \param v
     *     An arbitrary direction
     * \param m
     *     The microfacet normal
     */
    float smith_g1(const Vector3f &v, const Vector3f &m) const {
        float xy_alpha_2 = sqr(m_alpha_u * v.x()) + sqr(m_alpha_v * v.y()),
              tan_theta_alpha_2 = xy_alpha_2 / sqr(v.z()),
              result;

        if (m_type == MicrofacetType::Beckmann) {
            float a = rsqrt(tan_theta_alpha_2), a_sqr = sqr(a);
            /* Use a fast and accurate (<0.35% rel. error) rational
               approximation to the shadowing-masking function */
            result = (a >= 1.6f ? 
                      1.f:
                      (3.535f * a + 2.181f * a_sqr) / (1.f + 2.276f * a + 2.577f * a_sqr));
        } else {
            result = 2.f / (1.f + sqrt(1.f + tan_theta_alpha_2));
        }

        // Perpendicular incidence -- no shadowing/masking
        if (xy_alpha_2==0.f){
            result = 1.0f;
        }

        /* Ensure consistent orientation (can't see the back
           of the microfacet from the front and vice versa) */
        if(v.dot(m)*Frame::cosTheta(v)<=0.f){
            result = 0.f;
        }
        return result;
    }

    /// \brief Visible normal sampling code for the alpha=1 case
    Vector2f sample_visible_11(float cos_theta_i, Point2f sample) const {
        if (m_type == MicrofacetType::Beckmann) {
            //TODO
        } else {
            // Choose a projection direction and re-scale the sample
            Point2f p = square_to_uniform_disk_concentric(sample);

            float s = 0.5f * (1.f + cos_theta_i);
            p.y() = disney_lerp(safe_sqrt(1.f - sqr(p.x())), p.y(), s);

            // Project onto chosen side of the hemisphere
            float x = p.x(), y = p.y(),
                  z = safe_sqrt(1.f - p.squaredNorm());

            // Convert to slope
            float sin_theta_i = safe_sqrt(1.f - sqr(cos_theta_i));
            float norm = 1.f/(fmadd(sin_theta_i, y, cos_theta_i * z));
            return Vector2f(fmadd(cos_theta_i, y, -sin_theta_i * z), x) * norm;
        }
    }


protected:
    void configure() {
        m_alpha_u = std::max(m_alpha_u, 1e-4f);
        m_alpha_v = std::max(m_alpha_v, 1e-4f);
    }

    /// Compute the squared 1D roughness along direction \c v
    float project_roughness_2(const Vector3f &v) const {
        if (is_isotropic())
            return sqr(m_alpha_u);

        float sin_phi_2, cos_phi_2;
        sin_phi_2 = Frame::sinPhi2(v);
        cos_phi_2 = Frame::cosPhi2(v);

        return sin_phi_2 * sqr(m_alpha_v) + cos_phi_2 * sqr(m_alpha_u);
    }

protected:
    MicrofacetType m_type;
    float m_alpha_u, m_alpha_v;
    bool  m_sample_visible;
};







NORI_NAMESPACE_END
#endif