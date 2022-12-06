#include <kombu/medium.h>
#include <kombu/volume.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class HeterogeneousMedium : public Medium {
public:
    /// Possible free path sampling & transmittence estimation modes
    enum EIntegrationMethod {
        EDeltaTracking = 0,
        ERatioTracking,
        ERayMarching
    };
    
    HeterogeneousMedium(const PropertyList &propList){
        m_scale = propList.getFloat("scale", 1.f);
        string method_string = propList.getString("method", "delta");
        if(method_string == "delta")
            m_method = EDeltaTracking;
        else if(method_string == "ratio")
            m_method = ERatioTracking;
        else if(method_string == "march")
            m_method = ERayMarching;
        
    }

    bool sample_intersection(MediumQueryRecord &mRec, Sampler* sampler) const {        // sample distance
        Vector3f direction = -mRec.wi;
        Ray3f ray(mRec.ref, direction, 0, mRec.tMax);
        float pdf_failure = 1.0f;
        float pdf_success = 1.0f;
        Color3f transmittance(1.0f);
        float mint, maxt;
        if (!m_shape->getBoundingBox().rayIntersect(ray, mint, maxt))
            return false;
        mint = std::max(mint, ray.mint);
        maxt = std::min(maxt, ray.maxt);
        if(m_method==EDeltaTracking || m_method == ERatioTracking){
            float t = mint, densityAtT = 0;
            while (true) {
                t -= log(1-sampler->next1D()) * m_invMaxDensity;
                if (t >= maxt)
                    break;

                Point3f p = ray(t);
                densityAtT = m_sigma_t->lookup(m_shape, p) * m_scale;

                if (densityAtT * m_invMaxDensity > sampler->next1D()) {
                    mRec.p = p;
                    mRec.albedo = m_albedo->lookup(m_shape, p);
                    return true;
                }
            }
        } else if(m_method==ERayMarching){
            float m_stepSize = std::min(m_sigma_t->getStepSize(m_shape), m_albedo->getStepSize(m_shape));
            if (m_stepSize == std::numeric_limits<float>::infinity())
                throw NoriException("Unable to infer a suitable step size for deterministic "
                        "integration, please specify one manually using the 'stepSize' "
                        "parameter.");
            float int_density = 0.f, desired_density = -std::log(1-sampler->next1D());;
            float length = maxt - mint, max_comp = 0;
            Point3f p = ray(mint), p_last = ray(maxt);
            // Ignore degenerate path segments (max and min too close)
            for(int i = 0; i < 3; i++){
                max_comp = std::max(std::max(max_comp, std::abs(p[i])), std::abs(p_last[i]));
            }
            if(length < 1e-6f * max_comp){
                mRec.p = mRec.ref + mRec.tMax * direction;
                mRec.albedo = 1.f; 
                return false;
            }
            // Compute a suitable step size
            int n_steps = (int) std::ceil(length / (2 * m_stepSize));
            float step_size = length / n_steps,
                  multiplier = (1.f/6.f) * step_size * m_scale;
            Vector3f full_step = ray.d * step_size, half_step = full_step * 0.5f;

            // inverse integration
            float lo = m_sigma_t->lookup(m_shape, p);
            for(int i = 0; i < n_steps; i++){
                float mi = m_sigma_t->lookup(m_shape, p + half_step),
                      hi = m_sigma_t->lookup(m_shape, p + full_step),
                      new_density = int_density + multiplier * (lo + 4*mi + hi);
                if(new_density >= desired_density){
                    // cross over the line
                    float a = 0, b = step_size, x = a,
                          residual = int_density - desired_density,
                          step_size_pow = step_size * step_size,
                          temp = m_scale / step_size_pow;
                    int it = 1;
                    while(true){
                        // Lagrange polynomial
                        float d_residual = temp * (lo * step_size_pow
                            -(3 * lo - 4 * mi + hi) * step_size * x
                            + 2 * (lo - 2*mi + hi) * x * x
                        );
                        x -= residual / d_residual;
                        if(!(x <= a || x >=b || d_residual==0)){
                            x = 0.5f * (a + b);
                        }
                        float int_poly = int_density + temp * (1.0f / 6.0f) * ( x *
                            (6 * lo * step_size_pow - 3 * (3 * lo - 4 * mi + hi) * step_size * x
                            + 4 * (lo - 2 * mi + hi) * x * x));
                        residual = int_poly - desired_density;
                        if (std::abs(residual) < 1e-6f) {
                            float t = mint + step_size * i + x;
                            int_density = int_poly;
                            mRec.p = mRec.ref + t * direction;
                            mRec.albedo = m_albedo->lookup(m_shape, mRec.p);
                            return true;
                        } else if (++it > 30) {
                            mRec.p = mRec.ref + mRec.tMax * direction;
                            mRec.albedo = 1.f; 
                            return false;
                        }
                        if(residual > 0) b = x;
                        else a = x;
                    }
                }
                Point3f next = p + full_step;
                if(p==next) break;
                int_density = new_density;
                lo = hi;
                p = next;
            }
        }
        mRec.p = mRec.ref + mRec.tMax * direction;
        mRec.albedo = 1.f; 
        return false;
    }

    Color3f eval(const MediumQueryRecord &mRec) const {
        return m_albedo->lookup(m_shape, mRec.ref);
    }

    float pdf(const MediumQueryRecord &mRec) const override {
        PhaseFunctionQueryRecord pRec(mRec.wi, mRec.wo, mRec.measure);
        float pdf_w = m_phase -> eval(pRec);
        return pdf_w;
    }

    Color3f evalTransmittance(const MediumQueryRecord &mRec, Sampler* sampler) const override {
        Vector3f direction = (mRec.p - mRec.ref).normalized();
        float t = (mRec.p - mRec.ref).norm();
        Ray3f ray = Ray3f(mRec.ref, direction, 0, t);
        float mint, maxt;
        if (!m_shape->getBoundingBox().rayIntersect(ray, mint, maxt))
            return Color3f(1.f);
        mint = std::max(mint, ray.mint);
        maxt = std::min(maxt, ray.maxt);
        if(m_method==EDeltaTracking){
            int nSamples = 2; /// XXX make configurable
            float result = 0;
            for (int i=0; i<nSamples; ++i) {
                float t = mint;
                while (true) {
                    t -= log(1-sampler->next1D()) * m_invMaxDensity;
                    if (t >= maxt) {
                        result += 1;
                        break;
                    }
                    Point3f p = ray(t);
                    float density = m_sigma_t->lookup(m_shape, p) * m_scale;

                    if (density * m_invMaxDensity > sampler->next1D())
                        break;
                }
            }
            return Color3f(result/nSamples);
        } else if(m_method == ERatioTracking){
            float t = mint;
            Color3f Tr(1.f);
            while (true) {
                t -= log(1-sampler->next1D()) * m_invMaxDensity;
                if (t >= maxt) break;
                Point3f p = ray(t);
                float densityT = m_sigma_t->lookup(m_shape, p) * m_scale;
                Tr *= (1.f - densityT * m_invMaxDensity);
            }
            return Tr;
        } else if(m_method==ERayMarching){
            float m_stepSize = std::min(m_sigma_t->getStepSize(m_shape), m_albedo->getStepSize(m_shape));
            if (m_stepSize == std::numeric_limits<float>::infinity())
                throw NoriException("Unable to infer a suitable step size for deterministic "
                        "integration, please specify one manually using the 'stepSize' "
                        "parameter.");
            float length = maxt - mint, max_comp = 0;
            Point3f p = ray(mint), p_last = ray(maxt);
            // Ignore degenerate path segments (max and min too close)
            for(int i = 0; i < 3; i++){
                max_comp = std::max(std::max(max_comp, std::abs(p[i])), std::abs(p_last[i]));
            }
            if(length < 1e-6f * max_comp) return 1.f;
             // Compute a suitable step size
            int n_steps = (int) std::ceil(length / m_stepSize);
            n_steps += n_steps % 2;
            const float step_size = length/(float)n_steps;
            const Vector3f increment = ray.d * step_size;
            float int_density = m_sigma_t->lookup(m_shape, p) + m_sigma_t->lookup(m_shape, p_last);
            const float stop_density = -std::log(Epsilon);
            const float stop_value = stop_density * 3.f / (step_size * m_scale);
            
            // integration
            p += increment;
            
            float weight = 4.0;
            for(int i = 1; i < n_steps; ++i){
                int_density += weight * m_sigma_t->lookup(m_shape, p);
                weight = 6.0 - weight;
                if(int_density > stop_value) return 0.f;
                Point3f next = p + increment;
                if(p==next) break;
                p = next;
            }
            int_density = int_density * m_scale * step_size *( 1.f/ 3.f);
            return Color3f(expf(-int_density));
        }
    }
    
    
    bool contains(Point3f &p) {
        if (!m_shape)
			throw NoriException(
					"There is no shape attached to this medium!");
		return m_shape->getBoundingBox().contains(p, 0);
    }

    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case EPhaseFunction:
                if (m_phase)
                    throw NoriException(
                        "Medium: tried to register multiple Phase functions!");
                m_phase = static_cast<PhaseFunction *>(obj);
                break;
            case EVolume:
                // m_phase = static_cast<Volume *>(obj);
                if(obj->getIdName() == "albedo"){
                    m_albedo = static_cast<Volume *>(obj);
                } else if(obj->getIdName() == "sigma_t"){
                    m_sigma_t = static_cast<Volume *>(obj);
                    m_maxDensity = m_scale * m_sigma_t->getMaximumValue();
                    m_invMaxDensity = 1.0f/m_maxDensity;
                }
                break;
            default:
                throw NoriException("Medium::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }
    
    void setShape(Shape * shape) override{
        m_shape = shape; 
    }
    
    std::string toString() const override{
        return tfm::format(
                    "HeterogeneousMedium[\n"
                    "  scale = %s,\n"
                    "  sigmaT = %s,\n"
                    "  albedo = %s,\n"
                    
                    "]",
                    m_scale,
                    m_sigma_t->toString(),
                    m_albedo->toString());
    }
    
   
protected:
    
    Volume* m_sigma_t = nullptr;
    Volume* m_albedo = nullptr;
    float m_scale;
    EIntegrationMethod m_method;
    float m_maxDensity;
    float m_invMaxDensity;
    Transform m_worldToGrid;
};

NORI_REGISTER_CLASS(HeterogeneousMedium, "heterogeneous");
NORI_NAMESPACE_END