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
        if(m_method==EDeltaTracking){
            float pdf_failure = 1.0f;
            float pdf_success = 1.0f;
            Color3f transmittance(1.0f);

            float mint, maxt;
            
            if (!m_shape->getBoundingBox().rayIntersect(ray, mint, maxt))
                return false;
            mint = std::max(mint, ray.mint);
            maxt = std::min(maxt, ray.maxt);

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
        } else if(m_method == ERatioTracking){

        } else if(m_method==ERayMarching){

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
        if(m_method==EDeltaTracking){
            float mint, maxt;
            
            if (!m_shape->getBoundingBox().rayIntersect(ray, mint, maxt))
                return Color3f(1.f);
            mint = std::max(mint, ray.mint);
            maxt = std::min(maxt, ray.maxt);
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

        } else if(m_method==ERayMarching){

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
    Volume* m_sigma_t;
    Volume* m_albedo;
    float m_scale;
    EIntegrationMethod m_method;
    float m_maxDensity;
    float m_invMaxDensity;
    Transform m_worldToGrid;
};

NORI_REGISTER_CLASS(HeterogeneousMedium, "heterogeneous");
NORI_NAMESPACE_END