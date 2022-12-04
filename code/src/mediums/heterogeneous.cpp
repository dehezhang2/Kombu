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
        else if(method_string == "")
            m_method = ERayMarching;
        
    }

    bool sample_intersection(MediumQueryRecord &mRec, Sampler* sampler) const {        // sample distance
        float sigma_t = m_sigma_t->lookup(Vector3f(0.f));
        float albedo = m_albedo->lookup(Vector3f(0.f));
        // float sigma_s = albedo * sigma_t;
        Vector3f direction = -mRec.wi;
        float density = sigma_t;

        float t = -log(1.f - sampler->next1D()) / density;


        if(t < mRec.tMax){
            mRec.p = mRec.ref + t * direction;
            if (mRec.p == mRec.ref) return false;
            mRec.albedo = albedo; 
        } else {
            mRec.p = mRec.ref + mRec.tMax * direction;
            mRec.albedo =  1.f;            
        }
        return t < mRec.tMax;
    }

    Color3f eval(const MediumQueryRecord &mRec) const {
        return m_albedo->lookup(Vector3f(0.f));
    }

    float pdf(const MediumQueryRecord &mRec) const override {
        PhaseFunctionQueryRecord pRec(mRec.wi, mRec.wo, mRec.measure);
        float pdf_w = m_phase -> eval(pRec);
        // float pdf_x = evalTransmittance(mRec).getLuminance();
        // return pdf_x * pdf_w;
        return pdf_w;
    }

    Color3f evalTransmittance(const MediumQueryRecord &mRec, Sampler* sampler) const override {
        float t = (mRec.p - mRec.ref).norm();
        Color3f sigma_t_t = -m_sigma_t->lookup(Vector3f(0.f)) * t;
        return {exp(sigma_t_t[0]),exp(sigma_t_t[1]), exp(sigma_t_t[2]) };
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
                    m_albedo->setBoundingBox(m_shape->getBoundingBox());
                } else if(obj->getIdName() == "sigma_t"){
                    m_sigma_t = static_cast<Volume *>(obj);
                    m_sigma_t->setBoundingBox(m_shape->getBoundingBox());
                }
                break;
            default:
                throw NoriException("Medium::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }
    
    std::string toString() const override{
        return tfm::format(
                    "HeterogeneousMedium[\n"
                    "  scale = %s,\n"
                    "  sigmaT = %s,\n"
                    "  albedo = %s,\n"
                    
                    "]",
                    m_scale,
                    m_sigma_t,
                    m_albedo);
    }

protected:
    Volume* m_sigma_t;
    Volume* m_albedo;
    float m_scale;
    EIntegrationMethod m_method;
};

NORI_REGISTER_CLASS(HeterogeneousMedium, "heterogeneous");
NORI_NAMESPACE_END