#include <kombu/medium.h>
#include <nori/shape.h>
#include <boost/algorithm/string.hpp>
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
        m_albedo = propList.getColor("albedo", Color3f(0.5f));
        m_sigma_t = propList.getColor("sigma_t", Color3f(0.5f));
        m_scale = propList.getFloat("scale", 1.f);
        
        string method_string = boost::to_lower_copy(propList.getString("method", "delta"));
        if(method_string == "delta")
            m_method = EDeltaTracking;
        else if(method_string == "ratio")
            m_method = ERatioTracking;
        else if(method_string == "")
            m_method = ERayMarching;
        
    }

    bool sample_intersection(MediumQueryRecord &mRec, const float &sample) const {        // sample distance
        float t = -log(1.f - sample) / m_sigma_t.maxCoeff();
        if(t < mRec.tMax){
            mRec.p = mRec.ref + t * mRec.wi;
            if (mRec.p == mRec.ref) return false;
            mRec.albedo = m_albedo;            
            return true;
        }
        mRec.p = mRec.ref + mRec.tMax * mRec.wi;
        return false;
    }

    Color3f eval(const MediumQueryRecord &mRec) const {
        return m_albedo;
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
        Color3f sigma_t_t = -m_sigma_t * t;
        return {exp(sigma_t_t[0]),exp(sigma_t_t[1]), exp(sigma_t_t[2]) };
    }
    
    bool isHomogeneous() const override { return true; }
    
    bool contains(Point3f &p) {
        if (!m_shape)
			throw NoriException(
					"There is no shape attached to this medium!");
		return m_shape->getBoundingBox().contains(p, 0);
    }
    
    Color3f &getSigmaT() {return m_sigma_t;}
    Color3f &getAlbedo() {return m_albedo;}

    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case EPhaseFunction:
                if (m_phase)
                    throw NoriException(
                        "Medium: tried to register multiple Phase functions!");
                m_phase = static_cast<PhaseFunction *>(obj);
                break;

            default:
                throw NoriException("Medium::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }
    
    std::string toString() const override{
        return tfm::format(
                    "HomogeneousMedium[\n"
                    "  sigmaT = %s,\n"
                    "  albedo = %s,\n"
                    "]",
                    m_sigma_t,
                    m_albedo);
    }

protected:
    Color3f m_sigma_t;
    Color3f m_albedo;
    float m_scale;
    EIntegrationMethod m_method;
};

NORI_REGISTER_CLASS(HeterogeneousMedium, "heterogeneous");
NORI_NAMESPACE_END