#include <kombu/medium.h>
#include <nori/shape.h>
NORI_NAMESPACE_BEGIN
using namespace std;

class HomogeneousMedium : public Medium {
public:
    HomogeneousMedium(const PropertyList &propList){
        m_albedo = propList.getColor("albedo", Color3f(0.5f));
        m_sigma_t = propList.getColor("sigma_t", Color3f(0.5f));
        m_scale = propList.getFloat("scale", 1.f);
        m_sigma_t *= m_scale;
        m_sigma_s = m_albedo * m_sigma_t;
        m_sigma_a = m_sigma_t - m_sigma_s;
    }

    bool sample_intersection(MediumQueryRecord &mRec, Sampler* sampler) const {        // sample distance
        Vector3f direction = -mRec.wi;
        int channel = std::min((int) (sampler->next1D()* 3.f), 2);
        float density = m_sigma_t[channel];
        if(density<Epsilon){
            mRec.p = mRec.ref + mRec.tMax * direction;
            mRec.albedo = 1;
            return false; 
        }
        float t = -log(1.f - sampler->next1D()) / density;

        float pdf_failure = 0;
        float pdf_success = 0;
        float sampled_distance = t < mRec.tMax ? t :  mRec.tMax;
        for (int i=0; i<3; ++i) {
            float tmp = exp(-m_sigma_t[i] * sampled_distance);
            pdf_failure += tmp;
            pdf_success += m_sigma_t[i] * tmp;
        }
        pdf_success /= 3.f;
        pdf_failure /= 3.f;
        
        Color3f transmittance = (-m_sigma_t * sampled_distance).exp();
        if(t < mRec.tMax){
            mRec.p = mRec.ref + t * direction;
            if (mRec.p == mRec.ref) return false;
            mRec.albedo = m_sigma_s * transmittance / pdf_success;            

        } else {
            mRec.p = mRec.ref + mRec.tMax * direction;
            mRec.albedo =  transmittance / pdf_failure;            
        }
        return t < mRec.tMax;
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
    Color3f &getSigmaA() {return m_sigma_a;}
    Color3f &getSigmaS() {return m_sigma_s;}
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
                    "  sigmaA = %s,\n"
                    "  sigmaS = %s,\n"
                    "  sigmaT = %s,\n"
                    "  albedo = %s,\n"
                    "]",
                    m_sigma_a,
                    m_sigma_s,
                    m_sigma_t,
                    m_albedo);
    }

protected:
    Color3f m_sigma_a;
    Color3f m_sigma_s;
    Color3f m_sigma_t;
    Color3f m_albedo;
    float m_scale;
};

NORI_REGISTER_CLASS(HomogeneousMedium, "homogeneous");
NORI_NAMESPACE_END