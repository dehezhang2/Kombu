#include <kombu/medium.h>

NORI_NAMESPACE_BEGIN
class HomogeneousMedium : public Medium {
public:
    HomogeneousMedium(const PropertyList &propList){
        m_sigma_a = propList.getColor("sigma_a", Color3f(0.5f));
        m_sigma_s = propList.getColor("sigma_s", Color3f(0.5f));
        m_sigma_t = m_sigma_a + m_sigma_s;
        m_albedo =  m_sigma_s / m_sigma_t;
    }

    Color3f sample(MediumQueryRecord &mRec, const Point2f &sample){
        return Color3f(0.0f);
    }

    Color3f eval(const MediumQueryRecord &mRec){
        return Color3f(0.0f);
    }

    float pdf(const MediumQueryRecord &mRec){
        return 0.0f;
    }

    Color3f evalTransmittance(const MediumQueryRecord &mRec){
        return 0.0f;
    }
    
    bool isHomogeneous() const { return true; }
    
    bool contains(Point3f &p){
        if (!m_shape)
			throw NoriException(
					"There is no shape attached to this medium!");
		return m_shape->getBoundingBox().contains(p, 0);
    }
    Color3f &getSigmaA(){return m_sigma_a;}
    Color3f &getSigmaS(){return m_sigma_s;}
    Color3f &getSigmaT(){return m_sigma_t;}
    Color3f &getAlbedo(){return m_albedo;}

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
    std::string toString() const {
        return tfm::format(
                    "HomogeneousMedium[\n"
                    "  sigmaA = %s,\n"
                    "  sigmaS = %s,\n"
                    "  sigmaT = %s,\n"
                    "]",
                    m_sigma_a,
                    m_sigma_s,
                    m_sigma_t);
    }

protected:
    Color3f m_sigma_a;
    Color3f m_sigma_s;
    Color3f m_sigma_t;
    Color3f m_albedo;
};

NORI_REGISTER_CLASS(HomogeneousMedium, "homogeneous");
NORI_NAMESPACE_END