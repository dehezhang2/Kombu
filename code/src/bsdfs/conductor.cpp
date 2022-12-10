
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>
#include <Eigen/Geometry>
NORI_NAMESPACE_BEGIN

inline float safe_sqrt(float a){
    return sqrt(std::max(a,0.f));
}

inline float fresnel_conductor(float cos_theta_i, float eta_r, float eta_i ) {
    // Modified from "Optics" by K.D. Moeller, University Science Books, 1988
    float cos_theta_i_2 = cos_theta_i * cos_theta_i,
          sin_theta_i_2 = 1.f - cos_theta_i_2,
          sin_theta_i_4 = sin_theta_i_2 * sin_theta_i_2;


    float temp_1   = eta_r * eta_r - eta_i * eta_i - sin_theta_i_2,
          a_2_pb_2 = safe_sqrt(temp_1*temp_1 + 4.f * eta_i * eta_i * eta_r * eta_r),
          a        = safe_sqrt(.5f * (a_2_pb_2 + temp_1));

    float term_1 = a_2_pb_2 + cos_theta_i_2,
          term_2 = 2.f * cos_theta_i * a;

    float r_s = (term_1 - term_2) / (term_1 + term_2);

    float term_3 = a_2_pb_2 * cos_theta_i_2 + sin_theta_i_4,
          term_4 = term_2 * sin_theta_i_2;

    float r_p = r_s * (term_3 - term_4) / (term_3 + term_4);

    return 0.5f * (r_s + r_p);
}

class ConductorBSDF : public BSDF {
public:
    ConductorBSDF(const PropertyList &prop) { 
        if(prop.has("eta")) {
            PropertyList l;
            l.setColor("value", prop.getFloat("eta"));
            m_eta.reset(static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l)));
        }
        if(prop.has("k")) {
            PropertyList l;
            l.setColor("value", prop.getFloat("k"));
            m_k.reset(static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l)));
        }
        if(prop.has("specular_reflectance")) {
            PropertyList l;
            l.setColor("value", prop.getColor("specular_reflectance"));
            m_specular_reflectance.reset(static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l)));
        }
    }

    virtual void addChild(NoriObject *obj) override {
        switch (obj->getClassType()) {
            case ETexture:
                if(obj->getIdName() == "eta") {
                    if (m_eta.get())
                        throw NoriException("There is already an eta defined!");
                    m_eta.reset(static_cast<Texture<Color3f> *>(obj));
                }
                else if(obj->getIdName() == "k") {
                    if (m_k.get())
                        throw NoriException("There is already a k defined!");
                    m_k.reset(static_cast<Texture<Color3f> *>(obj));
                }
                else if(obj->getIdName() == "specular_reflectance") {
                    if (m_specular_reflectance.get())
                        throw NoriException("There is already an spec_tint defined!");
                    m_specular_reflectance.reset(static_cast<Texture<Color3f> *>(obj));
                }
                break;

            default:
                throw NoriException("ConductorBSDF::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    virtual void activate() override {
        if(!m_eta.get()) {
            PropertyList l;
            l.setColor("value", Color3f(1.15f));
            m_eta.reset(static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l)));
        }
        if(!m_k.get()) {
            PropertyList l;
            l.setColor("value", Color3f(2.504f));
            m_k.reset(static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l)));
        }
        if(!m_specular_reflectance.get()) {
            PropertyList l;
            l.setColor("value", Color3f(1.0f));
            m_k.reset(static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l)));
        }
    }

    virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
        return 0.f;
    }


    virtual float pdf(const BSDFQueryRecord &bRec) const override {
        return 0.f;
    }

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample2d) const override {
        bRec.measure = EDiscrete;
        float cos_theta_i = Frame::cosTheta(bRec.wi);
        bool active = cos_theta_i > 0.f;

        Color3f value(0.f);
        if (!active)
            return value;
        bRec.eta = 1.f;

        Vector3f wo = -bRec.wi;
        wo.z() = bRec.wi.z();
        bRec.wo  = wo;
        float eta = m_eta->eval(bRec.uv).getLuminance();
        float k = m_k->eval(bRec.uv).getLuminance();
        Color3f reflectance = m_specular_reflectance->eval(bRec.uv);
        value = reflectance * fresnel_conductor(cos_theta_i, eta, k);
        return value;
    }


    virtual std::string toString() const override {
        std::ostringstream oss;
        oss << "ConductorBSDF[" << std::endl
            << "  eta = " << indent(m_eta->toString(),2) << "," << std::endl
            << "  k = " << indent(m_k->toString(),2) << "," << std::endl
            << "  specular_reflectance = " << indent(m_specular_reflectance->toString(),2) << std::endl
            << "]";
        return oss.str();
    }




private:
    std::shared_ptr<Texture<Color3f>> m_eta;
    std::shared_ptr<Texture<Color3f>> m_k;
    std::shared_ptr<Texture<Color3f>> m_specular_reflectance;
};

NORI_REGISTER_CLASS(ConductorBSDF, "conductor");
NORI_NAMESPACE_END
