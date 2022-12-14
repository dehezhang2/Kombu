
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

class BlendBSDF : public BSDF {
public:
    BlendBSDF(const PropertyList &prop) { 
        if(prop.has("weight")) {
            PropertyList l;
            l.setColor("value", prop.getFloat("weight"));
            m_weight.reset(static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l)));
        }
        is_fresnel = prop.getBoolean("fresnel", false);
        if(is_fresnel){
            /* Interior IOR (default: BK7 borosilicate optical glass) */
            m_intIOR = prop.getFloat("intIOR", 1.450f);
            /* Exterior IOR (default: air) */
            m_extIOR = prop.getFloat("extIOR", 1.000277f);
        }
    }

    virtual void addChild(NoriObject *obj) override {
        switch (obj->getClassType()) {
            case ETexture:
                if (m_weight.get())
                    throw NoriException("There is already a weight defined!");
                m_weight.reset(static_cast<Texture<Color3f> *>(obj));
                break;
            case EBSDF:
                if(obj->getIdName() == "bsdf1") {
                    if (m_bsdf1.get())
                        throw NoriException("There is already an BSDF1 defined!");
                    m_bsdf1.reset(static_cast<BSDF*>(obj));
                }else if(obj->getIdName() == "bsdf2") {
                    if (m_bsdf2.get())
                        throw NoriException("There is already an BSDF2 defined!");
                    m_bsdf2.reset(static_cast<BSDF*>(obj));
                }
                break;
            default:
                throw NoriException("BlendBSDF::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    virtual void activate() override {
        if(!m_weight.get()) {
            PropertyList l;
            l.setColor("value", Color3f(0.5f));
            m_weight.reset(static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l)));
        }
        if(!m_bsdf1.get()) {
            PropertyList l;
            l.setColor("albedo", 0.5f);
            m_bsdf1.reset(static_cast<BSDF *>(NoriObjectFactory::createInstance("diffuse", l)));
            m_bsdf1->activate();
        }
        if(!m_bsdf2.get()) {
            PropertyList l;
            l.setColor("albedo", 0.5f);
            m_bsdf2.reset(static_cast<BSDF *>(NoriObjectFactory::createInstance("diffuse", l)));
            m_bsdf2->activate();
        }
    }

    virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
        float weight = is_fresnel ? fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR): m_weight->eval(bRec.uv).getLuminance();
        return m_bsdf1->eval(bRec) * (1 - weight) +
               m_bsdf2->eval(bRec) * weight;
    }


    virtual float pdf(const BSDFQueryRecord &bRec) const override {
        float weight = is_fresnel ? fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR): m_weight->eval(bRec.uv).getLuminance();
        return m_bsdf1->pdf(bRec) * (1 - weight) +
               m_bsdf2->pdf(bRec) * weight;
    }

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample2d) const override {
        float weight = is_fresnel ? fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR): m_weight->eval(bRec.uv).getLuminance();

        if (sample2d.x()>weight){
            Point2f new_sample((sample2d.x()-weight)/(1.f-weight),sample2d.y());
            return m_bsdf1->sample(bRec,new_sample);
        }else{
            Point2f new_sample(sample2d.x()/weight,sample2d.y());
            return m_bsdf2->sample(bRec,new_sample);
        }
    }


    virtual std::string toString() const override {
        std::ostringstream oss;
        oss << "BlendBSDF[" << std::endl
            << "  weight = " << indent(m_weight->toString(),2) << "," << std::endl
            << "  nested_bsdf1 = " << indent(m_bsdf1->toString(),2) << "," << std::endl
            << "  nested_bsdf2 = " << indent(m_bsdf2->toString(),2) << std::endl
            << "]";
        return oss.str();
    }




private:
    std::shared_ptr<BSDF> m_bsdf1;
    std::shared_ptr<BSDF> m_bsdf2;
    std::shared_ptr<Texture<Color3f>> m_weight;
    bool is_fresnel;
    float m_intIOR;
    float m_extIOR;
};

NORI_REGISTER_CLASS(BlendBSDF, "blendbsdf");
NORI_NAMESPACE_END
