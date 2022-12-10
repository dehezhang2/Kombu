
#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

class Bumpmap : public BSDF {
public:
    Bumpmap(const PropertyList &prop) { 
        m_enhance = prop.getFloat("enhance", 1.0f);
    }

    virtual void addChild(NoriObject *obj) override {
        switch (obj->getClassType()) {
            case ETexture:
                if (m_texture.get())
                    throw NoriException("There is already a Texture defined!");
                m_texture.reset(static_cast<Texture<Color3f> *>(obj));
                break;
            case EBSDF:
                if (m_bsdf.get())
                    throw NoriException("There is already an BSDF defined!");
                m_bsdf.reset(static_cast<BSDF*>(obj));
                break;
            default:
                throw NoriException("Bumpmap::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    virtual void activate() override {
        if(!m_texture.get()) {
            PropertyList l;
            l.setColor("value", Color3f(0.5f));
            l.setInteger("boundary",0);
            m_texture.reset(static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l)));
        }
        if(!m_bsdf.get()) {
            PropertyList l;
            l.setColor("albedo", 0.5f);
            m_bsdf.reset(static_cast<BSDF *>(NoriObjectFactory::createInstance("diffuse", l)));
            m_bsdf->activate();
        }
    }

    Frame getFrame(const Intersection &its) const {
        Color3f grad[2];
        m_texture->evalGradient(its.uv, grad);

        float dDispDu = grad[0].getLuminance()*m_enhance;
        float dDispDv = grad[1].getLuminance()*m_enhance;

        /* Build a perturbed frame -- ignores the usually
           negligible normal derivative term */
        Vector3f dpdu = its.dpdu + its.shFrame.n * (
                dDispDu - its.shFrame.n.dot( its.dpdu));
        Vector3f dpdv = its.dpdv + its.shFrame.n * (
                dDispDv - its.shFrame.n.dot( its.dpdv));

        Frame result;
        result.n = dpdu.cross(dpdv).normalized();
        result.s = (dpdu - result.n * result.n.dot(dpdu)).normalized();
        result.t = result.n.cross(result.s);

        if (result.n.dot(its.geoFrame.n) < 0)
            result.n *= -1;
        return result;
    }

    virtual Color3f eval(const BSDFQueryRecord &bRec) const override {
        // return m_bsdf->eval(bRec);
        const Intersection& its = bRec.its;
        Intersection perturbed(its);
        perturbed.shFrame = getFrame(its);

        BSDFQueryRecord perturbedQuery(perturbed.toLocal(its.toWorld(bRec.wi)),
                                       perturbed.toLocal(its.toWorld(bRec.wo)),
                                       bRec.measure);
        perturbedQuery.its = perturbed;

        if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
            return Color3f(0.f);
        perturbedQuery.eta = bRec.eta;
        perturbedQuery.uv = bRec.uv;
        perturbedQuery.p = bRec.p;
        return m_bsdf->eval(perturbedQuery);
        

    }


    virtual float pdf(const BSDFQueryRecord &bRec) const override {
        // return m_bsdf->pdf(bRec);

        const Intersection& its = bRec.its;
        Intersection perturbed(its);
        perturbed.shFrame = getFrame(its);

        BSDFQueryRecord perturbedQuery(perturbed.toLocal(its.toWorld(bRec.wi)),
                                perturbed.toLocal(its.toWorld(bRec.wo)),
                                bRec.measure);
        perturbedQuery.its = perturbed;
        if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
            return 0;
        perturbedQuery.eta = bRec.eta;
        perturbedQuery.uv = bRec.uv;
        perturbedQuery.p = bRec.p;
        return m_bsdf->pdf(perturbedQuery);
        

    }

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample2d) const override {
        // return m_bsdf->sample(bRec, sample2d);
        
        const Intersection& its = bRec.its;
        Intersection perturbed(its);
        perturbed.shFrame = getFrame(its);

        BSDFQueryRecord perturbedQuery(perturbed.toLocal(its.toWorld(bRec.wi)));
        perturbedQuery.its = perturbed;
        perturbedQuery.eta = bRec.eta;
        perturbedQuery.uv = bRec.uv;
        perturbedQuery.p = bRec.p;

        Color3f result = m_bsdf->sample(perturbedQuery, sample2d);
        if (!result.isZero()) {
            bRec.p = perturbedQuery.p;
            bRec.uv = perturbedQuery.uv;
            bRec.eta = perturbedQuery.eta;
            bRec.wo = its.toLocal(perturbed.toWorld(perturbedQuery.wo));
            if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(perturbedQuery.wo) <= 0)
                return Color3f(0.0f);
        }

        return result;
    }


    virtual std::string toString() const override {
        std::string ret = "Bumpmap[];\n";
        return ret + m_bsdf->toString();
    }




private:
    float m_enhance;
    std::shared_ptr<BSDF> m_bsdf;
    std::shared_ptr<Texture<Color3f>> m_texture;
};

NORI_REGISTER_CLASS(Bumpmap, "bumpmap");
NORI_NAMESPACE_END
