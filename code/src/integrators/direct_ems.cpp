#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectEMSIntegrator : public Integrator {
private:
    float ray_length;
public:
    DirectEMSIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        //https://www.pbr-book.org/3ed-2018/Light_Transport_I_Surface_Reflection/Direct_Lighting
        Intersection its;
        if (!scene->rayIntersect(ray, its)) return Color3f(0.0f);
        auto lights = scene -> getLights();
        
        Color3f Li_over_pdf, bsdf, Lo = 0;
        float cos_theta;

        // L_e
        if(its.mesh->isEmitter()){
            EmitterQueryRecord emitter_lRec(ray.o, its.p, its.shFrame.n); // intersection (camera), light, normal
            Lo = its.mesh->getEmitter()->eval(emitter_lRec);
        }
        // (1/n)*sum_n(bsdf * Li/pdf * cos theta)
        for(auto light : lights){
            EmitterQueryRecord lRec(its.p);
            Li_over_pdf = light -> sample(lRec, sampler->next2D());

            if(scene -> rayIntersect(lRec.shadowRay)) continue;
            cos_theta = its.shFrame.n.dot(lRec.wi);
            if(cos_theta < 0) cos_theta = -cos_theta;
            BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d), its.shFrame.toLocal(lRec.wi), ESolidAngle);
            bRec.uv = its.uv;
            bRec.p = its.p;
            bsdf = its.mesh->getBSDF()->eval(bRec);
            Lo += bsdf * Li_over_pdf * cos_theta;
        }
        return Lo;
    }
    
    std::string toString() const {
        return "DirectEMSIntegrator[]";
    }
};

NORI_REGISTER_CLASS(DirectEMSIntegrator, "direct_ems");
NORI_NAMESPACE_END