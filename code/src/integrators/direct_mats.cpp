#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectMATSIntegrator : public Integrator {
private:
    float ray_length;
public:
    DirectMATSIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its_surface;
        if (!scene->rayIntersect(ray, its_surface)) return Color3f(0.0f);

        Color3f bsdf_cos_theta_over_pdf, Le = 0, Li;

        // L_e
        if(its_surface.mesh->isEmitter()){
            EmitterQueryRecord emitter_lRec(ray.o, its_surface.p, its_surface.shFrame.n); // intersection (camera), light, normal
            Le = its_surface.mesh->getEmitter()->eval(emitter_lRec);
        }
        // (1/n)*sum_n(bsdf * Li/pdf * cos theta)
        // Notice here we give wo to sample wi (reciprocity)
        BSDFQueryRecord bRec(its_surface.shFrame.toLocal(-ray.d));
        bRec.uv = its_surface.uv;
        bRec.p = its_surface.p;
        bsdf_cos_theta_over_pdf = its_surface.mesh->getBSDF()->sample(bRec, sampler->next2D());

        Ray3f shadowRay(its_surface.p, its_surface.shFrame.toWorld(bRec.wo));
        
        Intersection its_emitter;
        if(!scene->rayIntersect(shadowRay, its_emitter)||!its_emitter.mesh->isEmitter()) return Le;
        EmitterQueryRecord lRec(its_surface.p, its_emitter.p, its_emitter.shFrame.n);
        Li = its_emitter.mesh->getEmitter()->eval(lRec);

        return Le + Li * bsdf_cos_theta_over_pdf;
    }
    
    std::string toString() const {
        return "DirectMATSIntegrator[]";
    }
};

NORI_REGISTER_CLASS(DirectMATSIntegrator, "direct_mats");
NORI_NAMESPACE_END