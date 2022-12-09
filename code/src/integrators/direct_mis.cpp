#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectMISIntegrator : public Integrator {
private:
    float ray_length;
public:
    DirectMISIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its_surface;
        if (!scene->rayIntersect(ray, its_surface)) return Color3f(0.0f);
        auto lights = scene -> getLights();

        Color3f bsdf_cos_theta_over_pdf, Li_over_pdf, bsdf, Le = 0, Li_mat, Lo_mat, Lo_em = 0;
        float pdf_mat, pdf_em;
        // L_e
        if(its_surface.mesh->isEmitter()){
            EmitterQueryRecord emitter_lRec(ray.o, its_surface.p, its_surface.shFrame.n); // intersection (camera), light, normal
            Le = its_surface.mesh->getEmitter()->eval(emitter_lRec);
        }
        // (1/n)*sum_n(bsdf * Li/pdf * cos theta)

        // Lo_MAT
        // Notice here we give wo to sample wi (reciprocity)
        {
            BSDFQueryRecord bRec(its_surface.shFrame.toLocal(-ray.d));
            bRec.uv = its_surface.uv;
            bRec.p = its_surface.p;
            bsdf_cos_theta_over_pdf = its_surface.mesh->getBSDF()->sample(bRec, sampler->next2D());
            pdf_mat = its_surface.mesh->getBSDF()->pdf(bRec);
    
            Ray3f shadowRay(its_surface.p, its_surface.shFrame.toWorld(bRec.wo));
            
            Intersection its_emitter;
            if(scene->rayIntersect(shadowRay, its_emitter)&&its_emitter.mesh->isEmitter()){
                EmitterQueryRecord lRec(its_surface.p, its_emitter.p, its_emitter.shFrame.n);
                Li_mat = its_emitter.mesh->getEmitter()->eval(lRec);
                pdf_em = its_emitter.mesh->getEmitter()->pdf(lRec);
            }

            float weight = (pdf_em + pdf_mat > 0 ? pdf_mat/(pdf_em + pdf_mat) : 0.f);
            Lo_mat = weight * Li_mat * bsdf_cos_theta_over_pdf;
        }
        // Lo_em
        for(auto light : lights){
            EmitterQueryRecord lRec(its_surface.p);
            if(light->isDirectional()){
                lRec.bSphere_center = scene->getBoundingBox().getCenter();
                lRec.bSphere_radius = (lRec.bSphere_center - scene->getBoundingBox().max).norm();
            }
            Li_over_pdf = light -> sample(lRec, sampler->next2D());
            pdf_em = light -> pdf(lRec);

            if(scene -> rayIntersect(lRec.shadowRay)) continue;
            float cos_theta_i = its_surface.shFrame.n.dot(lRec.wi);
            // float cos_theta_i = Frame::cosTheta(its_surface.shFrame.toLocal(lRec.wi));
            // if(cos_theta_i < 0) cos_theta_i = -cos_theta_i;

            // BSDFQueryRecord bRec(its_surface.shFrame.toLocal(lRec.wi), its_surface.shFrame.toLocal(-ray.d), ESolidAngle);
            BSDFQueryRecord bRec(its_surface.shFrame.toLocal(-ray.d), its_surface.shFrame.toLocal(lRec.wi), ESolidAngle);
            bRec.uv = its_surface.uv;
            bRec.p = its_surface.p;
            bsdf = its_surface.mesh->getBSDF()->eval(bRec);
            pdf_mat = its_surface.mesh->getBSDF()->pdf(bRec);

            float weight = (pdf_em + pdf_mat > 0 ? pdf_em/(pdf_em + pdf_mat) : 0.f);
            if(lRec.isDelta) weight = 1.f;
            Lo_em += weight * bsdf * Li_over_pdf * cos_theta_i;
        }
        return Le + Lo_mat + Lo_em;
    }
    
    std::string toString() const {
        return "DirectMISIntegrator[]";
    }
};

NORI_REGISTER_CLASS(DirectMISIntegrator, "direct_mis");
NORI_NAMESPACE_END