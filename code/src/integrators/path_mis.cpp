#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class PathMISIntegrator : public Integrator {
private:
    float ray_length;
public:
    PathMISIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Color3f Li = 0;
        // larger throughput => larger contribution => larger probability for russian
        Color3f throughput = 1.f;
        Ray3f incident_ray = ray;
        int bounce_cnt = 0;
        // The first ray just compute Le, therefore, no MIS => set w_mat to one
        float w_mat = 1.f, w_ems = 0.f, pdf_mat_mat;
        bool prev_discrete = true;
        while(true){
            Intersection its_surface;
            if (!scene->rayIntersect(incident_ray, its_surface)) break;
            // count Li
            if(its_surface.mesh->isEmitter()){
                EmitterQueryRecord emitter_lRec(incident_ray.o, its_surface.p, its_surface.shFrame.n); // intersection (camera), light, normal
                if(!prev_discrete){
                    float pdf_em_mat = its_surface.mesh->getEmitter()->pdf(emitter_lRec) / scene->getLights().size();
                    w_mat = (pdf_em_mat + pdf_mat_mat > 0 ? (pdf_mat_mat / (pdf_em_mat + pdf_mat_mat)) : 0.f);
                }
                Li += its_surface.mesh->getEmitter()->eval(emitter_lRec) * throughput * w_mat;
            }
            // Ruassian Roulette 
            if(++bounce_cnt > 3){
                float success_prob = std::min(throughput.getLuminance(), 0.99f);
                if(sampler->next1D() > success_prob) break;
                throughput /= success_prob;
            }
            // sample BSDF ray
            BSDFQueryRecord bRec(its_surface.shFrame.toLocal(-incident_ray.d));
            bRec.uv = its_surface.uv;
            bRec.p = its_surface.p;
            Color3f bsdf_cos_theta_over_pdf = its_surface.mesh->getBSDF()->sample(bRec, sampler->next2D());
            pdf_mat_mat = its_surface.mesh->getBSDF()->pdf(bRec);
            prev_discrete = (bRec.measure == EDiscrete);
            w_mat = (prev_discrete ? 1.f : w_mat);
            // sample light source
            // Note that this MIS corresponds to the BSDF sample in next iteration
            if(!prev_discrete){
                const Emitter* light = scene->getRandomEmitter(sampler->next1D());
                EmitterQueryRecord lRec(its_surface.p);
                if(light->isDirectional()){
                    lRec.bSphere_center = scene->getBoundingBox().getCenter();
                    lRec.bSphere_radius = (lRec.bSphere_center - scene->getBoundingBox().max).norm();
                }
                Color3f Li_over_pdf = light -> sample(lRec, sampler->next2D()) * scene->getLights().size();
                float pdf_em_em = light -> pdf(lRec) / scene->getLights().size();
                if(!scene -> rayIntersect(lRec.shadowRay)){
                    float cos_theta_i = max(0.f, its_surface.shFrame.n.dot(lRec.wi));
                    BSDFQueryRecord bRec_em(its_surface.shFrame.toLocal(-incident_ray.d), its_surface.shFrame.toLocal(lRec.wi), ESolidAngle);
                    bRec_em.uv = its_surface.uv;
                    bRec_em.p = its_surface.p;
                    Color3f bsdf = its_surface.mesh->getBSDF()->eval(bRec_em);
                    float pdf_mat_em = its_surface.mesh->getBSDF()->pdf(bRec_em);
                    w_ems = (pdf_em_em + pdf_mat_em > 0 ? pdf_em_em/(pdf_em_em + pdf_mat_em) : 0.f);
                    if(lRec.isDelta) w_ems = 1.f;
                    Li += throughput * w_ems * bsdf * Li_over_pdf * cos_theta_i;
                }
            }
            // update throughput and incident_ray
            incident_ray = Ray3f(its_surface.p, its_surface.shFrame.toWorld(bRec.wo));
            throughput *= bsdf_cos_theta_over_pdf;
        }
        return Li;
    }
    
    std::string toString() const {
        return "PathMISIntegrator[]";
    }
};

NORI_REGISTER_CLASS(PathMISIntegrator, "path_mis");
NORI_NAMESPACE_END