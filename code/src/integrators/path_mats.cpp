#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class PathMATSIntegrator : public Integrator {
private:
    float ray_length;
public:
    PathMATSIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Color3f Li = 0;
        // larger throughput => larger contribution => larger probability for russian
        Color3f throughput = 1.f;
        Ray3f incident_ray = ray;
        int bounce_cnt = 0;
        while(true){
            Intersection its_surface;
            if (!scene->rayIntersect(incident_ray, its_surface)) break;
            // count Li
            if(its_surface.mesh->isEmitter()){
                EmitterQueryRecord emitter_lRec(incident_ray.o, its_surface.p, its_surface.shFrame.n); // intersection (camera), light, normal
                Li += its_surface.mesh->getEmitter()->eval(emitter_lRec) * throughput;
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
            bRec.its = its_surface;
            Color3f bsdf_cos_theta_over_pdf = its_surface.mesh->getBSDF()->sample(bRec, sampler->next2D());
            incident_ray = Ray3f(its_surface.p, its_surface.shFrame.toWorld(bRec.wo));
            // update throughput
            throughput *= bsdf_cos_theta_over_pdf;
        }
        return Li;
    }
    
    std::string toString() const {
        return "PathMATSIntegrator[]";
    }
};

NORI_REGISTER_CLASS(PathMATSIntegrator, "path_mats");
NORI_NAMESPACE_END