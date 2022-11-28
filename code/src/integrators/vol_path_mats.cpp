#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class VolPathMATSIntegrator : public Integrator {
private:
    float ray_length;
public:
    VolPathMATSIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Color3f Li = 0;
        // larger throughput => larger contribution => larger probability for russian
        Color3f throughput = 1.f;
        Ray3f incident_ray = ray;
        int bounce_cnt = 0;
        Medium* current_medium = nullptr;
        Intersection its_surface;
        bool use_volume = false;
        while(true){
            if (!scene->rayIntersect(incident_ray, its_surface)) break;
            // Check whether the current point is in the medium
            current_medium = nullptr;
            auto mediums = scene->getMediums();
            // cout << mediums.size() << endl;
            for(auto medium : mediums){
                if(medium->contains(incident_ray.o)){
                    current_medium = medium;
                    break; 
                }
            }
            use_volume = false;
            if(current_medium){
                // volumetric rendering
                Frame incident_local_frame(incident_ray.d);
                MediumQueryRecord mRec(its_surface.p, Vector3f(0.f), its_surface.t);
                Color3f albedo = current_medium->sample_intersection(mRec, sampler->next1D());
                if(mRec.hitMedium){
                    // volumetric
                    use_volume = true;
                    // sample phasefunction ray
                    current_medium->sample_phase(mRec, sampler->next2D());
                    incident_ray = Ray3f(mRec.p, incident_local_frame.toWorld(mRec.wo));
                    // update throughput 
                    throughput *= albedo;
                } 
            } 
            if(!use_volume){
                // same as path_mats
                // count Li
                if(its_surface.mesh->isEmitter()){
                    EmitterQueryRecord emitter_lRec(incident_ray.o, its_surface.p, its_surface.shFrame.n); // intersection (camera), light, normal
                    Li += its_surface.mesh->getEmitter()->eval(emitter_lRec) * throughput;
                }  
                // sample BSDF ray
                BSDFQueryRecord bRec(its_surface.shFrame.toLocal(-incident_ray.d));
                bRec.uv = its_surface.uv;
                bRec.p = its_surface.p;
                Color3f bsdf_cos_theta_over_pdf = its_surface.mesh->getBSDF()->sample(bRec, sampler->next2D());
                incident_ray = Ray3f(its_surface.p, its_surface.shFrame.toWorld(bRec.wo));
                // update throughput
                throughput *= bsdf_cos_theta_over_pdf;   
            }
            // Ruassian Roulette 
            if(++bounce_cnt > 3){
                float success_prob = std::min(throughput.getLuminance(), 0.99f);
                if(sampler->next1D() > success_prob) break;
                throughput /= success_prob;
            }            
        }
        return Li;
    }
    
    std::string toString() const {
        return "VolPathMATSIntegrator[]";
    }
};

NORI_REGISTER_CLASS(VolPathMATSIntegrator, "vol_path_mats");
NORI_NAMESPACE_END