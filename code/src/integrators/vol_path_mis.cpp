#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class VolPathMISIntegrator : public Integrator {
private:
    float ray_length;
public:
    VolPathMISIntegrator(const PropertyList &props) {
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
        const Medium* current_medium = nullptr;
        Intersection its_surface;
        bool has_intersection = scene->rayIntersect(incident_ray, its_surface);
        while(true){
            MediumQueryRecord mRec(incident_ray.o, incident_ray.d, its_surface.t);
            if(current_medium && current_medium->sample_intersection(mRec, sampler->next1D())){
                // volumetric rendering
                throughput *= mRec.albedo;

                // sample light
                const Emitter* light = scene->getRandomEmitter(sampler->next1D());
                EmitterQueryRecord lRec(mRec.p);
                Color3f Li_over_pdf = light -> sample(lRec, sampler->next2D()) * scene->getLights().size();
                float pdf_em_em = light -> pdf(lRec) / scene->getLights().size();
                if(!scene -> rayIntersect(lRec.shadowRay)){
                    PhaseFunctionQueryRecord pRec(-incident_ray.d, lRec.wi, ESolidAngle);
                    float pdf_mat_em = current_medium->getPhaseFunction()->eval(pRec);
                    w_ems = (pdf_em_em + pdf_mat_em > 0 ? pdf_em_em/(pdf_em_em + pdf_mat_em) : 0.f);
                    Li += throughput * w_ems * current_medium->evalTransmittance(mRec) * pdf_mat_em * Li_over_pdf ;
                }
                prev_discrete = false;
                
                // sample phasefunction ray
                pdf_mat_mat = current_medium->sample_phase(mRec, sampler->next2D());
                Frame incident_local_frame(incident_ray.d);
                incident_ray = Ray3f(mRec.p, incident_local_frame.toWorld(mRec.wo));
                has_intersection = scene->rayIntersect(incident_ray, its_surface);
            } else {
                if(!has_intersection) break;
                // same as path_mits
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
                        Li += throughput * w_ems * bsdf * Li_over_pdf * cos_theta_i;
                    }
                }
                // update throughput and incident_ray
                incident_ray = Ray3f(its_surface.p, its_surface.shFrame.toWorld(bRec.wo));
                has_intersection = scene->rayIntersect(incident_ray, its_surface);
                current_medium = its_surface.mesh->getMedium();
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

NORI_REGISTER_CLASS(VolPathMISIntegrator, "vol_path_mis");
NORI_NAMESPACE_END