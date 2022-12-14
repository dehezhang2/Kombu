#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN
using namespace std;

class VolPathMISIntegrator : public Integrator {
protected:
    Color3f recursiveTransmittence(const Scene* scene, const Point3f &p1, bool p1_onsurface, const Point3f &p2, bool p2_onsurface, const Medium* medium, Sampler* sampler) const{
        Vector3f d = p2 - p1;
        float remaining = d.norm();
        d /= remaining;

        float length_factor = p2_onsurface ? (1 - 1e-3f) : 1;
        Ray3f ray(p1, d, p1_onsurface ? Epsilon : 0, remaining * length_factor);

        Color3f transmittance(1.0f);
        Intersection its;

        while(remaining > 0){
            bool has_intersect = scene->rayIntersect(ray, its);
            if(has_intersect && !its.mesh->isMedium()) return 0.f;
            if(medium)
                transmittance *= medium->evalTransmittance(MediumQueryRecord(ray.o, ray.o + ray.d * std::min(its.t, remaining)), sampler);
            if(!has_intersect || transmittance.maxCoeff() == 0.f) 
                break;
            if(its.mesh->isMedium()){
                if(its.shFrame.n.dot(ray.d) >= 0) medium = nullptr;
                else medium = its.mesh->getMedium();
            }
            ray.o = ray(its.t);
            remaining -= its.t;
            ray.maxt = remaining * length_factor;
            ray.mint = Epsilon;
        }
        return transmittance;
    }

    void recursiveEmitterChecking(const Scene* scene, const Medium* medium, Ray3f ray, Intersection &_its, Color3f& value, float& pdf_em, Sampler* sampler) const {
        Intersection its2, *its = &_its;
        Color3f transmittance(1.f);
        bool has_intersect = false;
        Point3f origin = its->p;
        while(true){
            has_intersect = scene->rayIntersect(ray, *its);
            if(medium)
                transmittance *= medium->evalTransmittance(MediumQueryRecord(ray.o, ray.o + ray.d * its->t), sampler);
            if(has_intersect &&( !its->mesh->isMedium() || its->mesh->isEmitter()) )
                break;
            if(!has_intersect)
                break;
            if(transmittance.maxCoeff()==0.f) 
                return;
            // if (its->mesh->isMedium())
            //     medium = its->mesh->getMedium();
            if(its->mesh->isMedium()){
                break;
                if(its->shFrame.n.dot(ray.d) >= 0) medium = nullptr;
                else medium = its->mesh->getMedium();
            }
            ray.o = ray(its->t);
            ray.mint = Epsilon;
            its = &its2;
        }
        if(has_intersect){
            if(its->mesh->isEmitter()){
                EmitterQueryRecord emitter_lRec(origin, its->p, its->shFrame.n);    
                value = transmittance * its->mesh->getEmitter()->eval(emitter_lRec);
                pdf_em = its->mesh->getEmitter()->pdf(emitter_lRec) / scene->getLights().size();
            }
        } 
    }
    
    bool m_inmedium;
public:
    VolPathMISIntegrator(const PropertyList &props) {
        m_inmedium = props.getBoolean("inmedium", false);
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
        bool has_intersection;
        bool mats_strategy = true;
        if(m_inmedium){
            auto mediums = scene->getMediums();
            for(auto medium : mediums){
                Point3f ray_o = ray.o;
                if(medium->contains(ray_o)){
                    current_medium = medium;
                    scene->rayIntersect(ray, its_surface);
                    // cout << "haha\n";
                    if(its_surface.mesh->isEmitter()){
                        EmitterQueryRecord emitter_lRec(ray.o, its_surface.p, its_surface.shFrame.n); 
                        Li += its_surface.mesh->getEmitter()->eval(emitter_lRec) * throughput * current_medium->evalTransmittance(MediumQueryRecord(ray.o, its_surface.p), sampler);
                    }
                }
            }
        }

        while(true){
            MediumQueryRecord mRec(incident_ray.o, -incident_ray.d, its_surface.t);
            if(current_medium && current_medium->sample_intersection(mRec, sampler)){
                // volumetric rendering
                throughput *= mRec.albedo;
                // Ruassian Roulette 
                if(++bounce_cnt > 3){
                    float success_prob = std::min(throughput.getLuminance(), 0.99f);
                    if(sampler->next1D() > success_prob) break;
                    throughput /= success_prob;
                }
                // sample light
                const Emitter* light = scene->getRandomEmitter(sampler->next1D());
                EmitterQueryRecord lRec(mRec.p);
                if(light->isInfiniteDistance()){
                    lRec.bSphere_center = scene->getBoundingBox().getCenter();
                    lRec.bSphere_radius = (lRec.bSphere_center - scene->getBoundingBox().max).norm();
                }
                
                Color3f Li_over_pdf = light -> sample(lRec, sampler->next2D()) * scene->getLights().size();
                Color3f transmittance = recursiveTransmittence(scene, lRec.shadowRay.o, false, lRec.p, light->onSurface(), current_medium, sampler);
                float pdf_em_em = light -> pdf(lRec) / scene->getLights().size();
                if(transmittance.maxCoeff()!=0){
                    PhaseFunctionQueryRecord pRec(-incident_ray.d, lRec.wi, ESolidAngle);
                    float pdf_mat_em = current_medium->getPhaseFunction()->eval(pRec);
                    w_ems = (pdf_em_em + pdf_mat_em > 0 ? pdf_em_em/(pdf_em_em + pdf_mat_em) : 0.f);
                    if(lRec.isDelta) w_ems = 1.f;
                    Li += throughput * w_ems * transmittance * pdf_mat_em * Li_over_pdf ;
                }
                prev_discrete = false;
                
                // sample phasefunction ray
                pdf_mat_mat = current_medium->sample_phase(mRec, sampler->next2D());
                incident_ray = Ray3f(mRec.p, mRec.wo);
                has_intersection = scene->rayIntersect(incident_ray, its_surface);
            } else {
                if(current_medium) throughput *= mRec.albedo;
                // same as path_mits
                // count Li
                Color3f Le_tr(0.f);
                float pdf_em_mat = 0.f;
                recursiveEmitterChecking(scene, current_medium, incident_ray, its_surface, Le_tr, pdf_em_mat, sampler);
                has_intersection = its_surface.valid;

                if(!has_intersection){
                    if (scene->getEnvLight() == nullptr) {
                        break;
                    } else {
                        EmitterQueryRecord lRec(incident_ray.o);
                        lRec.wi = incident_ray.d;
                        if(!prev_discrete){
                            float pdf_em_mat = scene->getEnvLight()->pdf(lRec) / scene->getLights().size();
                            w_mat = (pdf_em_mat + pdf_mat_mat > 0 ? (pdf_mat_mat / (pdf_em_mat + pdf_mat_mat)) : 0.f);
                        }
                        return Li + w_mat * throughput * scene->getEnvLight()->eval(lRec);
                    }
                }
                if(Le_tr.maxCoeff()!=0){
                    if(!prev_discrete){
                        w_mat = (pdf_em_mat + pdf_mat_mat > 0 ? (pdf_mat_mat / (pdf_em_mat + pdf_mat_mat)) : 0.f);
                    }
                    Li += throughput * Le_tr * w_mat;
                }
                // sample BSDF ray
                BSDFQueryRecord bRec(its_surface.shFrame.toLocal(-incident_ray.d));
                bRec.uv = its_surface.uv;
                bRec.p = its_surface.p;
                Color3f bsdf_cos_theta_over_pdf = its_surface.mesh->getBSDF()->sample(bRec, sampler->next2D());
                pdf_mat_mat = its_surface.mesh->getBSDF()->pdf(bRec);
                prev_discrete = (bRec.measure == EDiscrete);
                w_mat = (prev_discrete ? 1.f : w_mat);

                if(bRec.isEmissive){
                    Li += throughput * bRec.radiance;
                }
                
                // Ruassian Roulette 
                if(++bounce_cnt > 5){
                    float success_prob = std::min(throughput.getLuminance(), 0.99f);
                    if(sampler->next1D() > success_prob) break;
                    throughput /= success_prob;
                }
                

                // sample light source
                // Note that this MIS corresponds to the BSDF sample in next iteration
                if(!prev_discrete){
                    
                    if(its_surface.mesh->isMedium()) current_medium = its_surface.mesh->getMedium();
                    const Emitter* light = scene->getRandomEmitter(sampler->next1D());
                    EmitterQueryRecord lRec(its_surface.p);
                    if(light->isInfiniteDistance()){
                        lRec.bSphere_center = scene->getBoundingBox().getCenter();
                        lRec.bSphere_radius = (lRec.bSphere_center - scene->getBoundingBox().max).norm();
                    }
                    Color3f Li_over_pdf = light -> sample(lRec, sampler->next2D()) * scene->getLights().size();
                    const Medium* medium = its_surface.mesh->isMedium() ?its_surface.mesh->getMedium(): current_medium;
                    Color3f transmittance = recursiveTransmittence(scene, lRec.shadowRay.o, true, lRec.p, light->onSurface(), medium, sampler);

                    float pdf_em_em = light -> pdf(lRec) / scene->getLights().size();
                    if(transmittance.maxCoeff()!=0){
                        float cos_theta_i = std::max(0.f, its_surface.shFrame.n.dot(lRec.wi));
                        BSDFQueryRecord bRec_em(its_surface.shFrame.toLocal(-incident_ray.d), its_surface.shFrame.toLocal(lRec.wi), ESolidAngle);
                        bRec_em.uv = its_surface.uv;
                        bRec_em.p = its_surface.p;
                        Color3f bsdf = its_surface.mesh->getBSDF()->eval(bRec_em);
                        float pdf_mat_em = its_surface.mesh->getBSDF()->pdf(bRec_em);
                        w_ems = (pdf_em_em + pdf_mat_em > 0 ? pdf_em_em/(pdf_em_em + pdf_mat_em) : 0.f);
                        if(lRec.isDelta) w_ems = 1.f;
                        Li += throughput * w_ems * bsdf * Li_over_pdf * cos_theta_i * transmittance;
                    }
                }
                // update throughput and incident_ray
                incident_ray = Ray3f(its_surface.p, its_surface.shFrame.toWorld(bRec.wo));
                throughput *= bsdf_cos_theta_over_pdf;
                
                if(its_surface.mesh->isMedium()){
                    if(its_surface.shFrame.n.dot(incident_ray.d) >= 0) current_medium = nullptr;
                    else current_medium = its_surface.mesh->getMedium();
                    has_intersection = scene->rayIntersect(incident_ray, its_surface);
                }
                
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