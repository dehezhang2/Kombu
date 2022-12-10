// #include <nori/integrator.h>
// #include <nori/scene.h>
// #include <nori/bsdf.h>
// #include <nori/sampler.h>

// NORI_NAMESPACE_BEGIN
// using namespace std;

// class PathMISIntegrator : public Integrator {
// private:
//     float ray_length;
// public:
//     PathMISIntegrator(const PropertyList &props) {
//     }

//     Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
//         Color3f Li = 0;
//         // larger throughput => larger contribution => larger probability for russian
//         Color3f throughput = 1.f;
//         Ray3f incident_ray = ray;
//         int bounce_cnt = 0;
//         // The first ray just compute Le, therefore, no MIS => set w_mat to one
//         float w_mat = 1.f, w_ems = 0.f, pdf_mat_mat;
//         bool prev_discrete = true;
//         while(true){
//             Intersection its_surface;
//             if (!scene->rayIntersect(incident_ray, its_surface)) break;
//             // count Li
//             if(its_surface.mesh->isEmitter()){
//                 EmitterQueryRecord emitter_lRec(incident_ray.o, its_surface.p, its_surface.shFrame.n); // intersection (camera), light, normal
//                 if(!prev_discrete){
//                     float pdf_em_mat = its_surface.mesh->getEmitter()->pdf(emitter_lRec) / scene->getLights().size();
//                     w_mat = (pdf_em_mat + pdf_mat_mat > 0 ? (pdf_mat_mat / (pdf_em_mat + pdf_mat_mat)) : 0.f);
//                 }
//                 Li += its_surface.mesh->getEmitter()->eval(emitter_lRec) * throughput * w_mat;
//             }
//             // Ruassian Roulette 
//             if(++bounce_cnt > 3){
//                 float success_prob = std::min(throughput.getLuminance(), 0.99f);
//                 if(sampler->next1D() > success_prob) break;
//                 throughput /= success_prob;
//             }
//             // sample BSDF ray
//             BSDFQueryRecord bRec(its_surface.shFrame.toLocal(-incident_ray.d));
//             bRec.uv = its_surface.uv;
//             bRec.p = its_surface.p;
//             Color3f bsdf_cos_theta_over_pdf = its_surface.mesh->getBSDF()->sample(bRec, sampler->next2D());
//             pdf_mat_mat = its_surface.mesh->getBSDF()->pdf(bRec);
//             prev_discrete = (bRec.measure == EDiscrete);
//             w_mat = (prev_discrete ? 1.f : w_mat);
//             // sample light source
//             // Note that this MIS corresponds to the BSDF sample in next iteration
//             if(!prev_discrete){
//                 const Emitter* light = scene->getRandomEmitter(sampler->next1D());
//                 EmitterQueryRecord lRec(its_surface.p);
//                 Color3f Li_over_pdf = light -> sample(lRec, sampler->next2D()) * scene->getLights().size();
//                 float pdf_em_em = light -> pdf(lRec) / scene->getLights().size();
//                 if(!scene -> rayIntersect(lRec.shadowRay)){
//                     float cos_theta_i = max(0.f, its_surface.shFrame.n.dot(lRec.wi));
//                     BSDFQueryRecord bRec_em(its_surface.shFrame.toLocal(-incident_ray.d), its_surface.shFrame.toLocal(lRec.wi), ESolidAngle);
//                     bRec_em.uv = its_surface.uv;
//                     bRec_em.p = its_surface.p;
//                     Color3f bsdf = its_surface.mesh->getBSDF()->eval(bRec_em);
//                     float pdf_mat_em = its_surface.mesh->getBSDF()->pdf(bRec_em);
//                     w_ems = (pdf_em_em + pdf_mat_em > 0 ? pdf_em_em/(pdf_em_em + pdf_mat_em) : 0.f);
//                     Li += throughput * w_ems * bsdf * Li_over_pdf * cos_theta_i;
//                 }
//             }
//             // update throughput and incident_ray
//             incident_ray = Ray3f(its_surface.p, its_surface.shFrame.toWorld(bRec.wo));
//             throughput *= bsdf_cos_theta_over_pdf;
//         }
//         return Li;
//     }
    
//     std::string toString() const {
//         return "PathMISIntegrator[]";
//     }
// };

// NORI_REGISTER_CLASS(PathMISIntegrator, "path_mis");
// NORI_NAMESPACE_END



#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>



NORI_NAMESPACE_BEGIN

class PathMisIntegrator : public Integrator {
public:
	PathMisIntegrator(const PropertyList &props) {
		/* No parameters this time */
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		// Initial radiance and throughput
		Color3f Li = 0, t = 1;
		Ray3f rayR = ray;
		float prob = 1, w_mats = 1, w_ems = 1;
		Color3f f(1,1,1);

		while (true) {
			Intersection its;

			//if intersect
			if (!scene->rayIntersect(rayR, its)) {
				if (scene->getEnvLight() == nullptr) {
					return Li;
				} else {
					EmitterQueryRecord lRec;
					lRec.wi = rayR.d;
					return Li + w_mats * t * scene->getEnvLight()->eval(lRec);
				}
			}

			// Emitted
			Color3f Le = 0;
			if (its.mesh->isEmitter()) {
				EmitterQueryRecord lRecE(rayR.o, its.p, its.shFrame.n);
				Le = its.mesh->getEmitter()->eval(lRecE);
			}
			Li += t * w_mats * Le;

			//Russian roulette
			prob = std::min(t.maxCoeff(), .99f);
			if (sampler->next1D() >= prob)
				return Li;

			t /= prob;

			//emiter sampling
			Color3f L_ems = 0;

			// get a random emitter
			const Emitter * emitter = scene->getRandomEmitter(sampler->next1D());

			// reflected
			EmitterQueryRecord lRec_ems;
			lRec_ems.ref = its.p;
			Color3f Li_ems = emitter->sample(lRec_ems, sampler->next2D())*scene->getLights().size();
			float pdf_ems = emitter->pdf(lRec_ems);

			// BSDF 
			BSDFQueryRecord bRec_ems(its.shFrame.toLocal(-rayR.d), its.shFrame.toLocal(lRec_ems.wi), ESolidAngle);
			bRec_ems.uv = its.uv;
			Color3f f_ems = its.mesh->getBSDF()->eval(bRec_ems);
			if (pdf_ems + its.mesh->getBSDF()->pdf(bRec_ems) != 0)
				w_ems = pdf_ems / (pdf_ems + its.mesh->getBSDF()->pdf(bRec_ems));

			// check if shadow ray is occluded
			Intersection its_ems;
			if (!scene->rayIntersect(lRec_ems.shadowRay, its_ems))
				L_ems = f_ems * Li_ems * std::max(0.f, Frame::cosTheta(its.shFrame.toLocal(lRec_ems.wi)));

			Li += t * w_ems * L_ems;

			//BSDF sampling
			BSDFQueryRecord bRec(its.shFrame.toLocal(-rayR.d));
			Color3f f = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
			t *= f;
			// shoot next ray
			rayR = Ray3f(its.p, its.toWorld(bRec.wo));

			//next Le
			float pdf_mats = its.mesh->getBSDF()->pdf(bRec);

			Intersection itsR;
			if (scene->rayIntersect(rayR, itsR)) {
				if (itsR.mesh->isEmitter()) {
					EmitterQueryRecord lRec_mats = EmitterQueryRecord(its.p, itsR.p, itsR.shFrame.n);
					if (pdf_mats + itsR.mesh->getEmitter()->pdf(lRec_mats) != 0)
						w_mats = pdf_mats / (pdf_mats + itsR.mesh->getEmitter()->pdf(lRec_mats));
				}
			} else if (scene->getEnvLight() != nullptr) {
				EmitterQueryRecord lRec_mats;
				lRec_mats.wi = rayR.d;
				if (pdf_mats + scene->getEnvLight()->pdf(lRec_mats) != 0)
					w_mats = pdf_mats / (pdf_mats + scene->getEnvLight()->pdf(lRec_mats));
			}

			if (bRec.measure == EDiscrete) {
				w_mats = 1;
				w_ems = 0;
			}
		}
	} 

	std::string toString() const {
		return "PathMisIntegrator[]";
	}
};

NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");
NORI_NAMESPACE_END
