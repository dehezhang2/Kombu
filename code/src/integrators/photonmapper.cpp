/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/photon.h>

NORI_NAMESPACE_BEGIN

class PhotonMapper : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    PhotonMapper(const PropertyList &props) {
        /* Lookup parameters */
        m_photonCount  = props.getInteger("photonCount", 1000000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);
    }

    virtual void preprocess(const Scene *scene) override {
        cout << "Gathering " << m_photonCount << " photons .. ";
        cout.flush();

        /* Create a sample generator for the preprocess step */
        Sampler *sampler = static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList()));

        /* Allocate memory for the photon map */
        m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
        m_photonMap->reserve(m_photonCount);

		/* Estimate a default photon radius */
		if (m_photonRadius == 0)
			m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;

	

		/* How to add a photon?
		 * m_photonMap->push_back(Photon(
		 *	Point3f(0, 0, 0),  // Position
		 *	Vector3f(0, 0, 1), // Direction
		 *	Color3f(1, 2, 3)   // Power
		 * ));
		 */

		// put your code to trace photons here
        m_emittedCount = 0;
        int registeredCount;
        for(m_emittedCount = 0, registeredCount = 0; registeredCount < m_photonCount; m_emittedCount++){
            Ray3f incident_ray;
            Intersection its_surface;
            auto light = scene->getRandomEmitter(sampler->next1D());
            Color3f W = light->samplePhoton(incident_ray, sampler->next2D(), sampler->next2D()) * scene->getLights().size();
            int bounce_cnt = 0;
            while(true){
                Intersection its_surface;
                if(!scene->rayIntersect(incident_ray, its_surface)) break;
                // Store diffuse photons
                if(its_surface.mesh->getBSDF()->isDiffuse()){
                    m_photonMap->push_back(Photon(
                        its_surface.p,
                        -incident_ray.d.normalized(),
                        W
                    ));
                    registeredCount++;
                }
                // Ruassian Roulette 
                if(++bounce_cnt > 3){
                    float success_prob = std::min(W.getLuminance(), 0.99f);
                    if(sampler->next1D() > success_prob) break;
                    W /= success_prob;
                }
                // sample BSDF ray
                BSDFQueryRecord bRec(its_surface.shFrame.toLocal(-incident_ray.d));
                bRec.uv = its_surface.uv;
                bRec.p = its_surface.p;
                Color3f bsdf_cos_theta_over_pdf = its_surface.mesh->getBSDF()->sample(bRec, sampler->next2D());
                incident_ray = Ray3f(its_surface.p, its_surface.shFrame.toWorld(bRec.wo));
                // update throughput
                W *= bsdf_cos_theta_over_pdf;
            }
        }
		/* Build the photon map */
        m_photonMap->build();
    }

    virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const override {
    	
		/* How to find photons?
		 * std::vector<uint32_t> results;
		 * m_photonMap->search(Point3f(0, 0, 0), // lookup position
		 *                     m_photonRadius,   // search radius
		 *                     results);
		 *
		 * for (uint32_t i : results) {
		 *    const Photon &photon = (*m_photonMap)[i];
		 *    cout << "Found photon!" << endl;
		 *    cout << " Position  : " << photon.getPosition().toString() << endl;
		 *    cout << " Power     : " << photon.getPower().toString() << endl;
		 *    cout << " Direction : " << photon.getDirection().toString() << endl;
		 * }
		 */
		// put your code for path tracing with photon gathering here
        Color3f Li = 0;
        Color3f throughput = 1.f;
        Ray3f incident_ray = _ray;
        int bounce_cnt = 0;
        while(true){
            Intersection its_surface;
            if (!scene->rayIntersect(incident_ray, its_surface)) break;
            // count Li
            if(its_surface.mesh->isEmitter()){
                EmitterQueryRecord emitter_lRec(incident_ray.o, its_surface.p, its_surface.shFrame.n); // intersection (camera), light, normal
                Li += its_surface.mesh->getEmitter()->eval(emitter_lRec) * throughput;
            }
            if(its_surface.mesh->getBSDF()->isDiffuse()){
                std::vector<uint32_t> results;
		        m_photonMap->search(its_surface.p, // lookup position
		                            m_photonRadius,   // search radius
		                            results);
                
                Color3f diffuse_Li = 0.f;
		        for (uint32_t i : results) {
		            const Photon &photon = (*m_photonMap)[i];
                    BSDFQueryRecord bRec(its_surface.shFrame.toLocal(-incident_ray.d), its_surface.shFrame.toLocal(photon.getDirection()), ESolidAngle);
                    bRec.uv = its_surface.uv;
                    bRec.p = its_surface.p;
                    Color3f bsdf = its_surface.mesh->getBSDF()->eval(bRec);
                    diffuse_Li += (photon.getPower()) * bsdf;
		        }
                float area = m_photonRadius * m_photonRadius * M_PI;
                diffuse_Li /= area * m_emittedCount;
                Li += throughput * diffuse_Li;
                break;
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
            incident_ray = Ray3f(its_surface.p, its_surface.shFrame.toWorld(bRec.wo));
            // update throughput
            throughput *= bsdf_cos_theta_over_pdf;
        }
        return Li;
    }

    virtual std::string toString() const override {
        return tfm::format(
            "PhotonMapper[\n"
            "  photonCount = %i,\n"
            "  photonRadius = %f\n"
            "]",
            m_photonCount,
            m_photonRadius
        );
    }
private:
    /* 
     * Important: m_photonCount is the total number of photons deposited in the photon map,
     * NOT the number of emitted photons. You will need to keep track of those yourself.
     */ 
    int m_photonCount;
    int m_emittedCount;
    float m_photonRadius;
    std::unique_ptr<PhotonMap> m_photonMap;
};

NORI_REGISTER_CLASS(PhotonMapper, "photonmapper");
NORI_NAMESPACE_END
