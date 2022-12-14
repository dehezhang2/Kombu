#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
NORI_NAMESPACE_BEGIN

class DirectIntegrator : public Integrator {
private:
    float ray_length;
public:
    DirectIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its)){
            if (scene->getEnvLight() == nullptr) {
                return 0.f;
            } else {
                EmitterQueryRecord lRec;
                lRec.wi = ray.d;
                return scene->getEnvLight()->eval(lRec);
            }
        }
        auto lights = scene -> getLights();
        Vector2f sample;
        Color3f Ei, bsdf, Lo = 0;
        float cos_theta;
        for(auto light : lights){
            EmitterQueryRecord lRec(its.p);
            if(light->isInfiniteDistance()){
                    lRec.bSphere_center = scene->getBoundingBox().getCenter();
                    lRec.bSphere_radius = (lRec.bSphere_center - scene->getBoundingBox().max).norm();
            }
            Ei = light->sample(lRec, sample); // integral of Li d omega

            if(scene -> rayIntersect(lRec.shadowRay)) continue;
            cos_theta = its.shFrame.n.dot(lRec.wi);
            if(cos_theta < 0) cos_theta = -cos_theta;

            BSDFQueryRecord bRec(its.shFrame.toLocal(lRec.wi), its.shFrame.toLocal(-ray.d), ESolidAngle);
            bRec.uv = its.uv;
            bRec.p = its.p;
            bsdf = its.mesh->getBSDF()->eval(bRec);

            Lo += bsdf * Ei * cos_theta;
        }
        return Lo;
    }
    
    std::string toString() const {
        return "DirectIntegrator[]";
    }
};

NORI_REGISTER_CLASS(DirectIntegrator, "direct");
NORI_NAMESPACE_END