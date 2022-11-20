#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>
#include <nori/ray.h>

NORI_NAMESPACE_BEGIN

class AverageVisibility : public Integrator {
private:
    float ray_length;
public:
    AverageVisibility(const PropertyList &props) {
        ray_length = props.getFloat("length",1);
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its)) return Color3f(1.0f);
        Vector3f direction = Warp::sampleUniformHemisphere(sampler, its.shFrame.n);
        Ray3f sample_ray(its.p, direction, Epsilon, ray_length - Epsilon);
        if(scene->rayIntersect(sample_ray, its)) return Color3f(0.0f);
        return Color3f(1.0f);
    }

    std::string toString() const {
        return "AverageVisibility[]";
    }
};

NORI_REGISTER_CLASS(AverageVisibility, "av");
NORI_NAMESPACE_END