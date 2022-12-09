#include <nori/emitter.h>
NORI_NAMESPACE_BEGIN

class PointLight: public Emitter {
    public:
        PointLight(const PropertyList &propList) {
            /* No parameters this time */
            m_position = propList.getPoint3("position", Point3f());
            m_power = propList.getColor("power", Color3f());
        }
        Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample) const {
            // light source position
            lRec.p = m_position;
            // direction from light source
            lRec.wi = (lRec.p - lRec.ref).normalized();
            // vector to light source
            lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, (lRec.p - lRec.ref).norm() - Epsilon);
            lRec.pdf = 1.0;
            lRec.isDelta = true;
            return eval(lRec) / pdf(lRec);
        }
        Color3f eval(const EmitterQueryRecord &lRec) const {
            return m_power / (4.0 * M_PI * (lRec.p - lRec.ref).squaredNorm());
        }
        float pdf(const EmitterQueryRecord &lRec) const {
            return 1.0;
        }
        virtual std::string toString() const {
		    return "PointLight[]";
	    }
    protected:
        Point3f m_position;
        Color3f m_power;
};

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END