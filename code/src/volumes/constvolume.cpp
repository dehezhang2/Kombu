#include <kombu/volume.h>

NORI_NAMESPACE_BEGIN
class ConstantVolume : public Volume{
public:
    ConstantVolume(const PropertyList &props){
        m_value = props.getFloat("value", 0.f);
    }

    float lookup(Shape* shape, const Point3f &_p) {
        return m_value;
    }

    float getStepSize(Shape* shape){
        return std::numeric_limits<float>::infinity();
    }

    float getMaximumValue() const {
        return m_value;
    }

    std::string toString() const {
        return tfm::format(
                    "ConstantVolume[\n"
                    "  value = %s,\n"
                    "]",
                    m_value);
    }

protected:
    float m_value;
};
NORI_REGISTER_CLASS(ConstantVolume, "constvolume");
NORI_NAMESPACE_END