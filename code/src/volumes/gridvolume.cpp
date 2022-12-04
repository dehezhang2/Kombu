#include <kombu/volume.h>
#include <QFile>
#include <QDataStream>
#if defined(PLATFORM_LINUX) || defined(PLATFORM_MACOS)
#include <sys/mman.h>
#include <fcntl.h>
#endif

NORI_NAMESPACE_BEGIN
class GridVolume : public Volume{
public:
    GridVolume(const PropertyList &props){
        m_worldToLocal = props.getTransform("toWorld", Transform()).inverse();
                
    }

    float lookup(const Point3f &p) const {
        
    }

    float getStepSize() const{
        return std::numeric_limits<float>::infinity();
    }

    float getMaximumFloatValue() const {
        
    }

    std::string toString() const {
        // return tfm::format(
        //             "GridVolume[\n"
        //             "  value = %s,\n"
        //             "]",
        //             m_value);
    }


protected:
    Transform m_worldToLocal;
    QString m_filename;
	size_t m_fileSize;
	float *m_data;
};
NORI_REGISTER_CLASS(GridVolume, "gridvolume");
NORI_NAMESPACE_END