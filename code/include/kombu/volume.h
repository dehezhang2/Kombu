#if !defined(__NORI_VOLUME_H)
#define __NORI_VOLUME_H

#include <nori/object.h>
#include <nori/bbox.h>


NORI_NAMESPACE_BEGIN

/**
 * \brief Superclass of all volumes
 */

class Volume: public NoriObject{
public:
  

    inline const BoundingBox3f* getBoundingBox() const{
        return m_bbox;
    }

    void setBoundingBox(const BoundingBox3f& bbox) {
        m_bbox = &bbox;
    }

    virtual float lookup(const Point3f &p) const = 0;

    virtual float getStepSize() const = 0;

    virtual float getMaximumFloatValue() const = 0;
    
    virtual void addChild(NoriObject *obj) {}

    virtual std::string toString() const {}

    EClassType getClassType() const override { return EVolume; }

protected:
    const BoundingBox3f* m_bbox;
    
};

NORI_NAMESPACE_END
#endif