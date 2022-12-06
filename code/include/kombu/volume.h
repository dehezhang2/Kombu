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
    virtual float lookup(Shape* shape, const Point3f &_p) = 0;

    virtual float getStepSize(Shape* shape) = 0;

    virtual float getMaximumValue() const = 0;
    
    virtual void addChild(NoriObject *obj) {}

    virtual std::string toString() const {}

    EClassType getClassType() const override { return EVolume; }
    
};

NORI_NAMESPACE_END
#endif