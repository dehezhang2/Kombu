#include <kombu/reference.h>
#include <nori/bbox.h>

NORI_NAMESPACE_BEGIN

void Reference::addChild(NoriObject *obj){
    switch (obj->getClassType()) {
        case EMesh:
            if (m_shape)
                throw NoriException(
                    "Reference: tried to register multiple shapes for one reference!");
            m_shape = static_cast<Shape *>(obj);
            if(m_shape->isInstance())
                throw NoriException(
                    "Reference: child cannot be instance!");
            else if(m_shape->isInstance())
                throw NoriException(
                    "Reference: child cannot be reference!");
            break;
        default:
            throw NoriException("Shape::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}
NORI_REGISTER_CLASS(Reference, "ref");
NORI_NAMESPACE_END
