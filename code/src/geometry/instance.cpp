#include <kombu/instance.h>
#include <nori/bbox.h>
#include <nori/transform.h>

NORI_NAMESPACE_BEGIN

void Instance::activate(){
    m_worldToLocal = m_localToWorld.inverse();
}
void Instance::addChild(NoriObject *obj){
     throw NoriException("Instance::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
}

const BoundingBox3f &Instance::getBoundingBox() const {
    if(!m_reference){
        throw NoriException(
					"There is no reference attached to this instance!");
    }
    const BoundingBox3f local_bbox = m_reference->getBoundingBox();
    if(!local_bbox.isValid())
        return local_bbox;
    BoundingBox3f m_bbox(m_localToWorld * local_bbox.min, m_localToWorld * local_bbox.max);
    return m_bbox;
}

BoundingBox3f Instance::getBoundingBox(uint32_t index) const {
    if(!m_reference){
        throw NoriException(
					"There is no reference attached to this instance!");
    }
    
    BoundingBox3f local_bbox = m_reference->getBoundingBox(index);
    if(!local_bbox.isValid())
        return local_bbox;
    BoundingBox3f m_bbox(m_localToWorld * local_bbox.min, m_localToWorld * local_bbox.max);
    return m_bbox;
}

bool Instance::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const{
    if(!m_reference){
        throw NoriException(
					"There is no reference attached to this instance!");
    }
    Ray3f local_ray = m_worldToLocal * ray;
    return m_reference->rayIntersect(index, local_ray, u, v, t);
}

void Instance::setHitInformation(uint32_t index, const Ray3f &ray, Intersection & its) const {
    if(!m_reference){
        throw NoriException(
					"There is no reference attached to this instance!");
    }
    Ray3f local_ray = m_worldToLocal * ray;
    m_reference->setHitInformation(index, local_ray, its);

    its.p = m_localToWorld * its.p;
    Vector3f normalized_coor = (m_localToWorld * its.shFrame.n).normalized();
    its.geoFrame = Frame(normalized_coor);
	its.shFrame = Frame(normalized_coor);
}

Point3f Instance::getCentroid(uint32_t index) const {
    if(!m_reference){
        throw NoriException(
					"There is no reference attached to this instance!");
    }
    Point3f local_centroid = m_reference->getCentroid(index);
    return m_localToWorld * local_centroid;
}

void Instance::sampleSurface(ShapeQueryRecord & sRec, const Point2f & sample) const {
    if(!m_reference){
        throw NoriException(
					"There is no reference attached to this instance!");
    }
    m_reference->sampleSurface(sRec, sample);
    sRec.p = m_localToWorld * sRec.p;
    sRec.n = m_localToWorld * sRec.n;
}

float Instance::pdfSurface(const ShapeQueryRecord & sRec) const {
    if(!m_reference){
        throw NoriException(
					"There is no reference attached to this instance!");
    }
    ShapeQueryRecord local_sRec(m_worldToLocal*sRec.ref, m_worldToLocal*sRec.p);
    return m_reference->pdfSurface(local_sRec);
}

NORI_REGISTER_CLASS(Instance, "instance");
NORI_NAMESPACE_END
