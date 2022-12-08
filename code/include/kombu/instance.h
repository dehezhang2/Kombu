#if !defined(__NORI_INSTANCE_H)
#define __NORI_INSTANCE_H
#include <nori/object.h>
#include <nori/shape.h>
#include <nori/transform.h>
#include "reference.h"

NORI_NAMESPACE_BEGIN
class Instance : public Shape{
public:
    Instance(const PropertyList &propList) { 
        m_localToWorld = propList.getTransform("toWorld", Transform());
        m_ref_name = propList.getString("ref", "");
    }
    virtual void activate() override;
    
    virtual void addChild(NoriObject *obj) override;

    //// Return an axis-aligned bounding box of the entire mesh
    virtual const BoundingBox3f &getBoundingBox() const override;

    //// Return an axis-aligned bounding box containing the given triangle
    virtual BoundingBox3f getBoundingBox(uint32_t index) const override;

    //// Ray-Shape intersection test
    virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const override;

    /// Set the intersection information: hit point, shading frame, UVs, etc.
    virtual void setHitInformation(uint32_t index, const Ray3f &ray, Intersection & its) const override;

    //// Return the centroid of the given triangle
    virtual Point3f getCentroid(uint32_t index) const override;

    /**
     * \brief Sample a point on the surface (potentially using the point sRec.ref to importance sample)
     * This method should set sRec.p, sRec.n and sRec.pdf
     * Probability should be with respect to area
     * */
    virtual void sampleSurface(ShapeQueryRecord & sRec, const Point2f & sample) const override;
    /**
     * \brief Return the probability of sampling a point sRec.p by the sampleSurface() method (sRec.ref should be set before)
     * sRec.n and sRec.pdf are ignored
     * */
    virtual float pdfSurface(const ShapeQueryRecord & sRec) const override;
    
    // bool link(const std::vector<Reference *>& references);
    virtual void setRef(Reference * ref) { m_reference = ref; }
    virtual std::string getRefName() { return m_ref_name; }

    /// Is this mesh an area emitter?
    virtual bool isEmitter() const override{
        return m_reference->isEmitter();
    }

    /// Is this mesh a medium?
    virtual bool isMedium() const override{
        return m_reference->isMedium();
    }

    /// Return a pointer to an attached area emitter instance
    Emitter *getEmitter() override{
        return m_reference->getEmitter();
    }

    /// Return a pointer to an attached area emitter instance (const version)
    const Emitter *getEmitter() const override{
        return m_reference->getEmitter();
    }

    /// Return a pointer to an attached medium instance
    Medium *getMedium() override{
        return m_reference->getMedium();
    }

    /// Return a pointer to an attached medium instance (const version)
    const Medium *getMedium() const override{
        return m_reference->getMedium();
    }

    /// Return a pointer to the BSDF associated with this mesh
    const BSDF *getBSDF() const { return m_reference->getBSDF(); }


    /// Return the total number of primitives in this shape
    virtual uint32_t getPrimitiveCount() const override{
        return m_reference->getPrimitiveCount();
    }
    virtual bool isInstance() const override{return true;}
    virtual std::string toString() const override {
        return tfm::format(
                "Instance[\n"
                "  reference = %s,\n"
                "  transform = %s,\n"
                "  ref_name = %s,\n"
                "]",
                m_reference->toString(),
                m_localToWorld.toString(),
                m_ref_name
                );
    }
protected:
    Reference * m_reference;
    Transform m_worldToLocal;
    Transform m_localToWorld;
    std::string m_ref_name;
};
NORI_NAMESPACE_END
#endif
