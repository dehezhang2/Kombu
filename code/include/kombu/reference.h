#if !defined(__NORI_REFERENCE_H)
#define __NORI_REFERENCE_H
#include <nori/object.h>
#include <nori/shape.h>
NORI_NAMESPACE_BEGIN
class Reference : public Shape{
public:
    Reference(const PropertyList &propList) { }
    virtual void addChild(NoriObject *obj) override;

    virtual const BoundingBox3f &getBoundingBox() const override{
        return m_shape->getBoundingBox();
    }

    /// Is this mesh an area emitter?
    virtual bool isEmitter() const override{
        return m_shape->isEmitter();
    }

    /// Is this mesh a medium?
    virtual bool isMedium() const override{
        return m_shape->isMedium();
    }

    /// Return a pointer to an attached area emitter instance
    Emitter *getEmitter() override{
        return m_shape->getEmitter();
    }

    /// Return a pointer to an attached area emitter instance (const version)
    const Emitter *getEmitter() const override{
        return m_shape->getEmitter();
    }

    /// Return a pointer to an attached medium instance
    Medium *getMedium() override{
        return m_shape->getMedium();
    }

    /// Return a pointer to an attached medium instance (const version)
    const Medium *getMedium() const override{
        return m_shape->getMedium();
    }

    /// Return a pointer to the BSDF associated with this mesh
    const BSDF *getBSDF() const { return m_shape->getBSDF(); }

    /// Return the total number of primitives in this shape
    virtual uint32_t getPrimitiveCount() const override{
        return m_shape->getPrimitiveCount();
    }

    //// Return an axis-aligned bounding box containing the given triangle
    virtual BoundingBox3f getBoundingBox(uint32_t index) const override{
        return m_shape->getBoundingBox(index);
    }

    //// Return the centroid of the given triangle
    virtual Point3f getCentroid(uint32_t index) const override{
        return m_shape->getCentroid(index);
    }

    //// Ray-Shape intersection test
    virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const override{
        return m_shape->rayIntersect(index, ray, u, v, t);
    }

    /// Set the intersection information: hit point, shading frame, UVs, etc.
    virtual void setHitInformation(uint32_t index, const Ray3f &ray, Intersection & its) const override{
        m_shape->setHitInformation(index, ray, its);
    }

    /**
     * \brief Sample a point on the surface (potentially using the point sRec.ref to importance sample)
     * This method should set sRec.p, sRec.n and sRec.pdf
     * Probability should be with respect to area
     * */
    virtual void sampleSurface(ShapeQueryRecord & sRec, const Point2f & sample) const override {
        m_shape->sampleSurface(sRec, sample);
    }
    virtual float pdfSurface(const ShapeQueryRecord & sRec) const override{
        return m_shape->pdfSurface(sRec);
    }
    virtual bool isReference() const override{return true;}
    virtual std::string toString() const override {
        return tfm::format(
                "Reference[\n"
                "  m_shape = %s,\n"
                "]",
                m_shape->toString());
    }
protected:
    Shape * m_shape = nullptr;
};
NORI_NAMESPACE_END
#endif
