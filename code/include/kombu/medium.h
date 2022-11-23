#if !defined(__NORI_MEDIUM_H)
#define __NORI_MEDIUM_H

#include <nori/object.h>
#include <kombu/phase.h>

NORI_NAMESPACE_BEGIN
/**
 * \brief Data record for conveniently querying and sampling the
 * a point in the medium
 */
struct MediumQueryRecord {
    /// Origin point
    Point3f ref;
    /// Point sampled 
    Point3f p;
    /// incoming Direction
    Vector3f wi;
    /// Sampled direction
    Vector3f wo;
    /// Probability of sampling
    float pdf;  
    /// Max free path length
    float tMax;
    /// whether the sampled point hits the medium (or surface)
    bool hitMedium;

    /// Empty constructor
    MediumQueryRecord(){}
    /// Sample phase function and distance (free path)
    MediumQueryRecord(const Point3f &_ref, const Vector3f &_wi) : ref(_ref), wi(_wi){}
    /// Evaluate albedo
    MediumQueryRecord(const Point3f &_ref) : ref(_ref){}
    /// Evaluate transmittance
    MediumQueryRecord(const Point3f &_ref, const Point3f &_p) : ref(_ref), p(_p){
        wo = (p - ref).normalized();
    }
    /// Query probability density of sampling
    MediumQueryRecord(const Point3f &_ref, const Point3f &_p, const Vector3f &_wi) : ref(_ref), p(_p), wi(_wi){
        wo = (p - ref).normalized();
    }
};

/**
 * \brief Superclass of all mediums
 */
class Medium: public NoriObject{
public:
/**
	 * \brief Sample the medium distance
	 *
	 * \param mRec    A medium query record (ref, wi are needed)
	 * \param sample  A uniformly distributed sample on \f$[0,1]^2\f$
	 *
	 * \return The albedo, evaluated for each color channel.
	 *         A zero value means that sampling failed.
	 */
    virtual Color3f sample(MediumQueryRecord &mRec, const Point2f &sample) const = 0;
    
    /**
	 * \brief Evaluate the medium albedo
	 *
	 * \param mRec
	 *    A medium query record (only ref is needed)
	 * \return
	 *     The albedo, evaluated for each color channel.
	 */
    virtual Color3f eval(const MediumQueryRecord &mRec) const = 0;

    /**
	 * \brief Compute the probability of sampling \c mRec.p.
	 *
	 * This method provides access to the probability density that
	 * is realized by the \ref sample() method.
	 *
	 * \param mRec
	 *    A medium query record (ref, p, wi are needed)
	 *
	 * \return
	 *     A probability/density value
	 */
	virtual float pdf(const MediumQueryRecord &mRec) const = 0;

    /**
     * @brief Check whether the point is in the medium
     * 
     * @param p 
     *    The point you want to check
     * @return true : inside
     * @return false : outside
     */

    virtual bool contains(Point3f &p) const = 0;

     /**
	 * \brief Evaluate transmittance from \c mRec.ref to \c mRec.p.
	 *
	 * This method provides access to the transmittance bewteen two points
	 * this is used for ray marching. 
	 *
	 * \param mRec
	 *    A medium query record (ref, p are needed)
	 *
	 * \return
	 *    The transmittance between the two points
	 */
    virtual Color3f evalTransmittance(const MediumQueryRecord &mRec) const = 0;

    // Return the phase function of this medium
    inline const PhaseFunction *getPhaseFunction() const { return m_phase; }
    
    /**
     * \brief Return whether or not this medium is homogeneous.
     */
    virtual bool isHomogeneous() const { return false; }

    virtual Color3f &getSigmaA() const = 0;
    virtual Color3f &getSigmaS() const = 0;
    virtual Color3f &getSigmaT() const = 0;
    virtual Color3f &getAlbedo() const = 0;

    /**
     * \brief Set the shape if the medium is attached to a shape
     * */
    void setShape(Shape * shape) { m_shape = shape; }
    
    virtual void addChild(NoriObject *obj) {}

    virtual std::string toString() const = 0;

    EClassType getClassType() const override { return EMedium; }



protected:
    PhaseFunction * m_phase = nullptr;
    Shape * m_shape = nullptr;
    // Color3f m_sigmaA;
    // Color3f m_sigmaS;
    // Color3f m_sigmaT;
};

NORI_NAMESPACE_END
#endif