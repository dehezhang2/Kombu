#if !defined(__NORI_PHASE_H)
#define __NORI_PHASE_H

#include <nori/object.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN
/**
 * \brief Convenience data structure used to pass multiple
 * parameters to the evaluation and sampling routines in \ref PhaseFunction
 */

struct PhaseFunctionQueryRecord {
    /// Incident direction (in the local frame)
    Vector3f wi;

    /// Outgoing direction (in the local frame)
    Vector3f wo;

    /// pdf of the query
    float pdf;

    /// Measure associated with the sample
    EMeasure measure;

    /// Empty record
    PhaseFunctionQueryRecord(){}
    
    /// Create a new record for sampling the PhaseFunction
    PhaseFunctionQueryRecord(const Vector3f &wi)
        : wi(wi), measure(EUnknownMeasure) { }
    
    /// Create a new record for querying the PhaseFunction
    PhaseFunctionQueryRecord(const Vector3f &wi,
            const Vector3f &wo, EMeasure measure)
        : wi(wi), wo(wo), measure(measure) { }
};

class PhaseFunction : public NoriObject {
public:
    PhaseFunction();
    PhaseFunction(const PropertyList &props);
    /// Release all memory
    virtual ~PhaseFunction() {}

    virtual float sample(PhaseFunctionQueryRecord &pRec, const Point2f &sample) const = 0;

    // pdf is equal to eval
    virtual float eval(PhaseFunctionQueryRecord &pRec) const = 0;

    EClassType getClassType() const { return EPhaseFunction; }

protected:
};

NORI_NAMESPACE_END
#endif