
#include <kombu/phase.h>
#include <nori/warp.h>
NORI_NAMESPACE_BEGIN 
using namespace std;

class HenyeyGreensteinPhaseFunction : public PhaseFunction {
public:
    HenyeyGreensteinPhaseFunction(const PropertyList &props) {
    }

    float sample(PhaseFunctionQueryRecord &pRec, const Point2f &sample) const override{
        pRec.wo = Warp::squareToUniformSphere(sample);
        return INV_FOURPI;
    }

    // pdf is equal to eval
    float eval(PhaseFunctionQueryRecord &pRec) const override{
        return INV_FOURPI;
    }
    
    std::string toString() const override {
        return tfm::format("IsotropicPhaseFunction[]");
    }
};

NORI_REGISTER_CLASS(HenyeyGreensteinPhaseFunction, "hg");
NORI_NAMESPACE_END