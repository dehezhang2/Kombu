
#include <kombu/phase.h>
#include <nori/warp.h>
NORI_NAMESPACE_BEGIN 
using namespace std;

class IsotropicPhaseFunction : public PhaseFunction {
public:
    IsotropicPhaseFunction(const PropertyList &props) {
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

NORI_REGISTER_CLASS(IsotropicPhaseFunction, "isotropic");
NORI_NAMESPACE_END