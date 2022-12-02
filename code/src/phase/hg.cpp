
#include <kombu/phase.h>
#include <nori/warp.h>
#include <nori/frame.h>
NORI_NAMESPACE_BEGIN 
using namespace std;

class HenyeyGreensteinPhaseFunction : public PhaseFunction {
public:
    HenyeyGreensteinPhaseFunction(const PropertyList &props) {
        m_g = props.getFloat("g", 0.f);
    }

    float sample(PhaseFunctionQueryRecord &pRec, const Point2f &sample) const override{
        pRec.wo = Warp::squareToHenyeyGreenstein(sample, m_g);
        Frame incident_local_frame(-pRec.wi);
        pRec.wo = incident_local_frame.toWorld(pRec.wo);
        return eval(pRec);
    }

    // pdf is equal to eval
    float eval(PhaseFunctionQueryRecord &pRec) const override{
        Frame incident_local_frame(-pRec.wi);
        pRec.pdf = Warp::squareToHenyeyGreensteinPdf(incident_local_frame.toLocal(pRec.wo), m_g);
        return pRec.pdf;
    }
    
    std::string toString() const override {
        return tfm::format(
            "HenyeyGreensteinPhaseFunction[\n"
            "  g = %f\n"
            "]", m_g);
    }
    private:
        float m_g;
};

NORI_REGISTER_CLASS(HenyeyGreensteinPhaseFunction, "hg");
NORI_NAMESPACE_END