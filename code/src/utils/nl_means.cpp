#include <kombu/denoiser.h>

NORI_NAMESPACE_BEGIN
class NLMeansDenoiser : public Denoiser{
public:
    NLMeansDenoiser(const PropertyList & propList) {
        m_f = propList.getInteger("f", 5);
        m_r = propList.getInteger("r", 3);
        m_sigma_i = propList.getFloat("sigma_i", 12.0f);
        m_sigma_d = propList.getFloat("sigma_d", 16.0f);
        m_k = propList.getFloat("k", 1.f);

    }
    void denoise(Bitmap* bitmap) const override{
        Vector2i m_size(bitmap->cols(), bitmap->rows());
        Bitmap target(m_size);
        int height = m_size.y();
        int width = m_size.x();
        int half = m_f + m_r;
        for(int i = half; i < height - half; i++) {
            for(int j = half; j < width - half; j++) {
                applyFilter(bitmap, &target, i, j);
            }
        }
        for(int i = half; i < height - half; i++) {
            for(int j = half; j < width - half; j++) {
                bitmap->coeffRef(i, j) = target.coeff(i, j);
            }
        }
    }

    float computePatchDistance(const Bitmap* source, int x0, int y0, int x1, int y1) const {
        int m_windowSize = 2 * m_f + 1;
        int neighbor_x0 = 0;
        int neighbor_y0 = 0;
        int neighbor_x1 = 0;
        int neighbor_y1 = 0;
        double distance = 0.f;
        for(int i = 0; i < m_windowSize; i++){
            for(int j = 0; j < m_windowSize; j++){
                neighbor_x0 = x0 - (m_f - i);
                neighbor_y0 = y0 - (m_f - j);
                neighbor_x1 = x1 - (m_f - i);
                neighbor_y1 = y1 - (m_f - j);
                Color3f center0 = source->coeff(neighbor_x0, neighbor_y0), center1 = source->coeff(neighbor_x1, neighbor_y1);
                Color3f diff = center0 - center1;
                distance += diff.x() * diff.x() + diff.y() * diff.y() + diff.y() * diff.y();
            }
        }
        distance /=(float) (m_windowSize * m_windowSize);
        return distance;
    }

    void applyFilter(Bitmap* source, Bitmap* target, int x, int y) const override{
        Color3f iFiltered = 0;
        double wP = 0;
        int neighbor_x = 0;
        int neighbor_y = 0;
        int m_windowSize = 2 * m_r + 1;

        for(int i = 0; i < m_windowSize; i++) {
            for(int j = 0; j < m_windowSize; j++) {
                neighbor_x = x - (m_r - i);
                neighbor_y = y - (m_r - j);
                double distance_i = computePatchDistance(source, x, y, neighbor_x, neighbor_y);
                double gi = exp(
                    -( 
                        (distance_i - 2 * m_sigma_i * m_sigma_i )/ (Epsilon + m_k * m_k * 2 * m_sigma_i * m_sigma_i)
                    )
                );
                double gd = exp(
                    -distance_pos(x, y, neighbor_x, neighbor_y)/(2 * m_sigma_d * m_sigma_d + Epsilon)
                );
                double w = gi * gd;
                iFiltered = iFiltered + source->coeff(neighbor_x, neighbor_y) * w;
                wP = wP + w;
            }
        }
        iFiltered = iFiltered / wP;
        target->coeffRef(x, y) = iFiltered;
    }

    std::string toString() const override {
        return tfm::format(
            "NLMeansDenoiser[\n"
            "  sigmaI = %s,\n"
            "  sigmaD = %s,\n"
            "  k = %s,\n"
            "  f = %s,\n"
            "  r = %s,\n"
            "]",
            m_sigma_i,
            m_sigma_d,
            m_k,
            m_f,
            m_r
        );
    }
protected:
    float m_sigma_i;
    float m_sigma_d;
    float m_k;
    float m_f;
    float m_r;
};
NORI_REGISTER_CLASS(NLMeansDenoiser, "nl_means");
NORI_NAMESPACE_END
