#include <kombu/denoiser.h>

NORI_NAMESPACE_BEGIN
class BilaterialDenoiser : public Denoiser{
public:
    BilaterialDenoiser(const PropertyList & propList) {
        m_windowSize = propList.getInteger("window", 5);
        m_sigma_i = propList.getFloat("sigma_i", 12.0f);
        m_sigma_d = propList.getFloat("sigma_d", 16.0f);
        m_k = propList.getFloat("k", 1.f);

    }
    void denoise(Bitmap* bitmap) const override{
        Vector2i m_size(bitmap->cols(), bitmap->rows());
        Bitmap target(m_size);
        int height = m_size.y();
        int width = m_size.x();
        int half = m_windowSize / 2;
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

    void applyFilter(Bitmap* source, Bitmap* target, int x, int y) const override{
        Color3f iFiltered = 0;
        double wP = 0;
        int neighbor_x = 0;
        int neighbor_y = 0;
        int half = m_windowSize / 2;
        for(int i = 0; i < m_windowSize; i++) {
            for(int j = 0; j < m_windowSize; j++) {
                neighbor_x = x - (half - i);
                neighbor_y = y - (half - j);
                Color3f center = source->coeff(x, y), neighbor = source->coeff(neighbor_x, neighbor_y);
                Color3f diff = center - neighbor;
                double distance_i = diff.x() * diff.x() + diff.y() * diff.y() + diff.y() * diff.y();
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
            "BilaterialDenoiser[\n"
            "  sigmaI = %s,\n"
            "  sigmaD = %s,\n"
            "  window size = %s,\n"
            "]",
            m_sigma_i,
            m_sigma_d,
            m_windowSize
        );
    }
protected:
    float m_sigma_i;
    float m_sigma_d;
    float m_k;
};
NORI_REGISTER_CLASS(BilaterialDenoiser, "bilaterial");
NORI_NAMESPACE_END
