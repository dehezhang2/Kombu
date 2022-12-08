#if !defined(__NORI_DENOISER_H)
#define __NORI_DENOISER_H

#include <nori/object.h>
#include <nori/bitmap.h>
NORI_NAMESPACE_BEGIN

/**
 * \brief Superclass of all denoisers
 */
class Denoiser: public NoriObject{
public:
    virtual void denoise(Bitmap* bitmap) const = 0;
    virtual void applyFilter(Bitmap* source, Bitmap* target, int x, int y) const = 0;
    virtual std::string toString() const = 0;
    virtual EClassType getClassType() const override { return EDenoiser; }
    virtual float distance_pos(int x, int y, int i, int j) const{
       return float(std::pow(x - i, 2) + std::pow(y - j, 2));
    }
protected:
    int m_windowSize;
};

NORI_NAMESPACE_END
#endif