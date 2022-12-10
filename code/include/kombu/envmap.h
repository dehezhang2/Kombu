
#if !defined(__NORI_ENVMAP_H)
#define __NORI_ENVMAP_H

#include <nori/bitmap.h>
#include <nori/frame.h>
#include <nori/object.h>


NORI_NAMESPACE_BEGIN

inline Color3f colLerp(float t, Color3f v1, Color3f v2) {
	return (1 - t) * v1 + t * v2;
}

// bilinear interpolation
inline Color3f bilerp(float tx, float ty, Color3f v00, Color3f v01, Color3f v10, Color3f v11) {
	Color3f up = colLerp(tx, v00, v01);
	Color3f bottom = colLerp(tx, v10, v11);
	return colLerp(ty, up, bottom);
}
NORI_NAMESPACE_END

#endif /* __NORI_ENVMAP_H */
