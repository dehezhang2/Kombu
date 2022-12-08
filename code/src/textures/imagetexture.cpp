
#include <nori/bitmap.h>
#include <nori/object.h>
#include <nori/texture.h>
#include <kombu/envmap.h>
#include <filesystem/resolver.h>
NORI_NAMESPACE_BEGIN


// Color3f colLerp(float t, Color3f v1, Color3f v2) {
// 	return (1 - t) * v1 + t * v2;
// }

// // bilinear interpolation
// Color3f bilerp(float tx, float ty, Color3f v00, Color3f v01, Color3f v10, Color3f v11) {
// 	Color3f up = colLerp(tx, v00, v01);
// 	Color3f bottom = colLerp(tx, v10, v11);
// 	return colLerp(ty, up, bottom);
// }


//png images to rgb vector
class ImageTexture : public Texture<Color3f> {
public:

    ImageTexture(const PropertyList &props) {
        m_name = props.getString("fileName", "");
        if(m_name[0]=='/'){
            m_image = Bitmap(m_name);
        }
        else{
            filesystem::path filePath = getFileResolver()->resolve(m_name);
            m_image = Bitmap(filePath.str());
        }


        m_scale = props.getVector2("scale", Vector2f(1));
        
        m_width = m_image.cols();
        m_height = m_image.rows();
    }


    //this function should return a rgb value depending on the coordinates UV passed
    Color3f eval(const Point2f & uv) {
        

        float uf = (uv.y() * m_height) / m_scale.y(); 
        float vf = ((1-uv.x()) * m_width) / m_scale.x(); 

        int u = floor(uf);
        int v = floor(vf);


        // in corners and on borders
        if (u >= m_image.rows() - 1 || v >= m_image.cols() - 1) {
            return m_image(m_image.rows() - 1,m_image.cols() - 1);
        }

        float du = uf - u;
        float dv = vf - v;

        // return 1.f;

        Color3f v00 = m_image(u, v);
        Color3f v01 = m_image(u, v+1);
        Color3f v10 = m_image(u+1, v);
        Color3f v11 = m_image(u+1, v+1);

        return bilerp(du, dv, v00, v01, v10, v11);
    }

    std::string toString() const {
        return tfm::format(
                "ImageTexture[\n"
                "  m_name = %s,\n"
                "  scale = %s,\n"
                "]",
                m_name,
                m_scale.toString()
        );
    }

protected:
    Vector2f m_scale;
    std::string m_name;
    // std::vector<unsigned char> m_image;
    Bitmap m_image;

    unsigned m_width;
    unsigned m_height;
};

NORI_REGISTER_CLASS(ImageTexture, "imagetexture")
NORI_NAMESPACE_END