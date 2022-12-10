/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Romain Prévost

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/emitter.h>
#include <kombu/envmap.h>
#include <nori/frame.h>
#include <nori/shape.h>


NORI_NAMESPACE_BEGIN

class EnvironmentEmitter : public Emitter {
public:
    EnvironmentEmitter(const PropertyList &props) {
    	//load an exr file
    	m_envmap = Envmap(props.getString("path2map"));
    }

    virtual std::string toString() const override {
        return tfm::format(
                "EnvironmentLight[\n"
                "]"
        		);
    }

    virtual Color3f eval(const EmitterQueryRecord & lRec) const override {
    	return m_envmap.eval(lRec.wi);
    }

    virtual Color3f sample(EmitterQueryRecord & lRec, const Point2f & sample) const override {
    	Color3f col = m_envmap.sample(lRec.wi, sample);
    	lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, 10000);
    	return col;
    }

    virtual float pdf(const EmitterQueryRecord &lRec) const override {
    	return m_envmap.pdf(lRec.wi);
    }

    /*
    Point2f dirToPixel(const Vector3f &vec) const{
		float theta = acos(Frame::cosTheta(vec));
		float phi = acos(Frame::cosPhi(vec));;

		float u = theta * INV_PI * (m_envmap.rows() -1);
		float v = phi  * INV_PI * 0.5 * (m_envmap.cols() - 1);

		std::cout << theta << ", " << phi << "\n";

		return Point2f(u,v);
    }

    Vector3f pixelToDir(const Point2f &p) const {
    	float u = p.x();
    	float v = p.y();

    	float theta = u * M_PI / (m_envmap.rows() -1);
    	float phi = v * 2 * M_PI / (m_envmap.cols() - 1);

    	return Vector3f(sin(theta) * cos(phi), sin(theta)* sin(phi), cos(theta)).normalized();
    }
    */


private:
    Envmap m_envmap;
};

NORI_REGISTER_CLASS(EnvironmentEmitter, "environment");
NORI_NAMESPACE_END
