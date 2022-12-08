/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob, Romain Prévost

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

#include <nori/scene.h>
#include <nori/bitmap.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/camera.h>
#include <nori/emitter.h>
#include <kombu/instance.h>

NORI_NAMESPACE_BEGIN

Scene::Scene(const PropertyList &) {
    m_bvh = new BVH();
}

Scene::~Scene() {
    delete m_bvh;
    delete m_sampler;
    delete m_camera;
    delete m_integrator;
    for(auto e : m_emitters)
        delete e;
    m_emitters.clear();
}

void Scene::activate() {
    // link the references and the instances
    for(auto& instance : m_instances){
        bool found = false;
        for(auto ref: m_references){
            if(ref->getIdName()==instance->getRefName()){
                found = true;
                instance->setRef(ref);
            }
        }
        if(!found) throw NoriException("Reference not found!");
        m_bvh->addShape(instance);
        m_shapes.push_back(instance);
        if(instance->isEmitter())
            m_emitters.push_back(instance->getEmitter());
        if(instance->isMedium()){
            m_mediums.push_back(instance->getMedium());
        }
    }

    m_bvh->build();

    if (!m_integrator)
        throw NoriException("No integrator was specified!");
    if (!m_camera)
        throw NoriException("No camera was specified!");
    
    if (!m_sampler) {
        /* Create a default (independent) sampler */
        m_sampler = static_cast<Sampler*>(
            NoriObjectFactory::createInstance("independent", PropertyList()));
        m_sampler->activate();
    }

    cout << endl;
    cout << "Configuration: " << toString() << endl;
    cout << endl;
}

void Scene::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EMesh: {
                Shape *mesh = static_cast<Shape *>(obj);
                if(mesh->isReference()){
                    Reference *ref = static_cast<Reference*>(obj);
                    m_references.push_back(ref);
                } else if (mesh->isInstance()){
                    Instance *ins = static_cast<Instance*>(obj);
                    m_instances.push_back(ins);
                } else {
                    m_bvh->addShape(mesh);
                    m_shapes.push_back(mesh);
                    if(mesh->isEmitter())
                        m_emitters.push_back(mesh->getEmitter());
                    if(mesh->isMedium())
                        m_mediums.push_back(mesh->getMedium());
                }
            }
            break;
        
        case EEmitter:
            m_emitters.push_back(static_cast<Emitter *>(obj));
            break;

        case ESampler:
            if (m_sampler)
                throw NoriException("There can only be one sampler per scene!");
            m_sampler = static_cast<Sampler *>(obj);
            break;

        case ECamera:
            if (m_camera)
                throw NoriException("There can only be one camera per scene!");
            m_camera = static_cast<Camera *>(obj);
            break;
        
        case EIntegrator:
            if (m_integrator)
                throw NoriException("There can only be one integrator per scene!");
            m_integrator = static_cast<Integrator *>(obj);
            break;

        case EMedium:
            {
                m_mediums.push_back(static_cast<Medium *>(obj));
                break;
            }
            break;

        case EDenoiser:
            {
                if (m_denoiser)
                    throw NoriException("There can only be one denoiser per scene!");
                m_denoiser = static_cast<Denoiser *>(obj);
            }
            break;

        default:
            throw NoriException("Scene::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
    }
}

std::string Scene::toString() const {
    std::string shapes;
    for (size_t i=0; i<m_shapes.size(); ++i) {
        shapes += std::string("  ") + indent(m_shapes[i]->toString(), 2);
        if (i + 1 < m_shapes.size())
            shapes += ",";
        shapes += "\n";
    }

    std::string lights;
    for (size_t i=0; i<m_emitters.size(); ++i) {
        lights += std::string("  ") + indent(m_emitters[i]->toString(), 2);
        if (i + 1 < m_emitters.size())
            lights += ",";
        lights += "\n";
    }

    std::string mediums;
    for (size_t i=0; i<m_mediums.size(); ++i) {
        mediums += std::string("  ") + indent(m_mediums[i]->toString(), 2);
        if (i + 1 < m_mediums.size())
            mediums += ",";
        mediums += "\n";
    }

    return tfm::format(
        "Scene[\n"
        "  integrator = %s,\n"
        "  sampler = %s\n"
        "  camera = %s,\n"
        "  shapes = {\n"
        "  %s  }\n"
        "  emitters = {\n"
        "  %s  }\n"
        "  mediums = {\n"
        "  %s  }\n"
        "]",
        indent(m_integrator->toString()),
        indent(m_sampler->toString()),
        indent(m_camera->toString()),
        indent(shapes, 2),
        indent(lights,2),
        indent(mediums,2)
    );
}

NORI_REGISTER_CLASS(Scene, "scene");
NORI_NAMESPACE_END
