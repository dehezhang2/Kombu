# This is the main CMake build file used to compile Nori
cmake_minimum_required (VERSION 2.9)
project("kombu")

# This file included below takes care of compiling dependencies
# and setting compiler flags on different platforms
include("CMakeConfig.txt")

include_directories(ext)

# The following lines build the main executable. If you add a source
# code file to Nori, be sure to include it in this list.
add_executable(kombu

  # Header files
  include/nori/bbox.h
  include/nori/bitmap.h
  include/nori/block.h
  include/nori/bsdf.h
  include/nori/bvh.h
  include/nori/camera.h
  include/nori/color.h
  include/nori/common.h
  include/nori/dpdf.h
  include/nori/frame.h
  include/nori/gui.h
  include/nori/integrator.h
  include/nori/emitter.h
  include/nori/kdtree.h
  include/nori/mesh.h
  include/nori/object.h
  include/nori/parser.h
  include/nori/proplist.h
  include/nori/photon.h
  include/nori/ray.h
  include/nori/render.h
  include/nori/rfilter.h
  include/nori/sampler.h
  include/nori/scene.h
  include/nori/shape.h
  include/nori/texture.h
  include/nori/timer.h
  include/nori/transform.h
  include/nori/vector.h
  include/nori/warp.h

  include/kombu/medium.h 
  include/kombu/envmap.h
  include/kombu/disney_utils.h
  include/kombu/principledhelpers.h
  include/kombu/microfacet.h
  

  # Source code files
  src/bsdfs/dielectric.cpp
  src/bsdfs/diffuse.cpp
  src/bsdfs/mirror.cpp
  src/bsdfs/microfacet.cpp
  src/bsdfs/transparent.cpp
  src/bsdfs/disney_glass.cpp
  src/bsdfs/disney_principled.cpp
  src/bsdfs/bumpmap.cpp

  src/cameras/perspective.cpp
  
  src/data_structure/bvh.cpp
  src/data_structure/photon.cpp

  src/emitters/arealight.cpp
  src/emitters/pointlight.cpp
  src/emitters/envmap.cpp
  # src/emitters/envlight.cpp
  
  src/geometry/shape.cpp
  src/geometry/mesh.cpp
  src/geometry/sphere.cpp
  src/geometry/obj.cpp

  src/integrators/averagevisibility.cpp
  src/integrators/direct.cpp
  src/integrators/direct_ems.cpp
  src/integrators/direct_mats.cpp
  src/integrators/direct_mis.cpp
  src/integrators/normals.cpp
  src/integrators/path_mats.cpp
  src/integrators/path_mis.cpp
  src/integrators/photonmapper.cpp
  src/integrators/vol_path_mats.cpp
  src/integrators/vol_path_mis.cpp

  src/mediums/homogeneous.cpp

  src/phase/isotropic.cpp
  
  src/sample/independent.cpp
  src/sample/warp.cpp
  src/sample/stratified.cpp

  src/tests/ttest.cpp
  src/tests/chi2test.cpp

  src/textures/consttexture.cpp
  src/textures/checkerboard.cpp
  src/textures/imagetexture.cpp
  
  src/utils/bitmap.cpp
  src/utils/block.cpp
  src/utils/common.cpp
  src/utils/parser.cpp
  src/utils/proplist.cpp
  src/utils/rfilter.cpp

  src/core/gui.cpp
  src/core/main.cpp
  src/core/object.cpp
  src/core/render.cpp
  src/core/scene.cpp

)

# The following lines build the warping test application
add_executable(warptest
  include/nori/warp.h
  src/sample/warp.cpp
  src/tests/warptest.cpp
  src/bsdfs/microfacet.cpp
  src/core/object.cpp
  src/utils/proplist.cpp
  src/utils/common.cpp
)

add_executable(tonemapper
        include/nori/bitmap.h
        src/utils/bitmap.cpp
        src/utils/common.cpp
        src/utils/hdrToLdr.cpp)

# Nori depends on some libraries created in CMakeConfig.txt. The following two
# lines ensure that Nori is built *after* those libraries have been created.
add_dependencies(kombu OpenEXR_p)
add_dependencies(kombu nanogui_p)
add_dependencies(kombu tbb_p)
add_dependencies(kombu pugixml)
add_dependencies(warptest kombu)
add_dependencies(tonemapper kombu)

# Link to dependency libraries
target_link_libraries(kombu ${extra_libs})
target_link_libraries(warptest ${extra_libs})
target_link_libraries(tonemapper ${extra_libs})

# vim: set et ts=2 sw=2 ft=cmake nospell:
