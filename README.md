# The Kombu Renderer

[**Deheng Zhang**](https://github.com/dehezhang2), [**Ganlin Zhang**](https://github.com/zhangganlin)

## Overview

<img src="./Project_page/assets/earth_water_2.0.png" alt="earth_water_2.0" style="zoom:10%;" />

This is a physically based renderer based on Nori in C++11. The functions supported include path tracing with multiple importance sampling (MIS), photon mapping, volumetric path tracing with MIS, heterogeneous media, different distance sampling and transmittance estimation methods (ray marching, delta tracking, and ratio tracking), bilateral filter with uniform variance denoting, directional light, anisotropic phase function, object instancing, Disney BSDF, environment map and texture, stratified sampling, blend and conductor BSDF, bump mapping.  For more information, please view our [project website](). 

The renderer is developed for the rendering competition in [ETH Computer Graphics course 2022](https://cgl.ethz.ch/teaching/cg22/home.php). **Please DO NOT directly copy the code from this repository**, we only provide the code for reference.

## Installation

* All the codes are tested in the following environment:

  - Mac OS (Version 11.6), Linux (Ubuntu 22.04)

  * CMake (Version 3.20.2, Version 3.25.2)

* Install Qt4: You need to install Qt4 for the heterogeneous volume data reading:

  * Mac

    ```bash
    conda remove --force qt
    brew uninstall --force qt
    brew tap cartr/qt4
    brew tap-pin cartr/qt4
    brew install qt@4
    conda install -c pkgw-forge qt4
    ```

  * Linux: follow the instruction in this [link](https://ubuntuhandbook.org/index.php/2020/07/install-qt4-ubuntu-20-04/).

* Build the code (or you can directly use the compiled release version)

  ```bash
  mkdir build
  cd build
  cmake ..
  make -j 4
  ```

## Data Preparation

Please download the scene from [polybox](https://polybox.ethz.ch/index.php/s/FuiuW27Vocy9fjF), and the layout should looks like this:

```

└── final_scene
    ├── bounds.obj
    ├── earth_water_4.0.exr
    ├── earth_water_5.0.exr
    ├── earth_water_5.0.png
    ├── earth_waterdenoised_4.0.exr
    ├── earth_waterdenoised_5.0.exr
    ├── earth_waterdenoised.exr
    ├── earth_waterdenoised.png
    ├── earth_water.exr
    ├── earth_water.png
    ├── earth_water.xml
    ├── light.obj
    ├── meshes
    ├── obj_chimney_meshes
    ├── obj_meshes
    ├── smoke_exp
    ├── smoke.vol
    └── texture
```

## Demo

* Run the renderer using:

```bash
cd build
./kombu ../final_scene/earth_water.xml
```

## License

The code is released under the GPL-3.0 license.
