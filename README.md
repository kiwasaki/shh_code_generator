# shh_code_generator

A C++ code generator to compute Spherical Harmonics Hessian used in  
["Spherical Lighting with Spherical Harmonics Hessian"](https://dl.acm.org/doi/10.1145/3721238.3730689)
SIGGRAPH 2025 Tehnical Paper Conference Track.


If you find it useful, we are very happy to be cited.
```bibtex 
@inproceedings{10.1145/3721238.3730689,
author = {Iwasaki, Kei and Dobashi, Yoshinori},
title = {Spherical Lighting with Spherical Harmonics Hessian},
year = {2025},
isbn = {9798400715402},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3721238.3730689},
doi = {10.1145/3721238.3730689},
abstract = {In this paper, we introduce a second-order derivative of spherical harmonics, spherical harmonics Hessian, and solid spherical harmonics, a variant of spherical harmonics, to the computer graphics community. These mathematical tools are used to develop an analytical representation of the Hessian matrix of spherical harmonics coefficients for spherical lights. We apply our analytic representation of the Hessian matrix to grid-based SH lighting rendering applications with many spherical lights that store the incident light field as spherical harmonics coefficients and their spatial gradient at sparse grid. We develop a Hessian-based error metric, with which our method automatically and adaptively subdivides the grid whether the interpolation using the spatial gradient is appropriate. Our method can be easily incorporated into the grid-based precomputed radiance transfer (PRT) framework with small additional storage. We demonstrate that our adaptive grid subdivided by using the Hessian-based error metric can substantially improve the rendering quality in equal-time grid construction.},
booktitle = {Proceedings of the Special Interest Group on Computer Graphics and Interactive Techniques Conference Conference Papers},
articleno = {33},
numpages = {10},
keywords = {Spherical Harmonics, Solid Spherical Harmonics, Spherical Harmonics Hessian, Spherical Light},
location = {
},
series = {SIGGRAPH Conference Papers '25}
}
```