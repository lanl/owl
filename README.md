# Description
**OWL: Open Wave Library for seismic wave modeling and full-waveform inversion**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20668968.svg)](https://doi.org/10.5281/zenodo.20668968)

`OWL` is an open-source package for seismic wave simulation and full-waveform inversion (FWI) in the following types of media:

- 2D/3D isotropic acoustic media, with or without a flat acoustic free surface at the top.
- 2D/3D isotropic or simple anisotropic elastic media (see below for the definition), with or without a flat elastic free surface at the top.
- 2D/3D isotropic or general anisotropic elastic media (see below for the definition), with or without a flat or topographic elastic free surface at the top.

Algorithm features:

- For seismic wave modeling and FWI in 2D/3D isotropic acoustic media, the solvers are based on a high-order standard staggered-grid finite-difference scheme (SSG-FD). The medium is parameterized with $(V_p, \rho)$.
- For seismic wave modeling and FWI in 2D/3D isotropic elastic media, or in simple anisotropic elastic media where $C_{ij} = 0$ ($ij=14,15,16,24,25,26,34,35,36,45,46,56$), the solvers are based on a high-order standard staggered-grid finite-difference scheme (SSG-FD).
    - In the isotropic case, the elastic medium is parameterized with $(V_p, V_s, \rho)$.
    - In the anisotropic case, the elastic medium can be parameterized with one of three methods:
        - Thomsen parameters: $(V_p, V_s, \rho, \varepsilon, \delta, \theta, \phi)$, where $\theta$ and $\phi$ can only be $N\pi/2$ (where $N=0, 1, 2$); this can describe anisotropy up to VTI (transverse isotropy with a vertical symmetry axis) or HTI (TI with a horizontal axis).
        - Alkhalifah-Tsvankin parameters: $(V_p, V_s, \rho, \varepsilon, \eta, \theta, \phi)$, where $\theta$ and $\phi$ can only be $N\pi/2$ (where $N=0, 1, 2$); this can describe anisotropy up to VTI or HTI.
        - Elasticity constants: $(C_{11}, C_{12}, C_{13}, C_{22}, C_{23}, C_{33}, C_{44}, C_{55}, C_{66})$; this can describe anisotropy up to VTI, HTI, or orthorhombic (ORT) anisotropy.
    - All of the above parameters can be spatially heterogeneous.
- For seismic wave modeling and FWI in 2D/3D isotropic or general anisotropic elastic media, the solvers are based on a high-order fully staggered-grid finite-difference scheme (FSG-FD, a.k.a. the Lebedev scheme).
    - In the isotropic case, the elastic medium is parameterized with $(V_p, V_s, \rho)$.
    - In the anisotropic case, the elastic medium can be parameterized with one of three methods:
        - Thomsen parameters: $(V_p, V_s, \rho, \varepsilon, \delta, \theta, \phi)$, where $\theta$ and $\phi$ can be arbitrary; this can describe anisotropy up to TI with a tilted axis (TTI).
        - Alkhalifah-Tsvankin parameters: $(V_p, V_s, \rho, \varepsilon, \eta, \theta, \phi)$, where $\theta$ and $\phi$ can be arbitrary; this can describe anisotropy up to TTI.
        - Elasticity constants: $C_{ij}$ with $i,j=1,2,\cdots, 6$ and $j \geq i$; this can describe anisotropy up to general anisotropy, including TTI, rotated ORT, monoclinic, and triclinic anisotropies.
    - All of the above parameters can be spatially heterogeneous.
    - The model can have a flat or topographic elastic free surface.
        - In both cases, the mesh will be refined in the near-surface region to properly simulate surface waves.
        - In the topographic free-surface case, a curvilinear mesh with near-surface refinement and support for a high-resolution topography map (up to the grid resolution) will be automatically generated.

- Various misfit functions, including:
    - $L_2$-norm waveform misfit
    - Envelope misfit
    - Zero-lag cross-correlation misfit
    - Weighted normalized deconvolution (a.k.a. AWI -- adaptive waveform inversion)
    - Local AWI
    - Phase/amplitude misfits based on generalized dynamic time warping (developed in this work).
- Various inversion schemes for FWI, including:
    - SD: steepest descent
    - NCG: nonlinear conjugate gradient
    - _l_-BFGS: limited-memory Broyden-Fletcher-Goldfarb-Shanno algorithm
    - Adam: adaptive moment estimation algorithm
- Various gradient preconditioning and regularization schemes for FWI.
- Built-in data processing and interpolation functions for automated modeling/FWI.
- User-friendly parameter input.

Currently, `OWL` does not support wave modeling or FWI in:
- Coupled acoustic-elastic media (i.e., coupled fluid and solid media)
- Visco-acoustic, visco-elastic media
- Poro-elastic, poro-visco-elastic media
- Thermoelastic media
- Unstructured meshes as used by FEM/SEM methods

These features may be included in future releases.

This work was supported by Los Alamos National Laboratory (LANL) Laboratory Directed Research and Development (LDRD) project 20240322ER. LANL is operated by Triad National Security, LLC, for the National Nuclear Security Administration (NNSA) of the U.S. Department of Energy (DOE) under Contract No. 89233218CNA000001. The research used high-performance computing resources provided by LANL's Institutional Computing program.

The code is released under LANL open source approval reference O4921.

# Requirements
`OWL` depends on [FLIT](https://github.com/lanl/flit).

Some examples in [example](example) use [RGM](https://github.com/lanl/rgm) to generate random geological models and [pymplot](https://github.com/lanl/pymplot) for plotting.

The code is written in Fortran + MPI + OpenMP. Currently, it can be compiled only with Intel's compiler suite, which is freely available through the [Intel HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html).

# Usage
To install `OWL`:

```
cd src
ruby install.rb
```

The compiled `OWL` executables will be placed in `bin`.

To rebuild:

```
cd src
ruby install.rb clean
```

Several examples are included in the `example` directory.

# License
&copy; 2025-2026. Triad National Security, LLC. All rights reserved.

This program is open source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
- Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Author
Kai Gao, <kaigao@lanl.gov>

We welcome feedback, bug reports, and ideas for improving `OWL`.

If you use this package in your research and find it useful, please cite it as:

* Kai Gao, 2025, OWL: Open Wave Library for seismic wave modeling and full-waveform inversion in acoustic and elastic media, URL: [github.com/lanl/owl](https://github.com/lanl/owl)
* Kai Gao, Jackson W. Saftner, Ting Chen, Ryan T. Modrak, 2026, OWL: Open Wave Library for seismic wave modeling and full-waveform inversion in acoustic and elastic media, _under review_ with GJI.

# Examples
Below are some of the examples included in the paper under review (LA-UR-26-20025).

<p align="center">
  <img src="doc/Figures/valles_caldera/wave.jpg" alt="Description" width="400">
</p>
<p align="center"><strong>Wavefield snapshot in a 3D elastic model with a topographic free surface.</strong></p>

<p align="center">
  <img src="doc/Figures/topo/vp_gt.jpg" alt="Description" width="200">
  <img src="doc/Figures/topo/vs_gt.jpg" alt="Description" width="200"><br>
  <img src="doc/Figures/topo/data.jpg" alt="Description" width="200"><br>
  <img src="doc/Figures/topo/vp_init.jpg" alt="Description" width="200">
  <img src="doc/Figures/topo/vs_init.jpg" alt="Description" width="200"><br>
  <img src="doc/Figures/topo/vp_l2.jpg" alt="Description" width="200">
  <img src="doc/Figures/topo/vs_l2.jpg" alt="Description" width="200"><br>
  <img src="doc/Figures/topo/vp_dtw.jpg" alt="Description" width="200">
  <img src="doc/Figures/topo/vs_dtw.jpg" alt="Description" width="200">
</p>
<p align="center"><strong>FWI in an elastic medium with a topographic free surface. From top to bottom: ground-truth Vp and Vs models, simulated data (vz component), models inverted by L2-norm FWI, and models inverted by GWI.</strong></p>
