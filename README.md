# Description
**OWL: Open Wave Library for seismic wave modeling and full-waveform inversion**

`OWL` is an open-source package written mostly in Fortran for seismic wave simulation and performing full-waveform inversion (FWI). It supports waveform modeling and FWI in the following types of media: 

- 2D/3D isotropic acoustic media, with or without an acoustic free surface at the top. 
- 2D/3D isotropic elastic media, with or without an elastic free surface at the top. 
- 2D/3D general anisotropic elastic media, with or without a flat or topographic elastic free surface at the top. 

Algorithm features: 

- For seismic wave modeling and FWI in 2D/3D isotropic acoustic media, the solvers are based on high-order standard staggered-grid finite-difference scheme (SSG-FD). The medium is parameterized with $(V_p, \rho)$. 
- For seismic wave modeling and FWI in 2D/3D isotropic elastic media, or anisotropic elastic media where $C_{ij} = 0$ (where $ij=14,15,16,24,25,26,34,35,36,45,46,56$), the solvers are based on high-order standard staggered-grid finite-difference scheme (SSG-FD). 
    - In the isotropic case, the elastic medium is parameterized with $(V_p, V_s, \rho)$. 
    - In the anisotropic case, the elastic medium can be parameterized with one of the three methods:
        - Thomsen parameters: $(V_p, V_s, \rho, \varepsilon, \delta, \theta, \phi)$, where $\theta$ and $\phi$ can only be $N\pi/2$ (where $N=0, 1, 2$). 
        - Alkhalifah-Tsvankin parameters: $(V_p, V_s, \rho, \varepsilon, \eta, \theta, \phi)$, where $\theta$ and $\phi$ can only be $N\pi/2$ (where $N=0, 1, 2$). 
        - Elasticity constants: $(C_{11}, C_{12}, C_{13}, C_{22}, C_{23}, C_{33}, C_{44}, C_{55}, C_{66})$. 
    - All the above parameters can be spatially heterogeneous. 
- For seismic wave modeling and FWI in 2D/3D isotropic or general anisotropic media, the solvers are based on high-order fully staggered-grid finite-difference scheme (FSG-FD, a.k.a. Lebedev scheme). 
    - In the isotropic case, the elastic medium is parameterized with $(V_p, V_s, \rho)$. 
    - In the anisotropic case, the elastic medium can be parameterized with one of the three methods:
        - Thomsen parameters: $(V_p, V_s, \rho, \varepsilon, \delta, \theta, \phi)$, where $\theta$ and $\phi$ can be arbitrary. 
        - Alkhalifah-Tsvankin parameters: $(V_p, V_s, \rho, \varepsilon, \eta, \theta, \phi)$, where $\theta$ and $\phi$ can be arbitrary. 
        - Elasticity constants: $C_{ij}$ with $i,j=1,2,\cdots, 6$ and $j \geq i$. 
    - All the above parameters can be spatially heterogeneous. 
    - The model can have a flat or topographic elastic free surface. 
        - In both cases, the mesh will be refined in the near-surface region to properly simulate surface wave. 
        - In the topographic free surface case, a curvilinear mesh with refined near-surface mesh, and supports high-resolution topography map (up to grid resolution), will be automatically generated. 

- Various FWI inversion schemes, including steepest descent (SD), nonlinear conjugate gradient (NCG), limited-memory BFGS (l-BFGS), and Adam inversion schemes for FWI. 
- Various regularization schemes for FWI. 
- Various adjoint sources including l2-norm, envelope, zero-lag cross-correlation, weighted normalized deconvolution, phase, traveltime, and generalized dynamic time warping misfit functions. 
- Built-in data processing and interpolation functionalities for automated modeling/FWI. 
- User-friendly parameter input. 

Currently, `OWL' does not yet support:

- Coupled acoustic-elastic media.
- Viscoacoustic/viscoelastic media. 
- Unstructured mesh. 

These features may be included in future development pending separate approval. 

The work was supported by Los Alamos National Laboratory (LANL) Laboratory Directory Research and Development (LDRD) project 20240322ER. LANL is operated by Triad National Security, LLC, for the National Nuclear Security Administration (NNSA) of the U.S. Department of Energy (DOE) under Contract No. 89233218CNA000001. The research used high-performance computing resources provided by LANL's Institutional Computing program. 

The work is under LANL open source approval reference O4921.

# Reference
A manual will be released in the near future. 

# Requirement
`OWL` depends on [FLIT](https://github.com/lanl/flit). Some examples in [example](example) use [RGM](https://github.com/lanl/rgm) to generate random geological models. 

The code is written in Fortran + MPI. Currently, it can only be compiled with Intel's compiler suite, which is freely available at [Intel HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html). 

# Use
To install `OWL`, 

```
cd src
ruby install.rb
```

The compiled `OWL` executables will be at `bin`.

To remake, 

```
cd src
ruby install.rb clean
```

We include several simple examples to use `OWL` in [example](example). To run the tests, 

```
cd test
```

and the scripts to reproduce the examples in the mansucript are contained in subfolders. 

# License
&copy; 2025. Triad National Security, LLC. All rights reserved. 

This program is Open-Source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
- Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Author
Kai Gao, <kaigao@lanl.gov>

We welcome feedbacks, bug reports, and improvement ideas on `OWL`. 

If you use this package in your research and find it useful, please cite it as

* Kai Gao, 2025, OWL: Open Wave Library for seismic wave modeling and full-waveform inversion in acoustic and elastic media, url: [github.com/lanl/owl](https://github.com/lanl/owl)
* Kai Gao, Jackson W. Saftner, Ting Chen, Ryan T. Modrak, 2025, OWL: Open Wave Library for seismic wave modeling and full-waveform inversion in acoustic and elastic media, _under review_ with GJI. 
