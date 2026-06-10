
# Short Tutorial for OWL

This is a short tutorial on how to install and use `OWL` for forward wavefield modeling and FWI.

## Table of Contents
- [Installation](#installation)
- [Parameters](#parameters)
	- [Geometry](#geometry)
	- [Dimension](#dimension)
	- [Medium properties](#medium-properties)
	- [Data and data processing](#data-and-data-processing)
	- [Free surface and topography](#free-surface-and-topography)
	- [Inversion](#inversion)
	- [Other parameters](#other-parameters)
- [Examples](#examples)


## Installation

Installation of `OWL` is straightforward:

```bash
git clone https://github.com/lanl/owl.git
cd owl
ruby install.rb clean
```

You need to install [`FLIT`](https://github.com/lanl/flit) before installing `OWL`. To reproduce examples in the `example` directory, you need to install [`RGM`](https://github.com/lanl/rgm) and [`pymplot`](https://github.com/lanl/pymplot). 

## Parameters

### Geometry

To use `OWL` for forward modeling or FWI, a source-receiver geometry file must be provided. The format of `OWL`'s geometry file is:

<!-- #### > Overall geometry file -->

`OWL` requires an overall geometry file (say, `./geometry/geometry.txt`) in the following form:
 
```ruby
	shot_1_geometry.txt
	shot_2_geometry.txt
	shot_3_geometry.txt
	....
```

where each row is the name of one common-shot gather's geometry file. The names of the individual geometry files can be arbitrary, e.g., 

```ruby
	shot_1_geometry.txt
	source_222_geometry.txt
	event_33_sr.txt
	...
```

However, these geometry files should be distinct (otherwise, it does not make much sense to include two identical common-shot gathers in one survey) and should be present in the same directory as the overall geometry file, i.e., the path to an individual geometry file is `./geometry/shot_1_geometry.txt`, ...

<!-- #### > Individual geometry file -->

An individual geometry file describes the source-receiver distribution as well as source parameters. The form is:

```ruby
shot id

number of point sources
point source 1 location
point source 1 mechanism
point source 1 time function, frequency, amplitude, origin time
point source 1 time function processing (frequency filtering)
point source 2 location
point source 2 mechanism
point source 2 time function, frequency, amplitude, origin time
point source 2 time function processing (frequency filtering)
...

number of receivers
receiver 1 location, weight
receiver 2 location, weight
...
```

For example:
```ruby
1								# The unique id of this shot is 1

1								# The shot contains 1 point source
100.0 500.0 10.0				# The x, y, z location of this point source is (100, 500, 10) meters
force 45.0 45.0					# The source is a force vector with (polar, azimuth) = (45, 45) degrees
ricker 20.0 1.0e4 0.0			# The source time function is a Ricker with f0 = 20, A0 = 1e4, and t0 = 0
0.0 0.0							# No processing (frequency filtering) is applied to this source

2								# The shot contains 2 receivers
200.0 400.0 10.0 1.0			# The first receiver is located at (200, 400, 10) meters, with a weight of 1
1200.0 1400.0 100.0 1.0			# The second receiver is located at (1200, 1400, 100) meters, with a weight of 1
```

The source can be an explosion source (i.e., isotropic moment tensor): 
```ruby
	...
	explosion					# The source is an explosion source
	...
```

Or a general moment tensor source:
```ruby
	...
	mt 1.0 1.0 -1.0 -0.5 -0.2 0.1		# The notation convention is (Mxx, Myy, Mzz, Mxy, Mxz, Myz) 
	...
```

For source time function, the valid options are:
- `gaussian`
- `gaussian_deriv` (the first-order derivative of Gaussian)
- `gaussian_deriv_deriv`, `ricker` (the second-order derivative of Gaussian)
- `gaussian_deriv_deriv_deriv`, `ricker_deriv` (the third-order derivative of Gaussian)
- `ormsby` (approximation to sinc with four corner frequencies)
- `custom` (user-provided custom source time function)

To use a custom stf, the user must use the following form:
```ruby
	...
	custom 10.0 1.0 0.0			# The f0 here is used by ADE-CFS-MPML
	custom_wavelet.txt			# Name of the custom stf file, in ASCII format
	... 
```

The custom stf file must be in the following format:
```ruby
	t1 amp1
	t2 amp2
	t3 amp3
	...
	tn ampn
```

That is, each row contains `time amplitude` values for the stf, where the time is in seconds. For example:

```ruby
	0.0			0.0
	1.0e-3		2.0e-5
	2.0e-3		2.5e-5
	....
	1.0e-1		1.0e0
	...
	2.0e-1		0.0
```
The `time` values can be irregularly sampled, although in practice they are usually regularly sampled. The provided stf will be resampled during modeling or FWI to be consistent with the solver. Therefore, the sampling interval can be arbitrary, as long as it is consistent with the physics of the survey. 

For FWI applications, due to historical limitations of SU headers, we do not use SU to store source-receiver geometry information. Therefore, in the provided observed data, the SU files can have null headers (except time-related headers like sampling interval and number of samples), but the number of receivers must be consistent with the geometry file. Future versions of `OWL` may use a more flexible data format for I/O.  

The parameters for defining source-receiver geometry are summarized below: 

#### > Geometry file

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `ns` | integer | Total number of shot gathers in the geometry file | `1` | no |
| `file_geometry` | string | Path to the geometry index file (lists per-shot geometry files) | — | **yes** |

#### > Source/Shot Selection

All source-selection filters are applied after loading the geometry file.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `src_index` | integer list | Shot sequence range: `[first, step, last]` | `[1, 1, ns]` | no |
| `sid_min` | integer | Minimum source field ID to include | `0` | no |
| `sid_max` | integer | Maximum source field ID to include | `INT_MAX` | no |
| `sid_select` | integer list | Whitelist of source field IDs; overrides `src_index` | none | no |
| `sid_exclude` | integer list | Blacklist of source field IDs | none | no |
| `src_select` | integer list | Whitelist of source sequence numbers | none | no |
| `src_exclude` | integer list | Blacklist of source sequence numbers | none | no |
| `sx_min` | float | Minimum source x coordinate | `-∞` | no |
| `sx_max` | float | Maximum source x coordinate | `+∞` | no |
| `sy_min` | float | Minimum source y coordinate | `-∞` | no |
| `sy_max` | float | Maximum source y coordinate | `+∞` | no |
| `sz_min` | float | Minimum source z coordinate | `-∞` | no |
| `sz_max` | float | Maximum source z coordinate | `+∞` | no |

#### > Receiver Selection

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `rec_index` | integer list | Receiver sequence range within each gather: `[first, step, last]` | `[1, 1, INT_MAX]` | no |
| `rec_exclude` | integer list | Receiver sequence numbers to exclude | none | no |
| `rx_min` | float | Minimum receiver x coordinate | `-∞` | no |
| `rx_max` | float | Maximum receiver x coordinate | `+∞` | no |
| `ry_min` | float | Minimum receiver y coordinate | `-∞` | no |
| `ry_max` | float | Maximum receiver y coordinate | `+∞` | no |
| `rz_min` | float | Minimum receiver z coordinate | `-∞` | no |
| `rz_max` | float | Maximum receiver z coordinate | `+∞` | no |
| `offset_min` | float | Minimum source-receiver offset for first-arrival data | `0.0` | no |
| `offset_max` | float | Maximum source-receiver offset for first-arrival data | `+∞` | no |

#### > Geometry Rotation and Source Wavelet

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `geometry_rotate_x` | float | Rotate source and receiver coordinates about the x axis, in degrees | `0.0` | no |
| `geometry_rotate_y` | float | Rotate source and receiver coordinates about the y axis, in degrees | `0.0` | no |
| `geometry_rotate_z` | float | Rotate source and receiver coordinates about the z axis, in degrees | `0.0` | no |
| `geometry_ox` | float | x coordinate of the rotation origin | `0.0` | no |
| `geometry_oy` | float | y coordinate of the rotation origin | `0.0` | no |
| `geometry_oz` | float | z coordinate of the rotation origin | `0.0` | no |
| `f0_factor` | float | Origin-time factor for analytical source wavelets; the origin time is `f0_factor/f0` | `1.0` | no |
| `src_filt_freqs` | float list | Frequencies used to filter analytical source time functions | `[-1.0]` | no |
| `src_filt_coefs` | float list | Filter coefficients paired with `src_filt_freqs` | `[-1.0]` | no |

#### > MPI Decomposition and Time Axis

For 3-D runs, the total number of MPI processes should be

```ruby
ngroup * rankx * ranky * rankz
```

where `rankx`, `ranky`, and `rankz` define the MPI domain decomposition for one shot group, and `ngroup` is the number of shot groups run concurrently. Each MPI process can also use multiple OpenMP threads, for example by setting `OMP_NUM_THREADS` before running `OWL`. For instance, `ngroup = 2`, `rankx = 4`, `ranky = 3`, and `rankz = 2` requires `2*4*3*2 = 48` MPI processes; if `OMP_NUM_THREADS = 8`, each of those 48 processes uses 8 OpenMP threads.

For both 2-D and 3-D runs, the number of MPI ranks used by `mpirun -np` can be smaller than the number of shots `ns`. In that case, the shots are distributed approximately evenly among the available shot groups/ranks; if `ns` is not exactly divisible, some groups receive one more shot than others, and the total number of assigned shots remains `ns`.

For 2-D runs, current `OWL` does not use `rankx`, `ranky`, and `rankz` for domain decomposition. It can use at most `ns` MPI processes, one per shot group, although each process can still use multiple OpenMP threads.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `rankx` | integer | Number of MPI subdomains in x for 3-D runs | `1` | no |
| `ranky` | integer | Number of MPI subdomains in y for 3-D runs | `1` | no |
| `rankz` | integer | Number of MPI subdomains in z for 3-D runs | `1` | no |
| `ngroup` | integer | Number of shot groups run concurrently; for 3-D runs, total MPI processes are `ngroup*rankx*ranky*rankz` | auto | no |
| `dt` | float | Modeling time step | `0.0` | no |
| `tmax` | float | Maximum modeling time; `nt = nint(tmax/dt + 1)` | `1.0` | no |
| `data_dt` | float | Data sampling interval for recorded or synthetic seismograms | `dt` | no |
| `data_tmax` | float | Maximum output data time | `tmax` | no |
| `cc_step_interval` | integer | Wavefield cross-correlation interval for adjoint-state calculations | `1` | no |

### Dimension

`OWL` currently supports (1) regularly sampled grids or (2) curvilinear grids for elastic-wave modeling, FWI, or MT inversion. For convenience, `OWL` supports resampling of the input models from `(nx, ny, nz)` to `(nnx, nny, nnz)`; correspondingly, the grid size can be resampled from `(dx, dy, dz)` to `(ddx, ddy, ddz)`. This is done internally using linear interpolation. For modeling, the outputs (e.g., wavefield snapshots) will have a size of `(nnx, nny, nnz)`; for FWI, the outputs (e.g., updated medium parameter models) will have a size of `(nnx, nny, nnz)`. 

The parameters for defining model dimensions and grid properties are summarized below:

#### > Model Grid (original)

The *original* grid describes the files on disk. A *target* grid (see below) can differ if resampling is needed.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `nx` | integer | Number of grid points in x | — | **yes** |
| `nz` | integer | Number of grid points in z | — | **yes** |
| `ny` | integer | Number of grid points in y (set 1 for 2-D) | `1` | no |
| `dx` | float | Grid spacing in x (same units as coordinates) | — | **yes** |
| `dz` | float | Grid spacing in z | — | **yes** |
| `dy` | float | Grid spacing in y | `1.0` | no |
| `ox` | float | Origin coordinate in x | `0.0` | no |
| `oz` | float | Origin coordinate in z | `0.0` | no |
| `oy` | float | Origin coordinate in y | `0.0` | no |


#### > Model Grid (target, or resampled)

If any target parameter differs from the original, the model is interpolated onto the target grid before use.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `nnx` | integer | Target number of grid points in x | same as `nx` | no |
| `nnz` | integer | Target number of grid points in z | same as `nz` | no |
| `nny` | integer | Target number of grid points in y | same as `ny` | no |
| `ddx` | float | Target grid spacing in x | same as `dx` | no |
| `ddz` | float | Target grid spacing in z | same as `dz` | no |
| `ddy` | float | Target grid spacing in y | same as `dy` | no |
| `oox` | float | Target origin in x | same as `ox` | no |
| `ooz` | float | Target origin in z | same as `oz` | no |
| `ooy` | float | Target origin in y | same as `oy` | no |


#### > Model Domain Restriction

Restrict the active model domain to a sub-region. Points outside are ignored.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `xmin` | float | Minimum x of the active domain | `ox` | no |
| `xmax` | float | Maximum x of the active domain | `ox + (nx-1)*dx` | no |
| `ymin` | float | Minimum y of the active domain | `oy` | no |
| `ymax` | float | Maximum y of the active domain | `oy + (ny-1)*dy` | no |
| `zmin` | float | Minimum z of the active domain | `oz` | no |
| `zmax` | float | Maximum z of the active domain | `oz + (nz-1)*dz` | no |

#### > Adaptive Model Range (per-shot subvolume)

`OWL` provides several parameters to reduce computation by trimming the domain to the source-receiver coverage of each shot.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `yn_adpx` | logical | Enable adaptive x range | `.false.` | no |
| `yn_adpy` | logical | Enable adaptive y range | `.false.` | no |
| `yn_adpz` | logical | Enable adaptive z range | `.false.` | no |
| `adp_extrax` | float | Extra padding beyond source-receiver extent in x | `0.0` | no |
| `adp_extray` | float | Extra padding in y | `0.0` | no |
| `adp_extraz` | float | Extra padding in z | `0.0` | no |


### Medium properties

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `which_medium` | string | Medium type: `acoustic-iso`, `acoustic-tti`, `elastic-iso`, `elastic-vhtiort`, or `elastic-tti` | `acoustic-iso` | no |
| `anisotropy_type` | string | Anisotropy parameterization for anisotropic elastic media: `iso` (isotropic), `vhtiort` (VTI, HTI, or orthorhombic anisotropy), `tti` (TTI), or `cij` (general anisotropy) | `iso` | no |
| `model_name` | string list | Model parameter files required by forward modeling | — | **yes** for forward modeling |
| `file_<name>` | string | Binary file for the model/source parameter `<name>`, where `<name>` is any entry in `model_name`, `model_update`, or `model_aux`; examples include `file_vp`, `file_vs`, `file_rho`, `file_c11`, and `file_mt` | `''` | normally **yes** for each supplied forward or auxiliary model |

In the above parameters, `model_name` is used by forward modeling to list the model files that must be loaded. Each model name is paired with a `file_<name>` parameter. For example, `model_name = vp, rho` expects `file_vp` and `file_rho`. Medium files are raw binary real arrays on the model grid. When model interpolation is enabled, OWL reads the original-grid binary file and interpolates it internally.

The names are literal and case-sensitive in the parameter file. `OWL` can use the following model/source parameter names, depending on `which_medium` and `anisotropy_type`:

| Parameterization | Names |
|------------------|-------|
| Acoustic isotropic | `vp`, optional `rho` |
| Elastic isotropic | `vp`, `vs`, optional `rho` |
| Thomsen anisotropy (`anisotropy_type = thomsen`) | `vp`, `vs`, `epsilon`, `delta`, optional `gamma`, optional tilt angles `theta`, `phi`, optional `rho` |
| Alkhalifah-Tsvankin anisotropy (`anisotropy_type = a-t`) | `vp`, `vs`, `epsilon`, `eta`, optional `gamma`, optional tilt angles `theta`, `phi`, optional `rho` |
| Elastic constants, 2-D VTI/HTI/ORT style (`elastic-vhtiort`) | `c11`, `c13`, `c33`, `c55`, optional `rho` |
| Elastic constants, 2-D TTI/general style (`elastic-tti`) | `c11`, `c13`, `c15`, `c33`, `c35`, `c55`, optional `rho` |
| Elastic constants, 3-D VTI/HTI/ORT style (`elastic-vhtiort`) | `c11`, `c12`, `c13`, `c22`, `c23`, `c33`, `c44`, `c55`, `c66`, optional `rho` |
| Elastic constants, 3-D TTI/general style (`elastic-tti`) | `c11`, `c12`, `c13`, `c14`, `c15`, `c16`, `c22`, `c23`, `c24`, `c25`, `c26`, `c33`, `c34`, `c35`, `c36`, `c44`, `c45`, `c46`, `c55`, `c56`, `c66`, optional `rho` |

### Data and data processing

`OWL` reads and writes one SU file per shot and component, using names such as `shot_1_seismogram_z.su`. The active components are controlled by `data_name`; acoustic media default to pressure (`p`), 2-D elastic media to `x z`, and 3-D elastic media to `x y z`.

#### > Data Components and Snapshots

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `data_name` | string list | Data components to model or invert: `x`, `y`, `z`, and/or `p` | medium dependent | no |
| `snaps` | float list | Wavefield snapshot schedule in the form `snapshot_begin, snapshot_interval, snapshot_end`; intermediate snapshot times are computed automatically. `[-1.0]` disables requested snapshots | `[-1.0]` | no |

#### > Processing Lists

Processing is specified as an ordered list of step names. Each step is applied to every selected component and shot.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `process_record` | string list | Processing applied once to observed records before FWI iterations | `['']` | no |
| `process_record_encoded` | string list | Processing applied to encoded observed records | `['']` | no |
| `process_synthetic` | string list | Processing applied to synthetic seismograms; available in forward modeling and FWI | `['']` | no |
| `process_adjsrc` | string list | Processing applied to adjoint sources | `['']` | no |

Available data-processing steps are:

| Step | Description |
|------|-------------|
| `stf_correction` | Estimate a matching filter between synthetic and observed records and convolve the synthetic data |
| `time_shift` | Shift traces by `dp_t_shift` |
| `time_deriv` | Differentiate traces in time |
| `time_integ` | Integrate traces in time |
| `top_mute` | Zero samples before a top-mute time and taper the kept data |
| `bottom_mute` | Zero samples after a bottom-mute time and taper the kept data |
| `surgical_mute` | Mute or keep a time window between two moveout curves |
| `subtract_base` | Subtract matching SU files from `dir_base` |
| `t_power` | Multiply samples by `t**dp_t_power` |
| `t_balance` | Apply moving-window RMS/AGC balancing |
| `freq_filt` | Apply frequency-domain filtering |
| `rms_balance` | Normalize each trace by its L2 norm |
| `offset_balance` | Scale traces by decreasing source-receiver offset weight |
| `max_balance` | Normalize each trace by its peak absolute amplitude |
| `dip_filt` | Apply dip-domain filtering |
| `andf_filt` | Apply anisotropic-diffusion filtering |
| `remove_nan` | Replace NaN or Inf samples with zero |
| `rotate` | Rotate three-component data using `dp_rotate_x`, `dp_rotate_y`, and `dp_rotate_z` |

#### > Mute Parameters

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `dp_dir_top_mute_time` | string | Directory containing per-shot `shot_<sid>_traveltime_p.bin` files for top mute | `''` | no |
| `dp_top_mute_vel` | float | Constant velocity used to compute top-mute moveout from offset | `-1.0` | no |
| `dp_top_mute_width` | float | Time width added to the top-mute moveout | `0.1*tmax` | no |
| `dp_top_mute_taper` | float | Hann taper length after top mute | `0.1*tmax` | no |
| `dp_top_mute_shift` | float | Additional top-mute time shift | `0.0` | no |
| `dp_top_mute_const` | float | If positive, use a constant top-mute time instead of moveout | `-tmax` | no |
| `dp_dir_bottom_mute_time` | string | Directory containing bottom-mute traveltime files | `''` | no |
| `dp_bottom_mute_vel` | float | Constant velocity used to compute bottom-mute moveout; must be smaller than `dp_top_mute_vel` when both are active | `-2.0` | no |
| `dp_bottom_mute_width` | float | Reserved bottom-mute width parameter | `0.1*tmax` | no |
| `dp_bottom_mute_taper` | float | Hann taper length before bottom mute | `0.1*tmax` | no |
| `dp_bottom_mute_from_top` | float | Build bottom mute by shifting the top-mute time by this amount | `0.0` | no |
| `dp_bottom_mute_const` | float | If positive, use a constant bottom-mute time instead of moveout | `-tmax` | no |
| `dp_dir_surgical_mute_time` | string | Directory containing two-window traveltime files, with `2*nr` floats per shot | `''` | no |
| `dp_surgical_mute_vel` | float list | Two velocities defining surgical-mute moveout when no traveltime file is provided | `[-1.0, -1.0]` | no |
| `dp_surgical_mute_width` | float | Time width added to the second surgical-mute curve | `0.1*tmax` | no |
| `dp_surgical_mute_taper` | float list | Two taper lengths for the surgical mute | `[0.1*tmax, 0.1*tmax]` | no |
| `dp_surgical_mute_inverse` | logical | Keep the protected surgical-mute window instead of muting it | `.false.` | no |

#### > Filter, Balance, and Rotation Parameters

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `dp_t_shift` | float | Time shift used by `time_shift` | `0.0` | no |
| `dp_t_power` | float | Exponent used by `t_power` | `0.0` | no |
| `dp_t_balance_window` | float | Moving time window for `t_balance` | `tmax/20` | no |
| `dp_freq_filt_freqs` | float list | Frequencies for `freq_filt` | `[-1.0]` | no |
| `dp_freq_filt_coefs` | float list | Filter coefficients paired with `dp_freq_filt_freqs` | `[-1.0]` | no |
| `dp_dip_filt_dips` | float list | Dip samples for `dip_filt` | `[-1000.0, 0.0, 1000.0]` | no |
| `dp_dip_filt_coefs` | float list | Dip-filter coefficients | `[0.0, 1.0, 0.0]` | no |
| `dp_andf_smoothx` | float | Trace-direction smoothing for `andf_filt` | `2.0` | no |
| `dp_andf_smootht` | float | Time-direction smoothing for `andf_filt` | `tmax/20` | no |
| `dp_andf_powerm` | float | ANDF power-law exponent | `1.0` | no |
| `dp_andf_t` | integer | Number of ANDF iterations | `10` | no |
| `dp_andf_sigma` | float | ANDF diffusion threshold | `10.0` | no |
| `dp_andf_alpha` | float | ANDF alpha/lambda1 parameter | `0.0` | no |
| `dp_andf_beta` | float | ANDF beta/lambda2 parameter | `1.0` | no |
| `dp_rotate_x` | float | x-axis rotation angle in degrees for the `rotate` processing step | `0.0` | no |
| `dp_rotate_y` | float | y-axis rotation angle in degrees for the `rotate` processing step | `0.0` | no |
| `dp_rotate_z` | float | z-axis rotation angle in degrees for the `rotate` processing step | `0.0` | no |

### Free surface and topography

For elastic modeling, FWI, and MT inversion, `OWL` can use a free-surface boundary condition for the top surface. The free surface can be topographic, with the topography specified by an external file. 

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `yn_free_surface` | logical | Whether to assume free surface boundary condition for the top surface | `.false.` | no |
| `free_surface_dz_refine` | float | Near-surface vertical mesh-refinement factor | `4.0` | no |
| `file_topo` | string | Topography file for a topographic free surface | `''` | no |
| `topo_interp` | string | Interpolation method for the topography file | `cubic` | no |
| `measure_source_depth_from_surface` | logical | Treat source z values as depths below the free surface | `.false.` | no |
| `measure_receiver_depth_from_surface` | logical | Treat receiver z values as depths below the free surface | `.false.` | no |
| `source_vertical_to_surface` | logical | Project source depth vertically to the local surface | `.false.` | no |
| `receiver_vertical_to_surface` | logical | Project receiver depth vertically to the local surface | `.false.` | no |
| `yn_save_mesh` | logical | Save generated free-surface or topographic mesh files | `.false.` | no |


### Inversion

#### > Inversion Control

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `niter_max` | integer | Maximum number of inversion iterations | `100` | no |
| `misfit_type` | string | Misfit formulation. Available values are listed in the misfit section below | `waveform` | no |
| `yn_continue` | logical | Resume an interrupted inversion from the last completed iteration | `.false.` | no |
| `resume_from_iter` | integer | Iteration number to resume from (if `yn_continue = .false.`) | `1` | no |
| `yn_flat_stop` | logical | Stop when three consecutive iterations have identical misfit | `.false.` | no |
| `yn_shared_model_processing` | logical | Use shared model-processing settings across updated parameters | `.true.` | no |
| `yn_energy_precond` | logical | Apply energy-based gradient preconditioning | `.false.` | no |

#### > Medium properties 

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `model_update` | string list | Model/source parameters to invert for, such as `vp`, `vs`, `rho`, anisotropy parameters, `cij`, or `mt` | `vp` | no |
| `model_aux` | string list | Auxiliary (fixed) model parameter names | none | no |
| `file_<name>` | string | Binary file for the model/source parameter `<name>`, where `<name>` is any entry in `model_name`, `model_update`, or `model_aux`; examples include `file_vp`, `file_vs`, `file_rho`, `file_c11`, and `file_mt` | `''` | normally **yes** for each supplied forward or auxiliary model |
| `min_vpvsratio` | float | Minimum allowed Vp/Vs ratio (elastic media) | `1.1` | no |
| `max_vpvsratio` | float | Maximum allowed Vp/Vs ratio | `9.0` | no |
| `vpvsratio_smoothx` | float | Gaussian smoothing length of Vp/Vs ratio in x (0 = off) | `0.0` | no |
| `vpvsratio_smoothy` | float | Gaussian smoothing length of Vp/Vs ratio in y | `0.0` | no |
| `vpvsratio_smoothz` | float | Gaussian smoothing length of Vp/Vs ratio in z | `0.0` | no |

In the above parameters, `model_update` lists the parameters to update, while `model_aux` lists additional fixed parameters that are needed by the chosen parameterization but should not be updated. For example, an elastic inversion for `vp` only may still need `vs` and `rho` in `model_aux` so the solver can build the full medium.

Each model/source name is paired with a `file_<name>` parameter. For example, `model_update = vp, vs` with `model_aux = rho` expects initial/update files such as `file_vp` and `file_vs`, plus the fixed auxiliary file `file_rho`. Medium files are raw binary real arrays on the model grid; source files such as `file_mt` use source-indexed arrays. When model interpolation is enabled, OWL reads the original-grid binary file and interpolates it internally.

Additionally, for the moment tensor, the component order follows the geometry-file convention: `Mxx`, `Myy`, `Mzz`, `Mxy`, `Mxz`, `Myz`. The parameterization for this inversion is:

| Parameterization | Names |
|------------------|-------|
| Source inversion | `mt` for six moment-tensor components per shot; `stf` is recognized as a source-update name in the inversion bookkeeping |

Parameters with defaults, such as `rho`, are loaded from `model_update` or `model_aux` when provided; otherwise the code may use the default shown above for that solver path. 

#### > Misfit Parameters

`OWL` computes adjoint sources from the selected `misfit_type`. Before computing the misfit, traces are skipped when the receiver weight is zero, the observed or synthetic trace norm is zero, or either trace norm is below `trace_discard_threshold` times the maximum trace norm for that component.

Available misfit types are:

| `misfit_type` | Description | Main Parameters |
|----------------|-------------|-----------------|
| `waveform` | L2 waveform-difference misfit, `sum((d_obs - d_syn)**2)` | `adj_nt`, `trace_discard_threshold` |
| `corr` | Zero-lag normalized correlation misfit | `adj_nt`, `trace_discard_threshold` |
| `envelope` | Envelope-difference misfit | `adj_nt`, `trace_discard_threshold` |
| `phase` | Instantaneous phase misfit | `adj_nt`, `trace_discard_threshold` |
| `adaptive` | Adaptive waveform inversion (AWI), implemented as weighted normalized deconvolution | `tlag_max`, `penalty_method`, `penalty_power`, `deconv_eps`, `adj_nt`, `trace_discard_threshold` |
| `local-adaptive` | Localized AWI (LAWI), using local/Gabor-windowed matching filters | `tlag_max`, `penalty_method`, `penalty_power`, `deconv_eps`, `lawi_sigma`, `adj_nt`, `trace_discard_threshold` |
| `adaptive-spacetime` | Space-time adaptive misfit using a local trace window around each receiver | `adaptive_half_window`, `tlag_max`, `penalty_method`, `penalty_power`, `deconv_eps`, `adj_nt`, `trace_discard_threshold` |
| `dtw` | Generalized dynamic-time-warping misfit, computed component by component for each shot gather | all `dtw_*` parameters, `tlag_max`, `adj_nt`, `trace_discard_threshold` |
| `dtw-vector` | Vector generalized dynamic-time-warping misfit, computed with all components together | all `dtw_*` parameters, `tlag_max`, `adj_nt`, `trace_discard_threshold` |

Unknown `misfit_type` values fall back to `waveform` inside the component-by-component branch.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `tlag_max` | float | Maximum time lag for lag-based misfits | `0.1*tmax` | no |
| `adj_nt` | integer | Number of samples used internally for adjoint-source/misfit calculations; if positive and different from the data length, traces are resampled for the calculation and the adjoint source is resampled back | `nt` | no |
| `envelope_power` | float | Envelope exponent read by the parameter module; currently not used by the adjoint-source routines listed here | `2.0` | no |
| `penalty_method` | string | Penalty function for AWI/LAWI/adaptive-space-time matching-filter lags: `linear`, `power`, `gaussian`, or `exp` | `linear` | no |
| `penalty_power` | float | Power used by `power`, `gaussian`, and `exp` lag penalties | `1.0` | no |
| `deconv_eps` | float | Stabilization epsilon for deconvolution/AWI-style misfits | `0.1` | no |
| `lawi_sigma` | float | Gaussian local-window sigma for `local-adaptive`, in seconds | `0.25` | no |
| `adaptive_half_window` | integer | Half-window, in receiver traces, used by `adaptive-spacetime` | `3` | no |
| `yn_average_misfit` | logical | Read by the parameter module; not applied in `module_inversion_adjoint_source.f90` | `.false.` | no |
| `trace_discard_threshold` | float | Discard extremely small-amplitude traces below this threshold | `1.0e-6` | no |

For DTW-based misfit, there are several dedicated parameters. The DTW misfits estimate a time-shift field `tau` that warps observed data toward synthetic data. For `misfit_type = dtw`, each component is processed independently. For `misfit_type = dtw-vector`, all selected components are processed together with a shared time-shift field. When adjoint sources are saved, `dtw` also writes `shot_<sid>_time_shift_<component>.su`; `dtw-vector` writes `shot_<sid>_time_shift.su`.

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `dtw_niter` | integer | Maximum number of DTW solver iterations | `5` | no |
| `dtw_rinst` | float | Instantaneous regularization radius/weight | `0.1` | no |
| `dtw_rcuml` | float | Cumulative regularization radius/weight | `1.0` | no |
| `dtw_epsabs` | float | Absolute convergence tolerance | `1.0e-2` | no |
| `dtw_epsrel` | float | Relative convergence tolerance | `1.0e-3` | no |
| `dtw_loss` | string | DTW loss norm | `l0.5` | no |
| `dtw_form` | string | DTW adjoint/objective form: `phase`, `amp`, or `phase+amp` | `phase` | no |
| `dtw_amp_weight` | float | Amplitude-term weight for `dtw_form = phase+amp` | `1.0` | no |
| `dtw_trc` | integer | Number of neighboring traces on each side included in each DTW solve; `0` uses only the current trace | `0` | no |
| `dtw_smooth_median` | float | Median smoothing width, in trace samples, for DTW time shifts and weights | `0.0` | no |
| `dtw_smooth_gaussian` | float | Gaussian smoothing width, in trace samples, for DTW time shifts and weights | `0.0` | no |

#### > Model Bounds and Step Limits

For each model parameter `<name>`:

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `min_<name>` | float | Lower bound for model values during update | `0.0` (velocity), `0.0` (anisotropy), `0.0` (elastic constant) | no |
| `max_<name>` | float | Upper bound for model values | `1e5` (velocity), `0.5` (anisotropy), `1e9` (elastic constant) | no |
| `step_max_<name>` | float | Maximum absolute perturbation allowed for one update of parameter `<name>`; used to set the initial model step and to reject oversized trial steps | `100.0` for `vp`, `vs`, `rho`; `0.1` for `epsilon`, `delta`, `gamma`, `eta`; `0.1*pi/2` for `theta`, `phi`; `1.0e9` for `c*`/`C*`; `1.0` for `mt` | no |

#### > Search Direction

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `search_method` | string | Global inversion search direction. Supported values are `SD`/`sd`/`steepest-descent`, `CG`/`cg`/`conjugate-gradient`, and `L-BFGS`/`l-bfgs`/`l-BFGS` | `cg` | no |
| `search_method_<name>` | string | Per-parameter override of `search_method`; supports the same values | same as `search_method` | no |

#### > Step Size

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `step_size_method` | string | Step-size computation method: `linear`, `quadratic`, or `line-search` | `linear` | no |
| `jumpout_factor` | float | Misfit acceptance relaxation factor. Trial steps are accepted when the trial misfit is below `jumpout_factor` times the current misfit | `1.0`; `1.05` when `trigger_jumpout` is active | no |
| `yn_enforce_update` | logical | For `linear` step-size search, force one accepted update trial even if the misfit increases | `.false.` | no |
| `step_max_<name>` | float | Per-parameter maximum perturbation used by all step-size methods | See model-bounds table above | no |

The model update is `m = clip(m + model_step*search_direction*step_scaling_factor, min_<name>, max_<name>)`, followed by Vp/Vs-ratio clipping where applicable. The initial `model_step` is computed from `step_max_<name>/max(abs(search_direction))`; trial `step_scaling_factor` values are halved until every parameter perturbation is within its `step_max_<name>` range.

Available step-size methods are:

| Method | Description |
|--------|-------------|
| `linear` | Perturbation-based step estimate. OWL first computes an initial trial step (`0.2`, halved until suitable), runs one trial forward/misfit calculation to estimate linearized coefficients, then backtracks by halving until the trial misfit is acceptable. If no acceptable trial is found within the internal search limit, the step is set to zero. |
| `quadratic` | Three-point quadratic fit. OWL evaluates the current model, `0.1*step_scaling_factor`, and `step_scaling_factor`, then fits a quadratic step and evaluates that step. |
| `line-search` | Quadratic-plus-bisection line search. OWL brackets the search between zero and the largest suitable initial step, tries the right endpoint and midpoint, then alternates quadratic interpolation and bisection until the internal search limit or a small interval is reached. If the best trial misfit remains above `jumpout_factor` times the current misfit, the step is set to zero. |

All three methods save `updated_<name>.bin` files in the iteration model directory after applying the final step. Trial synthetics may be written in scratch and moved into the iteration synthetic directory when the trial is accepted.

#### > Regularization

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `model_regularization_method` | string list | Regularization methods applied to model updates; empty disables model regularization | `['']` | no |
| `source_regularization_method` | string list | Regularization methods applied to source updates; empty disables source regularization | `['']` | no |

Available model regularization methods are applied in the order listed:

| Method | Associated Parameters | Description |
|--------|-----------------------|-------------|
| `tikhonov`, `Tikhonov` | `reg_tikhonov_lambda` (`10.0`) | Tikhonov denoising/damping |
| `smooth` | `reg_smoothx`, `reg_smoothy` (3-D), `reg_smoothz`; or per-parameter `reg_smoothx_<name>`, `reg_smoothy_<name>` (3-D), `reg_smoothz_<name>` when the global value is negative | Gaussian smoothing regularization. Defaults are `-1.0` for global values and one grid spacing for per-parameter values |
| `tgpv`, `TGpV` | `reg_tv_mu_<name>`, `reg_tv_lambda1`, `reg_tv_lambda2`, `reg_tv_norm`, `reg_tv_niter`; 2-D also supports per-parameter `reg_tv_lambda1_<name>` and `reg_tv_lambda2_<name>` when global values are `-1.0` | Total generalized p-variation/TV-style denoising. Defaults: `reg_tv_norm = 0.5`, `reg_tv_niter = 50`, `reg_tv_lambda1 = reg_tv_lambda2 = 1.0` in 3-D and per-parameter fallback in 2-D |
| `structure` | `reg_andf_alpha`, `reg_andf_beta`, `reg_andf_gamma` (3-D), `reg_andf_smoothx`, `reg_andf_smoothy` (3-D), `reg_andf_smoothz`, `reg_andf_t`, `reg_andf_sigma`, `reg_andf_powerm`, `reg_andf_aux`, `reg_andf_coh` (2-D) | Structure-oriented anisotropic-diffusion regularization. Defaults are `alpha = 0.001`, `beta = 1.0`, `gamma = 1.0` in 3-D, `smoothx = 2.0`, `smoothy = 2.0` in 3-D, `smoothz = 8.0`, `sigma = 10.0`, `powerm = 4.0`; `reg_andf_t` defaults to `10` in 2-D and `5` in 3-D |

Many regularization parameters are read with iteration-aware readers, so values can be changed by iteration when supported by the parameter-file syntax.

#### > Gradient and Model Processing

Processing is specified as a list of steps applied in order.

#### > Per-Shot Gradient Processing

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `process_shot_grad` | string list | Processing steps applied to each shot's gradient before stacking; internally read as `process_shot_grad` | `['']` (none) | no |

The processing parameter name is `process_shot_<name>`, where `<name>` is commonly `grad` for `process_shot_grad`. Per-shot processing uses shot-local grid spacing (`mdx`, `mdy`, `mdz`) and crops masks or auxiliary fields to the active shot domain.

Available per-shot gradient-processing steps:

| Step | 2-D Parameters | 3-D Parameters | Description |
|------|----------------|----------------|-------------|
| `smooth` | `shot_<name>_smooth_x`, `shot_<name>_smooth_z` (`3*mdx`, `3*mdz`) | `shot_<name>_smoothx`, `shot_<name>_smoothy`, `shot_<name>_smoothz` (`3*mdx`, `3*mdy`, `3*mdz`) | Gaussian smoothing |
| `max_balance` | none | none | Normalize by the maximum value |
| `rms_balance` | none | none | Normalize by mean/RMS-style image energy |
| `moving_balance` | `shot_<name>_moving_balance_x`, `shot_<name>_moving_balance_z` (`3*mdx`, `3*mdz`) | `shot_<name>_movingbalx`, `shot_<name>_movingbaly`, `shot_<name>_movingbalz` (`6*mdx`, `6*mdy`, `6*mdz`) | Moving-window amplitude balancing |
| `median_filt` | `shot_<name>_median_filt_x`, `shot_<name>_median_filt_z` (`mdx`, `mdz`) | `shot_<name>_medianfiltx`, `shot_<name>_medianfilty`, `shot_<name>_medianfiltz` (`mdx`, `mdy`, `mdz`) | Median filtering |
| `dip_filt` | `shot_<name>_dip_filt_zx`, `shot_<name>_dip_filt_zx_coefs` | `shot_<name>_dipfiltzx`, `shot_<name>_dipfiltzx_amps`, `shot_<name>_dipfiltzy`, `shot_<name>_dipfiltzy_amps`, `shot_<name>_dipfiltyx`, `shot_<name>_dipfiltyx_amps` | Dip-domain filtering |
| `remove_nan` | none | none | Replace NaN/Inf values with finite values |
| `laplace_filt` | not available | none | Apply a Laplacian filter |
| `andf_filt` | `shot_<name>_andf_smooth_x`, `shot_<name>_andf_smooth_z`, `shot_<name>_andf_powerm`, `shot_<name>_andf_t`, `shot_<name>_andf_sigma`, `shot_<name>_andf_alpha`, `shot_<name>_andf_beta`, `shot_<name>_andf_aux`, `shot_<name>_andf_coh` | `shot_<name>_andf_smoothx`, `shot_<name>_andf_smoothy`, `shot_<name>_andf_smoothz`, `shot_<name>_andf_powerm`, `shot_<name>_andf_t`, `shot_<name>_andf_sigma`, `shot_<name>_andf_alpha`, `shot_<name>_andf_beta`, `shot_<name>_andf_gamma`, `shot_<name>_andf_aux`, `shot_<name>_andf_coh` | Structure-oriented anisotropic-diffusion filtering |
| `wavenumber_filt` | `shot_<name>_wavenumber_filt_x`, `shot_<name>_wavenumber_filt_x_coefs`, `shot_<name>_wavenumber_filt_z`, `shot_<name>_wavenumber_filt_z_coefs` | `shot_<name>_wavenumx`, `shot_<name>_wavenumx_amps`, `shot_<name>_wavenumy`, `shot_<name>_wavenumy_amps`, `shot_<name>_wavenumz`, `shot_<name>_wavenumz_amps` | Wavenumber-domain filtering |
| `taper` | `shot_<name>_taper_x`, `shot_<name>_taper_z` (`[0.0, 0.0]`) | `shot_<name>_taperx`, `shot_<name>_tapery`, `shot_<name>_taperz` (`[0.0, 0.0]`) | Blackman taper; one value is expanded to both sides |
| `mask` | `dir_shot_<name>_mask` or `shot_<name>_mask` | same | Multiply by a mask. Directory masks are read as `<dir_scratch>/shot_<sid>_mask.bin` under the supplied directory |
| `adaptive_mute` | `shot_<name>_adaptive_mute_x`, `shot_<name>_adaptive_mute_z` | not available | Taper outside the source-receiver aperture in x/z |
| `cone_mute` | `shot_<name>_cone_mute_x`, `shot_<name>_cone_mute_z`, `shot_<name>_cone_mute_power`, `shot_<name>_cone_mute_taper` | not available | Cone-shaped mute below the source |

#### > Global Gradient Processing

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `process_grad` | string list | Processing steps applied to the stacked gradient | `['']` (none) | no |
| `process_srch` | string list | Processing steps applied to the search direction | `['']` (none) | no |

The processing parameter name is `process_<name>`, where `<name>` is commonly `grad` or `srch`. If the routine is called for a model parameter, `<param_name>_update_iter` can restrict updates to an iteration range; a single value means `[value, niter_max]`.

Available global gradient/search-direction processing steps:

| Step | 2-D Parameters | 3-D Parameters | Description |
|------|----------------|----------------|-------------|
| `scale` | `<name>_scale` (`1.0`) | same | Multiply by a scalar |
| `max_balance` | none | not available | Normalize by maximum value |
| `rms_balance` | none | none | Normalize by mean/RMS-style energy |
| `rms_balance_x` | `<name>_rms_balance_x` (`dx`) | `<name>_rmsbalx` (`dx`) | Sliding RMS normalization along x |
| `rms_balance_y` | not available | `<name>_rmsbaly` (`dy`) | Sliding RMS normalization along y |
| `rms_balance_xy` | not available | `<name>_rmsbalx`, `<name>_rmsbaly` (`dx`, `dy`) | Sliding RMS normalization in x-y planes |
| `rms_balance_z` | `<name>_rms_balance_z` (`dz`) | `<name>_rmsbalz` (`dz`) | Sliding RMS normalization along z |
| `moving_balance` | `<name>_moving_balance_x`, `<name>_moving_balance_z` (`6*dx`, `6*dz`) | `<name>_movingbalx`, `<name>_movingbaly`, `<name>_movingbalz` (`3*dx`, `3*dy`, `3*dz`) | Moving-window amplitude balancing |
| `taper` | `<name>_taper_x`, `<name>_taper_z` (`[0.0, 0.0]`) | `<name>_taperx`, `<name>_tapery`, `<name>_taperz` (`[0.0, 0.0]`) | Blackman taper |
| `smooth` | `<name>_smooth_x`, `<name>_smooth_z` (`3*dx`, `3*dz`) | `<name>_smoothx`, `<name>_smoothy`, `<name>_smoothz` (`3*dx`, `3*dy`, `3*dz`) | Gaussian smoothing |
| `andf_filt` | `<name>_andf_smooth_x`, `<name>_andf_smooth_z`, `<name>_andf_powerm`, `<name>_andf_t`, `<name>_andf_sigma`, `<name>_andf_alpha`, `<name>_andf_beta`, `<name>_andf_aux`, `<name>_andf_coh` | `<name>_andf_smoothx`, `<name>_andf_smoothy`, `<name>_andf_smoothz`, `<name>_andf_powerm`, `<name>_andf_t`, `<name>_andf_sigma`, `<name>_andf_alpha`, `<name>_andf_beta`, `<name>_andf_gamma`, `<name>_andf_aux`, `<name>_andf_coh`, `<name>_andf_rankx`, `<name>_andf_ranky`, `<name>_andf_rankz` | Structure-oriented anisotropic-diffusion filtering; 3-D uses MPI-aware filtering |
| `median_filt` | `<name>_median_filt_x`, `<name>_median_filt_z` (`dx`, `dz`) | `<name>_medianfiltx`, `<name>_medianfilty`, `<name>_medianfiltz` (`dx`, `dy`, `dz`) | Median filtering |
| `mask` | `<name>_mask` | same | Multiply by a model-sized mask |

Many global processing parameters are read with iteration-aware readers, so values can be changed by iteration when supported by the parameter-file syntax.


### Other parameters

`OWL` uses the following parameters to create the relevant directories: 

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `dir_working` | string | Working directory for FWI output and processed observed data | `./test` | no |
| `dir_record` | string | Directory containing observed SU records | `./data` | no |
| `dir_synthetic` | string | Output directory for synthetic SU data | `./data_synthetic` | no |
| `dir_base` | string | Directory containing baseline SU data for `subtract_base` processing | `./data_base` | no |
| `dir_processed` | string | Output directory for processed synthetic data | `./data_processed` | no |
| `dir_snapshot` | string | Directory for wavefield snapshots | `./snapshot` | no |
| `verbose` | logical | Print extra progress and diagnostic information | `.false.` | no |

`OWL` saves misfits to the following files: 

| Parameter | Type | Description | Default | Required |
|-----------|------|-------------|---------|----------|
| `file_data_misfit` | string | ASCII file to record per-iteration data misfit | `<dir_working>/data_misfit.txt` | no |
| `file_shot_misfit` | string | Binary file to record per-shot misfit at each iteration; the dimension of this file is `N1*N2 = ns*niter` where `ns` is the number of shots and `niter` is the number of iterations already accomplished plus 1  | `<dir_working>/shot_misfit.bin` | no |



## Examples

We provide several reproducible examples in the `example` directory. Each example contains a parameter file and a small driver script that prepares the input files, runs `OWL`, and, when useful, plots the results.

As a simple starting point, consider `example/forward_2d/acoustic`. Its parameter file, `param_modeling.rb`, is:

```ruby
nx = 200
nz = 100

dx = 10
dz = 10

dt = 1.5e-3
data_dt = 1.0e-3
tmax = 2

ns = 1
file_geometry = ./geometry/geometry.txt

model_name = vp, rho
file_vp = model/vp.bin
file_rho = model/rho.bin

snaps = 0.0, 0.1, 1.5

verbose = y
```

This file defines a 2-D acoustic modeling run. The grid has `nx = 200` samples in x and `nz = 100` samples in z, with 10 m spacing in both directions. The time axis is controlled by `dt`, `data_dt`, and `tmax`: the solver advances with a 1.5 ms internal time step, saves seismograms at 1.0 ms sampling, and models 2 s of wave propagation. The geometry is read from `./geometry/geometry.txt`, which is created by the example setup program. The model is acoustic isotropic by default, so the listed model parameters are `vp` and `rho`, read from `model/vp.bin` and `model/rho.bin`. The `snaps` syntax is `snapbeg, snapinterval, snapend`; for example, `snaps = 0.0, 0.1, 1.5` requests snapshots from 0.0 s to 1.5 s every 0.1 s, and the intermediate snapshot times are computed automatically.

The example driver `test.rb` shows the usual workflow:

```ruby
system "x_runf90 create_test.f90"
system "mpirun -np 1 owl_modeling2 param_modeling.rb"
```

The first command builds the model and geometry files for the example. The second command runs the 2-D forward-modeling executable with `param_modeling.rb` as the parameter file. In other examples the executable and parameter file may change, for example `owl_inversion2 param_fwi.rb` for a 2-D inversion, but the idea is the same: the last command-line argument is the file from which `OWL` reads its run settings.

For an FWI example, see `example/inverse_2d/topo`. This example first generates synthetic observed data with `param_modeling.rb`, then runs several inversions from the shared base file `param_fwi.rb`. The core inversion settings are:

```ruby
nx = 301
nz = 101

dx = 10
dz = 10

dt = 1.25e-4
data_dt = 1.0e-3
tmax = 2.5

ns = 60
file_geometry = ./geometry/geometry.txt

which_medium = elastic-tti
anisotropy_type = iso
model_update = vp, vs
model_aux = rho
file_vp = model/vp_init.bin
file_vs = model/vs_init.bin
file_rho = model/rho.bin

min_vp = 2600
max_vp = 4800
min_vs = 1500
max_vs = 2800

step_max_vp = 100
step_max_vs = 50

yn_free_surface = y
measure_source_depth_from_surface = y
measure_receiver_depth_from_surface = y
free_surface_dz_refine = 2
file_topo = ftopo.txt

process_grad = smooth, mask
grad_smooth_x = 20
grad_smooth_z = 10
grad_mask = model/mask.bin

dir_record = data
niter_max = 200
yn_energy_precond = y
verbose = y
```

The first block again defines the computational grid, time axis, shot count, and geometry. The medium is run with the elastic TTI solver, but `anisotropy_type = iso` selects an isotropic parameterization inside that solver. The inversion updates `vp` and `vs`, while `rho` is listed as an auxiliary model: it is needed to build the elastic medium but is not updated. The files `model/vp_init.bin` and `model/vs_init.bin` are therefore the starting models for FWI, while `model/rho.bin` remains fixed.

The bounds and step limits control the model update. `min_vp`, `max_vp`, `min_vs`, and `max_vs` clip the updated models into physically reasonable ranges. `step_max_vp` and `step_max_vs` limit the maximum perturbation applied to each parameter during one update. The free-surface block enables topography, interprets source and receiver depths relative to the surface, refines the near-surface mesh in z, and reads the topographic surface from `ftopo.txt`.

The gradient-processing block regularizes the inversion before the model update. Here `process_grad = smooth, mask` applies Gaussian smoothing to the gradient and then multiplies it by `model/mask.bin`, so updates can be tapered or suppressed in selected regions. `dir_record = data` tells FWI where to find the observed records, `niter_max = 200` sets the maximum number of inversion iterations, and `yn_energy_precond = y` enables source-energy preconditioning.

The driver script shows a practical way to reuse the same base FWI settings for several objective functions:

```bash
p=param_waveform.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = waveform ' >> $p
echo 'dir_working = test_waveform ' >> $p
mpirun -np 60 owl_fwi2 $p

p=param_adaptive.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = adaptive ' >> $p
echo 'dir_working = test_adaptive ' >> $p
echo 'tlag_max = 1.0' >> $p
echo 'deconv_eps = 0.2' >> $p
mpirun -np 60 owl_fwi2 $p
```

In this pattern, `param_fwi.rb` holds the common grid, model, topography, and update controls. The script copies it to a run-specific file, appends the misfit parameters, and runs `owl_fwi2`. For waveform inversion, only `misfit_type` and `dir_working` are added. For adaptive waveform inversion, the script also sets `tlag_max` and `deconv_eps`. The same idea is used for the DTW-vector run, where parameters such as `adj_nt`, `dtw_form`, `dtw_smooth_median`, `dtw_smooth_gaussian`, `dtw_rinst`, `dtw_rcuml`, `dtw_loss`, and `jumpout_factor` tune the DTW adjoint source and line-search acceptance.

Parameter files use FLIT-style parameter reading. In practice, a parameter is written as

```ruby
parameter_name = value
```

Scalar values are written directly, lists are written as comma-separated values, and file names are written as plain paths. For example, `model_name = vp, rho` tells `OWL` to read two model parameters, while `snaps = 0.0, 0.1, 1.5` gives a snapshot schedule of begin time, interval, and end time. Parameters can be ordered by topic for readability; the code reads them by name and uses defaults for optional entries that are not present.

One useful FLIT convention is that reading stops at a line containing `exit`. Any parameters written after `exit` are ignored by `OWL`. This is convenient for keeping notes, old settings, or alternative parameter blocks in the same file without affecting the current run:

```ruby
nx = 200
nz = 100
dt = 1.5e-3

exit

# These lines are kept for reference only and will not be read.
nx = 400
dt = 0.75e-3
```

When using a driver script that appends parameters with `echo ... >> $p`, place `exit` only after all active settings. If `exit` appears in the copied base file before appended settings such as `misfit_type`, `dir_working`, or `tlag_max`, those appended settings will not be read.

When modifying an example, it is usually safest to start from the supplied parameter file, change one group of settings at a time, and keep the model dimensions, binary model files, and geometry mutually consistent.
