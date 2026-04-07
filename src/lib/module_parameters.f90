!
! © 2025. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!


module mod_parameters

    use libflit

    implicit none

    !==========================================================================
    ! General

    ! Parameter file
    character(len=1024) :: file_parameter

    ! Paths
    character(len=1024) :: dir_record
    character(len=1024) :: dir_working
    character(len=1024) :: dir_scratch
    character(len=1024) :: dir_synthetic
    character(len=1024) :: dir_synthetic_processed
    character(len=1024) :: dir_snapshot
    character(len=1024) :: dir_base
    character(len=1024) :: dir_muted
    character(len=1024) :: dir_adjoint
    character(len=1024) :: dir_base_encoded
    character(len=1024) :: dir_misfit

    ! Print progress and other information
    logical :: verbose = .false.

    !==========================================================================
    ! Dimension

    ! Time steps
    integer :: nt

    ! Time step size
    real :: dt

    ! Max time
    real :: tmax

    !==========================================================================
    ! Geometry

    ! Overall geometry file
    character(len=1204) :: file_geometry

    ! Number of shots
    integer :: ns = 1

    ! Select or exclude sources/receivers
    integer, allocatable, dimension(:, :) :: shot_in_group
    character(len=1024) :: shot_prefix
    integer :: ishot
    integer, allocatable, dimension(:) :: shot_index, rec_index
    integer, allocatable, dimension(:) :: src_exclude, sid_exclude, rec_exclude
    integer, allocatable, dimension(:) :: sid_select, src_select

    ! Windowing geometry to save time
    real :: sxmin = -float_huge
    real :: sxmax = +float_huge
    real :: symin = -float_huge
    real :: symax = +float_huge
    real :: szmin = -float_huge
    real :: szmax = +float_huge
    real :: rxmin = -float_huge
    real :: rxmax = +float_huge
    real :: rymin = -float_huge
    real :: rymax = +float_huge
    real :: rzmin = -float_huge
    real :: rzmax = +float_huge
    real :: offset_min = 0.0d0
    real :: offset_max = float_huge
    integer :: sid_min = 1
    integer :: sid_max = int4_huge
    integer :: shot_min = 1
    integer :: shot_max = 1
    integer :: shot_every = 1
    integer :: rec_min = 1
    integer :: rec_max = int4_huge
    integer :: rec_every = 1

    ! Rotate geometry
    real :: grotx, groty, grotz

    ! Origin for geometry rotation
    real :: gox, goy, goz

    real :: drotx, droty, drotz

    ! Time factor in generating analytical source time function
    ! f0_factor/f0 is the source origin time
    real :: f0_factor = 1.0

    ! Source time function filtering
    real, allocatable, dimension(:) :: src_filt_freqs, src_filt_coefs

    !==========================================================================
    ! Model

    ! Medium type
    character(len=24) :: which_medium

    ! Whether to use adaptive model range for acceleartion
    logical :: yn_adpx
    logical :: yn_adpy
    logical :: yn_adpz

    real :: adp_extrax = 0.0
    real :: adp_extray = 0.0
    real :: adp_extraz = 0.0

    character(len=24) :: aniso_param = 'thomsen'

    real :: min_vpvsratio, max_vpvsratio

    character(len=24), allocatable, dimension(:) :: model_name, model_name_aux
    integer :: nmodel, nmodel_aux
    real, allocatable, dimension(:) :: model_step, model_step_max
    real, allocatable, dimension(:) :: model_min, model_max

    logical :: require_model_interp = .false.

    !==========================================================================
    ! Free surface

    logical :: yn_free_surface = .false.
    character(len=1024) :: file_topo = ''
    character(len=12) :: topo_interp = 'cubic'
    logical :: measure_source_depth_from_surface = .false.
    logical :: measure_receiver_depth_from_surface = .false.
    logical :: source_vertical_to_surface = .false.
    logical :: receiver_vertical_to_surface = .false.
    real :: free_surface_dz_refine = 4.0
    logical :: yn_save_mesh = .false.

    !==========================================================================
    ! Data

    ! Resampling
    real :: data_dt
    real :: data_tmax

    ! Component
    logical :: yn_compx = .true.
    logical :: yn_compy = .true.
    logical :: yn_compz = .true.
    logical :: yn_compp = .false.

    character(len=24), allocatable, dimension(:) :: data_name
    integer :: ndata

    ! Snapshot times
    real, allocatable, dimension(:) :: snaps

    !==========================================================================
    ! Data processing

    character(len=32), allocatable, dimension(:) :: record_processing
    character(len=32), allocatable, dimension(:) :: adjoint_source_processing
    character(len=32), allocatable, dimension(:) :: encoded_record_processing
    character(len=32), allocatable, dimension(:) :: synthetic_processing
    logical :: synthetic_processed

    ! Mute
    character(len=1024) :: dir_topmute_time, dir_btmmute_time
    real :: topmute_vel, btmmute_vel
    real :: topmute_width, topmute_taper
    real :: btmmute_width, btmmute_taper
    real :: btmmute_fromtop = 0.0
    real :: topmute_shift
    real :: topmute_const, btmmute_const
    character(len=1024) :: dir_surgicalmute_time
    real, allocatable, dimension(:) :: surgicalmute_vel, surgicalmute_taper
    real :: surgicalmute_width
    logical :: surgicalmute_inverse = .false.

    ! T-power
    real :: tpow

    ! Time shift
    real :: tshift = 0.0

    ! Frequency filtering
    real, allocatable, dimension(:) :: freqs, amps

    ! Dip filtering
    real, allocatable, dimension(:) :: dips, dipamps

    ! Anisotropic-diffusion-based filtering
    type(andf_param) :: proc_andfparam

    ! Moving-window RMS (i.e., AGC)
    real :: movingbalwin

    !==========================================================================
    ! Inversion related

    ! Six components in a moment tensor
    integer :: nc_mt = 6

    ! Do wavefield cross-correlation every cc_step_interval steps
    integer :: cc_step_interval = 1

    ! Only compute misfit
    logical :: yn_misfit_only

    ! Search direction computation method in inversion
    character(len=32), allocatable, dimension(:) :: model_search_method
    character(len=32) :: search_method

    ! Tikhonov regularization
    integer :: tikhonov_order

    ! Regularization for inversion
    character(len=32), allocatable, dimension(:) :: model_regularization_method
    character(len=32), allocatable, dimension(:) :: source_regularization_method
    logical :: yn_regularize_model = .false.
    logical :: yn_regularize_source = .false.
    real :: reg_lambda

    integer :: iter
    integer :: niter_max
    character(len=32) :: misfit_type

    ! Energy precondition
    logical :: yn_energy_precond

    ! vp/vs ratio smoothing
    real :: vpvsratio_smoothx = 0
    real :: vpvsratio_smoothy = 0
    real :: vpvsratio_smoothz = 0

    character(len=32) :: step_size_method
    real :: step_scaling_factor
    real :: step_max_scale_factor = 1.0

    real :: val_misfit

    logical :: yn_save_adjsrc

    logical :: yn_reconstruct

    real :: tlag_max
    character(len=32) :: penalty_method = 'linear'
    real :: penalty_power = 1.0
    real :: deconv_eps = 0.1
    real :: lawi_sigma = 0.25
    logical :: yn_average_misfit = .false.
    integer :: adj_nt

    integer :: dtw_niter = 5
    real :: dtw_rinst = 0.1
    real :: dtw_rcuml = 1
    real :: dtw_epsabs = 1.0e-2
    real :: dtw_epsrel = 1.0e-3
    character(len=12) :: dtw_loss = 'l2'
    real :: dtw_smooth_median = 0.0
    real :: dtw_smooth_gaussian = 0.0
    character(len=12) :: dtw_form = 'phase'
    integer :: dtw_trc = 0
    real :: dtw_amp_weight = 0.01

    real, allocatable, dimension(:) :: step_misfit
    real, allocatable, dimension(:) :: data_misfit
    real, allocatable, dimension(:, :) :: shot_misfit

    real :: jumpout_factor = 1.0

    character(len=1024) :: file_data_misfit, file_shot_misfit

    logical :: yn_enforce_update = .false.
    logical :: trigger_jumpout = .false.
    logical :: put_synthetic_in_scratch = .false.

    logical :: yn_shared_model_processing = .true.
    character(len=32), allocatable, dimension(:) :: gradient_processing
    character(len=32), allocatable, dimension(:) :: search_direction_processing
    character(len=32), allocatable, dimension(:) :: process_shot_grad, process_grad

    character(len=128) :: kernel_v, kernel_a, prev_kernel_v, prev_kernel_a
    logical :: kernel_type_changed = .false.

    integer :: resume_from_iter
    logical :: yn_flat_stop = .false.
    logical :: yn_continue_inv = .false.

    real :: envelope_p = 2.0

    ! Some data may have extremely unbalanced trace amplitudes
    ! In this case, it may be good to discard small-amplitude traces
    real :: trace_discard_threshold = 1.0e-6

    ! In the Adam optimizer, β_1 and β_2 are exponential decay rates for computing
    ! the moving averages of the gradient and the squared gradient, respectively.
    real :: adam_beta1 = 0.9
    real :: adam_beta2 = 0.999
    real :: adam_eps = 1.0e-8

    integer :: htlen = 30

    logical :: yn_update_medium = .true.
    logical :: yn_update_source = .false.

contains

    !
    !> Read modeling or FWI parameters from parameter file
    !
    subroutine read_parameters

        character(len=1024) :: temporary_parameter_file

#ifdef _fwi_
        integer :: i
        real :: valmin, valmax
#endif

        if (command_argument_count() == 0) then
            if (rankid == 0) then
                call warn('')
                call warn(date_time_compact()//' Error: The program needs a parameter file. Exiting. ')
                call warn('')
            end if
            call mpibarrier
            call mpistop
        end if

        ! get name of the parameter file
        call get_command_argument(1, file_parameter)

        ! if parameter file does not exist, exit
        if (.not. file_exists(file_parameter)) then
            call warn(date_time_compact()//' Error: Parameter file = '//tidy(file_parameter)//' not found. Exiting. ')
            call mpibarrier
            call mpistop
        end if

        ! copy parameter file to working directory for log purpose
#ifdef _forward_
        call readpar_string(file_parameter, 'dir_synthetic', dir_synthetic, './data_synthetic')
        call make_directory(dir_synthetic)
        if (rankid == 0) then
            temporary_parameter_file = tidy(dir_synthetic)//'/parameters.forward.'//date_time_string()
            call copy_file(file_parameter, temporary_parameter_file)
            file_parameter = temporary_parameter_file
            call warn(' Parameter file: '//tidy(file_parameter))
        end if
#endif

#ifdef _fwi_
        call readpar_string(file_parameter, 'dir_working', dir_working, './test')
        call make_directory(dir_working)
        if (rankid == 0) then
            temporary_parameter_file = tidy(dir_working)//'/parameters.fwi.'//date_time_string()
            call copy_file(file_parameter, temporary_parameter_file)
            file_parameter = temporary_parameter_file
            call warn(' Parameter file: '//tidy(file_parameter))
        end if
#endif

        call mpibarrier

        call bcast(file_parameter)

        ! 2D does not support domain decomposition yet, for the purpose of simplicity
#ifdef _dim2_
        ngroup = nrank
        rank1_group = 1
        rank2_group = 1
        rank3_group = 1
#endif

        ! 3D supports domain decomposition, to make modeling/FWI in large models faster
#ifdef _dim3_
        call readpar_int(file_parameter, 'rankx', rank3_group, 1)
        call readpar_int(file_parameter, 'ranky', rank2_group, 1)
        call readpar_int(file_parameter, 'rankz', rank1_group, 1)
        call readpar_int(file_parameter, 'ngroup', ngroup, max(floor(nrank*1.0/(rank1_group*rank2_group*rank3_group)), 1))
#endif

        call readpar_float(file_parameter, 'dt', dt, 0.0)
        call readpar_float(file_parameter, 'tmax', tmax, 1.0)
        nt = nint(tmax/dt + 1)

        call readpar_float(file_parameter, 'data_dt', data_dt, dt)
        call readpar_float(file_parameter, 'data_tmax', data_tmax, tmax)

        call readpar_int(file_parameter, 'cc_step_interval', cc_step_interval, 1)

        call readpar_string(file_parameter, 'which_medium', which_medium, 'acoustic-iso')
        call readpar_string(file_parameter, 'anisotropy_type', aniso_param, 'iso')

        call readpar_float(file_parameter, 'min_vpvsratio', min_vpvsratio, 1.1)
        call readpar_float(file_parameter, 'max_vpvsratio', max_vpvsratio, 9.0)
        call readpar_float(file_parameter, 'vpvsratio_smoothx', vpvsratio_smoothx, 0.0)
        call readpar_float(file_parameter, 'vpvsratio_smoothy', vpvsratio_smoothy, 0.0)
        call readpar_float(file_parameter, 'vpvsratio_smoothz', vpvsratio_smoothz, 0.0)

        !================================================================================

        call readpar_int(file_parameter, 'ns', ns, 1)
        call readpar_string(file_parameter, 'file_geometry', file_geometry, '', required=.true.)
        call readpar_nint(file_parameter, 'src_index', shot_index, [1, 1, ns])
        shot_min = max(shot_index(1), 1)
        shot_every = max(shot_index(2), 1)
        shot_max = min(shot_index(3), ns)
        call readpar_nint(file_parameter, 'rec_index', rec_index, [1, 1, int4_huge])
        rec_min = max(rec_index(1), 1)
        rec_every = max(rec_index(2), 1)
        rec_max = min(rec_index(3), int4_huge)
        call readpar_nint(file_parameter, 'rec_exclude', rec_exclude, [0])
        call readpar_float(file_parameter, 'offset_min', offset_min, 0.0)
        call readpar_float(file_parameter, 'offset_max', offset_max, float_huge)
        call readpar_nint(file_parameter, 'sid_select', sid_select, [0])
        call readpar_nint(file_parameter, 'src_select', src_select, [0])
        call readpar_nint(file_parameter, 'sid_exclude', sid_exclude, [0])
        call readpar_nint(file_parameter, 'src_exclude', src_exclude, [0])

        call readpar_float(file_parameter, 'geometry_rotate_x', grotx, 0.0)
        call readpar_float(file_parameter, 'geometry_rotate_y', groty, 0.0)
        call readpar_float(file_parameter, 'geometry_rotate_z', grotz, 0.0)
        call readpar_float(file_parameter, 'geometry_ox', gox, 0.0)
        call readpar_float(file_parameter, 'geometry_oy', goy, 0.0)
        call readpar_float(file_parameter, 'geometry_oz', goz, 0.0)

        call readpar_int(file_parameter, 'sid_min', sid_min, 0)
        call readpar_int(file_parameter, 'sid_max', sid_max, int4_huge)
        call readpar_float(file_parameter, 'sx_min', sxmin, -float_huge)
        call readpar_float(file_parameter, 'sx_max', sxmax, +float_huge)
        call readpar_float(file_parameter, 'sy_min', symin, -float_huge)
        call readpar_float(file_parameter, 'sy_max', symax, +float_huge)
        call readpar_float(file_parameter, 'sz_min', szmin, -float_huge)
        call readpar_float(file_parameter, 'sz_max', szmax, +float_huge)
        call readpar_float(file_parameter, 'rx_min', rxmin, -float_huge)
        call readpar_float(file_parameter, 'rx_max', rxmax, +float_huge)
        call readpar_float(file_parameter, 'ry_min', rymin, -float_huge)
        call readpar_float(file_parameter, 'ry_max', rymax, +float_huge)
        call readpar_float(file_parameter, 'rz_min', rzmin, -float_huge)
        call readpar_float(file_parameter, 'rz_max', rzmax, +float_huge)

        call readpar_float(file_parameter, 'f0_factor', f0_factor, 1.0)
        call readpar_nfloat(file_parameter, 'src_filt_freqs', src_filt_freqs, [-1.0])
        call readpar_nfloat(file_parameter, 'src_filt_coefs', src_filt_coefs, [-1.0])
        call assert(size(src_filt_freqs) == size(src_filt_coefs), 'Error: size(src_filt_freqs) must = size(src_filt_coefs)')

        !================================================================================
        call readpar_string(file_parameter, 'dir_working', dir_working, './test')
        call readpar_string(file_parameter, 'dir_record', dir_record, './data')
        call readpar_string(file_parameter, 'dir_synthetic', dir_synthetic, './data_synthetic')
        call readpar_string(file_parameter, 'dir_base', dir_base, './data_base')
        call readpar_string(file_parameter, 'dir_processed', dir_synthetic_processed, './data_processed')
        call readpar_string(file_parameter, 'dir_snapshot', dir_snapshot, './snapshot')

        call readpar_nfloat(file_parameter, 'snaps', snaps, [-1.0])

        call readpar_logical(file_parameter, 'yn_adpx', yn_adpx, .false.)
        call readpar_logical(file_parameter, 'yn_adpy', yn_adpy, .false.)
        call readpar_logical(file_parameter, 'yn_adpz', yn_adpz, .false.)
        call readpar_float(file_parameter, 'adp_extrax', adp_extrax, 0.0)
        call readpar_float(file_parameter, 'adp_extray', adp_extray, 0.0)
        call readpar_float(file_parameter, 'adp_extraz', adp_extraz, 0.0)

        call readpar_logical(file_parameter, 'verbose', verbose, .false.)

        call readpar_logical(file_parameter, 'yn_free_surface', yn_free_surface, .false.)
        call readpar_float(file_parameter, 'free_surface_dz_refine', free_surface_dz_refine, 4.0)
        call readpar_string(file_parameter, 'file_topo', file_topo, '')
        call readpar_string(file_parameter, 'topo_interp', topo_interp, 'cubic')
        call readpar_logical(file_parameter, 'measure_source_depth_from_surface', measure_source_depth_from_surface, .false.)
        call readpar_logical(file_parameter, 'measure_receiver_depth_from_surface', measure_receiver_depth_from_surface, .false.)
        call readpar_logical(file_parameter, 'source_vertical_to_surface', source_vertical_to_surface, .false.)
        call readpar_logical(file_parameter, 'receiver_vertical_to_surface', receiver_vertical_to_surface, .false.)
        call readpar_logical(file_parameter, 'yn_save_mesh', yn_save_mesh, .false.)

        call readpar_string(file_parameter, 'dp_dir_top_mute_time', dir_topmute_time, '')
        call readpar_float(file_parameter, 'dp_top_mute_vel', topmute_vel, -1.0)
        call readpar_float(file_parameter, 'dp_top_mute_width', topmute_width, 0.1*tmax)
        call readpar_float(file_parameter, 'dp_top_mute_taper', topmute_taper, 0.1*tmax)
        call readpar_float(file_parameter, 'dp_top_mute_shift', topmute_shift, 0.0)
        call readpar_float(file_parameter, 'dp_top_mute_const', topmute_const, -tmax)

        call readpar_string(file_parameter, 'dp_dir_bottom_mute_time', dir_btmmute_time, '')
        call readpar_float(file_parameter, 'dp_bottom_mute_vel', btmmute_vel, -2.0)
        call assert(topmute_vel > btmmute_vel, ' <read_parameter> Error: dp_top_mute_vel must > dp_bottom_mute_vel.')
        call readpar_float(file_parameter, 'dp_bottom_mute_width', btmmute_width, 0.1*tmax)
        call readpar_float(file_parameter, 'dp_bottom_mute_taper', btmmute_taper, 0.1*tmax)
        call readpar_float(file_parameter, 'dp_bottom_mute_from_top', btmmute_fromtop, 0.0)
        call readpar_float(file_parameter, 'dp_bottom_mute_const', btmmute_const, -tmax)

        call readpar_string(file_parameter, 'dp_dir_surgical_mute_time', dir_surgicalmute_time, '')
        call readpar_nfloat(file_parameter, 'dp_surgical_mute_vel', surgicalmute_vel, [-1.0, -1.0])
        call readpar_float(file_parameter, 'dp_surgical_mute_width', surgicalmute_width, 0.1*tmax)
        call readpar_nfloat(file_parameter, 'dp_surgical_mute_taper', surgicalmute_taper, [0.1*tmax, 0.1*tmax])
        call assert(size(surgicalmute_vel) == 2, ' <read_parameter> Error: size(surgical_mute_vel) must = 2')
        call assert(size(surgicalmute_taper) == 2, ' <read_parameter> Error: size(surgical_mute_taper) must = 2')
        call readpar_logical(file_parameter, 'dp_surgical_mute_inverse', surgicalmute_inverse, .false.)

        call readpar_nfloat(file_parameter, 'dp_dip_filt_dips', dips, [-1000.0, 0.0, 1000.0])
        call readpar_nfloat(file_parameter, 'dp_dip_filt_coefs', dipamps, [0.0, 1.0, 0.0])
        call readpar_float(file_parameter, 'dp_andf_smoothx', proc_andfparam%smooth2, 2.0)
        call readpar_float(file_parameter, 'dp_andf_smootht', proc_andfparam%smooth1, tmax/20.0)
        call readpar_float(file_parameter, 'dp_andf_powerm', proc_andfparam%powerm, 1.0)
        call readpar_int(file_parameter, 'dp_andf_t', proc_andfparam%niter, 10)
        call readpar_float(file_parameter, 'dp_andf_sigma', proc_andfparam%sigma, 10.0)
        call readpar_float(file_parameter, 'dp_andf_alpha', proc_andfparam%lambda1, 0.0)
        call readpar_float(file_parameter, 'dp_andf_beta', proc_andfparam%lambda2, 1.0)
        call readpar_float(file_parameter, 'dp_t_balance_window', movingbalwin, tmax/20.0)
        call readpar_float(file_parameter, 'dp_t_shift', tshift, 0.0)
        call readpar_float(file_parameter, 'dp_t_power', tpow, 0.0)
        call readpar_nfloat(file_parameter, 'dp_freq_filt_freqs', freqs, [-1.0])
        call readpar_nfloat(file_parameter, 'dp_freq_filt_coefs', amps, [-1.0])
        call readpar_float(file_parameter, 'dp_rotate_x', drotx, 0.0)
        call readpar_float(file_parameter, 'dp_rotate_y', droty, 0.0)
        call readpar_float(file_parameter, 'dp_rotate_z', drotz, 0.0)

#ifdef _forward_

        call readpar_nstring(file_parameter, 'model_name', model_name, [''], required=.true.)
        call assert(model_name(1) /= '', ' <read_parameter> Error: model_name cannot be empty or start with null')
        nmodel = size(model_name)

        call readpar_nstring(file_parameter, 'process_synthetic', synthetic_processing, [''])

#endif

#ifdef _fwi_

        !================================================================================
        call readpar_nstring(file_parameter, 'model_regularization_method', model_regularization_method, [''])
        call readpar_nstring(file_parameter, 'source_regularization_method', source_regularization_method, [''])
        if (model_regularization_method(1) /= '') then
            yn_regularize_model = .true.
        end if
        if (source_regularization_method(1) /= '') then
            yn_regularize_source = .true.
        end if

        call readpar_int(file_parameter, 'niter_max', niter_max, 100)

        !=======================================================================
        ! Misfit parameters

        call readpar_string(file_parameter, 'misfit_type', misfit_type, 'waveform')
        call readpar_float(file_parameter, 'tlag_max', tlag_max, 0.1*tmax)
        call readpar_int(file_parameter, 'adj_nt', adj_nt, nt)
        call readpar_float(file_parameter, 'envelope_power', envelope_p, 2.0)
        call readpar_string(file_parameter, 'penalty_method', penalty_method, 'linear')
        call readpar_float(file_parameter, 'deconv_eps', deconv_eps, 0.1)
        call readpar_float(file_parameter, 'lawi_sigma', lawi_sigma, 0.25)
        call readpar_float(file_parameter, 'penalty_power', penalty_power, 1.0)
        call readpar_logical(file_parameter, 'yn_average_misfit', yn_average_misfit, .false.)
        call readpar_float(file_parameter, 'dtw_smooth_median', dtw_smooth_median, 0.0)
        call readpar_float(file_parameter, 'dtw_smooth_gaussian', dtw_smooth_gaussian, 0.0)
        call readpar_int(file_parameter, 'dtw_niter', dtw_niter, 5)
        call readpar_float(file_parameter, 'dtw_rinst', dtw_rinst, 0.1)
        call readpar_float(file_parameter, 'dtw_rcuml', dtw_rcuml, 1.0)
        call readpar_float(file_parameter, 'dtw_epsabs', dtw_epsabs, 1.0e-2)
        call readpar_float(file_parameter, 'dtw_epsrel', dtw_epsrel, 1.0e-3)
        call readpar_string(file_parameter, 'dtw_loss', dtw_loss, 'l0.5')
        call readpar_string(file_parameter, 'dtw_form', dtw_form, 'phase')
        call readpar_float(file_parameter, 'dtw_amp_weight', dtw_amp_weight, 1.0)
        call readpar_int(file_parameter, 'dtw_trc', dtw_trc, 0)
        call readpar_float(file_parameter, 'trace_discard_threshold', trace_discard_threshold, 1.0e-6)

        !=======================================================================
        call readpar_string(file_parameter, 'search_method', search_method, 'cg')

        call readpar_float(file_parameter, 'adam_beta1', adam_beta1, 0.9)
        call readpar_float(file_parameter, 'adam_beta2', adam_beta2, 0.999)
        call readpar_float(file_parameter, 'adam_eps', adam_eps, 1.0e-8)

        call readpar_string(file_parameter, 'step_size_method', step_size_method, 'linear')

        call readpar_string(file_parameter, 'file_data_misfit', file_data_misfit, tidy(dir_working)//'/data_misfit.txt')
        call readpar_string(file_parameter, 'file_shot_misfit', file_shot_misfit, tidy(dir_working)//'/shot_misfit.bin')

        call readpar_logical(file_parameter, 'yn_continue', yn_continue_inv, .false.)
        if (yn_continue_inv) then
            resume_from_iter = max(1, count_nonempty_lines(file_data_misfit) - 1)
        else
            call readpar_int(file_parameter, 'resume_from_iter', resume_from_iter, 1)
        end if

        call readpar_nstring(file_parameter, 'process_record', record_processing, [''])
        call readpar_nstring(file_parameter, 'process_record_encoded', encoded_record_processing, [''])
        call readpar_nstring(file_parameter, 'process_synthetic', synthetic_processing, [''])
        call readpar_nstring(file_parameter, 'process_adjsrc', adjoint_source_processing, [''])
        call readpar_nstring(file_parameter, 'process_shot_grad', process_shot_grad, [''])
        call readpar_nstring(file_parameter, 'process_grad', gradient_processing, [''])
        call readpar_nstring(file_parameter, 'process_srch', search_direction_processing, [''])

        call readpar_logical(file_parameter, 'yn_shared_model_processing', yn_shared_model_processing, .true.)

        call readpar_logical(file_parameter, 'yn_energy_precond', yn_energy_precond, .false.)

        call readpar_logical(file_parameter, 'yn_flat_stop', yn_flat_stop, .false.)

        !=======================================================================
        ! Model update parameters

        call readpar_nstring(file_parameter, 'model_update', model_name, ['vp'])
        call assert(model_name(1) /= '', ' <read_parameter> Error: model_update cannot be empty or start with null. ')
        nmodel = size(model_name)

        call readpar_nstring(file_parameter, 'model_aux', model_name_aux, [''])
        if (model_name_aux(1) == '') then
            nmodel_aux = 0
        else
            nmodel_aux = size(model_name_aux)
        end if

        model_min = zeros(nmodel)
        model_max = zeros(nmodel)
        model_step_max = zeros(nmodel)
        model_step = zeros(nmodel)
        allocate(model_search_method(1:nmodel))
        do i = 1, nmodel

            select case (remove_string_after(model_name(i), ['c', 'C']))
                case ('vp', 'vs', 'rho')
                    valmin = 0.0
                    valmax = 1.0e5
                case ('epsilon', 'delta', 'gamma', 'eta')
                    valmin = 0.0
                    valmax = 0.5
                case ('theta', 'phi')
                    valmin = 0.0
                    valmax = const_pi
                case ('c', 'C')
                    valmin = 0.0
                    valmax = 1.0e9
                case ('qp', 'qs')
                    valmin = 10.0
                    valmax = float_large
                case ('mt')
                    valmin = -1.0e9
                    valmax = 1.0e9
            end select
            call readpar_float(file_parameter, 'min_'//tidy(model_name(i)), model_min(i), valmin)
            call readpar_float(file_parameter, 'max_'//tidy(model_name(i)), model_max(i), valmax)

            ! The program allows using different search method for different parameters
            call readpar_string(file_parameter, 'search_method_'//tidy(model_name(i)), model_search_method(i), search_method)

        end do

        if (any(model_name /= 'mt' .and. model_name /= 'stf')) then
            yn_update_medium = .true.
        else
            yn_update_medium = .false.
        end if

        if (any(model_name == 'mt' .or. model_name == 'stf')) then
            yn_update_source = .true.
        else
            yn_update_source = .false.
        end if

#endif

        select case (which_medium)

            case ('acoustic-iso', 'acoustic-tti')
                call readpar_nstring(file_parameter, 'data_name', data_name, ['p'])

            case ('elastic-iso', 'elastic-vhtiort', 'elastic-tti')
#ifdef _dim2_
                call readpar_nstring(file_parameter, 'data_name', data_name, ['x', 'z'])
#endif

#ifdef _dim3_
                call readpar_nstring(file_parameter, 'data_name', data_name, ['x', 'y', 'z'])
#endif

                if (any(data_name == 'x')) then
                    yn_compx = .true.
                else
                    yn_compx = .false.
                end if

                if (any(data_name == 'y')) then
                    yn_compy = .true.
                else
                    yn_compy = .false.
                end if

                if (any(data_name == 'z')) then
                    yn_compz = .true.
                else
                    yn_compz = .false.
                end if

        end select

        ndata = size(data_name)

        call assert(ndata >= 1, ' <read_parameter> Error: data_name is empty. ')
        call assert(all_in(data_name, ['x', 'y', 'z', 'p']), &
            ' <read_parameter> Error: data_name must be in x, y, z, p. ')

    end subroutine read_parameters

end module
