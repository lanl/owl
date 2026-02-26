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


module gradient

    use mod_parameters
    use mod_model
    use mod_source_receiver
    use mod_data_processing
    use vars
    use inversion_adjoint_source
    use acoustic_iso_2d
    use elastic_vhtiort_2d
    use elastic_tti_2d

    implicit none

contains

    !
    !> Initialize gradient arrays
    !
    subroutine zero_gradient

        integer :: i

        do i = 1, nmodel
            select case (model_name(i))
                case default
                    model_grad(i)%array = zeros(nz, nx)
                case ('mt')
                    model_grad(i)%array = zeros(nc_mt, ns)
            end select
        end do

    end subroutine zero_gradient

    !
    !> Compute gradients shot by shot and sum
    !
    subroutine compute_gradient_shots

        character(len=1024) :: dir_field, dir_from, dir_to
        integer :: i
        type(grid2) :: grd
        type(wave_solver_acoustic_iso_2d) :: solver_acoustic_iso
        type(wave_solver_elastic_vhtiort_2d) :: solver_elastic_vhtiort
        type(wave_solver_elastic_tti_2d) :: solver_elastic_tti

        dir_scratch = tidy(dir_working)//'/scratch'
        dir_synthetic = dir_iter_synthetic(iter)
        dir_synthetic_processed = dir_iter_synthetic_processed(iter)
        dir_field = dir_iter_record(iter)
        dir_adjoint = dir_iter_adjoint_source(iter)
        if (rankid == 0) then
            call make_directory(dir_synthetic)
            call make_directory(dir_scratch)
        end if
        call mpibarrier

        ! Set current step misfit to zero
        step_misfit = 0.0d0

        ! Enforce Vp/Vs ratio in an appropriate range
        call clip_vpvsratio

        ! Get kernel type of medium parameters; can change in different iterations
        call readpar_xstring(file_parameter, 'kernel_v', kernel_v, 'full', iter*1.0)
        call readpar_xstring(file_parameter, 'kernel_a', kernel_a, 'full', iter*1.0)

        if (iter > 1) then
            call readpar_xstring(file_parameter, 'kernel_v', prev_kernel_v, 'full', iter - 1.0)
            call readpar_xstring(file_parameter, 'kernel_a', prev_kernel_a, 'full', iter - 1.0)
            if (prev_kernel_v /= kernel_v .or. prev_kernel_a /= kernel_v) then
                kernel_type_changed = .true.
            else
                kernel_type_changed = .false.
            end if
        end if

        ! Show the medium parameter statistics
        if (rankid == 0) then
            do i = 1, nmodel
                call plot_histogram(model_m(i)%array, &
                    label=date_time_compact()//' '//tidy(model_name(i))//' distribution ')
            end do
        end if
        call mpibarrier

        ! FWI gradient computation, shot paralellization with hybrid OpenMP + MPI for each shot
        do ishot = shot_in_group(groupid, 1), shot_in_group(groupid, 2)

            ! Shot prefix
            shot_prefix = tidy(dir_scratch)//'/shot_'//num2str(set_srcid(ishot))

            ! Setup geometry
            call set_adaptive_range(gmtr(ishot))

            !            if (yn_source_encoding) then
            !                call sgmtr%merge( &
                !                    gmtr(shot_in_super(ishot, 1):shot_in_super(ishot, 2)), &
                !                    encoding_id=ishot, &
                !                    encoding_beg=shot_in_super(ishot, 1), &
                !                    encoding_amp=amplitude_matrix, &
                !                    encoding_phase=delay_matrix)
            !            else
            !                sgmtr = gmtr(ishot)
            !            end if

            select case (which_medium)

                case ('acoustic-iso')
                    solver_acoustic_iso%nx = shot_nx
                    solver_acoustic_iso%nz = shot_nz
                    solver_acoustic_iso%dx = dx
                    solver_acoustic_iso%dz = dz
                    solver_acoustic_iso%ox = shot_xbeg
                    solver_acoustic_iso%oz = shot_zbeg
                    solver_acoustic_iso%dt = dt
                    solver_acoustic_iso%tmax = tmax
                    solver_acoustic_iso%data_dt = data_dt
                    solver_acoustic_iso%data_tmax = data_tmax
                    solver_acoustic_iso%dir_synthetic = tidy(dir_synthetic)
                    solver_acoustic_iso%dir_working = tidy(dir_scratch)
                    solver_acoustic_iso%pml = pml
                    solver_acoustic_iso%free_surface = yn_free_surface
                    solver_acoustic_iso%verbose = verbose
                    solver_acoustic_iso%gmtr = gmtr(ishot)
                    solver_acoustic_iso%reconstruct = yn_reconstruct
                    solver_acoustic_iso%vp = get_model('vp')
                    solver_acoustic_iso%rho = get_model('rho', 1.0)

                    call solver_acoustic_iso%forward

                case ('elastic-iso', 'elastic-vhtiort')
                    solver_elastic_vhtiort%nx = shot_nx
                    solver_elastic_vhtiort%nz = shot_nz
                    solver_elastic_vhtiort%dx = dx
                    solver_elastic_vhtiort%dz = dz
                    solver_elastic_vhtiort%ox = shot_xbeg
                    solver_elastic_vhtiort%oz = shot_zbeg
                    solver_elastic_vhtiort%dt = dt
                    solver_elastic_vhtiort%tmax = tmax
                    solver_elastic_vhtiort%data_dt = data_dt
                    solver_elastic_vhtiort%data_tmax = data_tmax
                    solver_elastic_vhtiort%dir_synthetic = tidy(dir_synthetic)
                    solver_elastic_vhtiort%dir_working = tidy(dir_scratch)
                    if (sum(snaps) > 0) then
                        solver_elastic_vhtiort%dir_snapshot = tidy(dir_snapshot)
                        solver_elastic_vhtiort%snaps = regspace(snaps(1), snaps(2), snaps(3))
                    end if
                    solver_elastic_vhtiort%pml = pml
                    solver_elastic_vhtiort%free_surface = yn_free_surface
                    solver_elastic_vhtiort%free_surface_dz_refine = free_surface_dz_refine
                    solver_elastic_vhtiort%dz_max = dz_max
                    solver_elastic_vhtiort%verbose = verbose
                    solver_elastic_vhtiort%gmtr = gmtr(ishot)
                    solver_elastic_vhtiort%reconstruct = yn_reconstruct
                    solver_elastic_vhtiort%compx = yn_compx
                    solver_elastic_vhtiort%compz = yn_compz
                    solver_elastic_vhtiort%yn_update_medium = yn_update_medium
                    solver_elastic_vhtiort%yn_update_source = yn_update_source
                    solver_elastic_vhtiort%mt = slice(get_model('mt', 0.0), dim=2, index=ishot)

                    select case (aniso_param)

                        case ('iso')
                            solver_elastic_vhtiort%anisotropy_type = 'iso'
                            solver_elastic_vhtiort%vp = get_model('vp')
                            solver_elastic_vhtiort%vs = get_model('vs')
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                        case ('thomsen')
                            solver_elastic_vhtiort%anisotropy_type = 'thomsen'
                            solver_elastic_vhtiort%vp = get_model('vp')
                            solver_elastic_vhtiort%vs = get_model('vs')
                            solver_elastic_vhtiort%tieps = get_model('epsilon', 0.0)
                            solver_elastic_vhtiort%tidel = get_model('delta', 0.0)
                            solver_elastic_vhtiort%tithe = get_model('theta', 0.0)
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                        case ('a-t')
                            solver_elastic_vhtiort%anisotropy_type = 'a-t'
                            solver_elastic_vhtiort%vp = get_model('vp')
                            solver_elastic_vhtiort%vs = get_model('vs')
                            solver_elastic_vhtiort%tieps = get_model('epsilon', 0.0)
                            solver_elastic_vhtiort%tieta = get_model('eta', 0.0)
                            solver_elastic_vhtiort%tithe = get_model('theta', 0.0)
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                        case ('cij')
                            solver_elastic_vhtiort%anisotropy_type = 'cij'
                            solver_elastic_vhtiort%c11 = get_model('c11')
                            solver_elastic_vhtiort%c13 = get_model('c13')
                            solver_elastic_vhtiort%c33 = get_model('c33')
                            solver_elastic_vhtiort%c55 = get_model('c55')
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                    end select

                    call solver_elastic_vhtiort%forward

                case ('elastic-tti')
                    solver_elastic_tti%nx = shot_nx
                    solver_elastic_tti%nz = shot_nz
                    solver_elastic_tti%dx = dx
                    solver_elastic_tti%dz = dz
                    solver_elastic_tti%ox = shot_xbeg
                    solver_elastic_tti%oz = shot_zbeg
                    solver_elastic_tti%dt = dt
                    solver_elastic_tti%tmax = tmax
                    solver_elastic_tti%data_dt = data_dt
                    solver_elastic_tti%data_tmax = data_tmax
                    solver_elastic_tti%dir_synthetic = tidy(dir_synthetic)
                    solver_elastic_tti%dir_working = tidy(dir_scratch)
                    if (sum(snaps) > 0) then
                        solver_elastic_tti%dir_snapshot = tidy(dir_snapshot)
                        solver_elastic_tti%snaps = regspace(snaps(1), snaps(2), snaps(3))
                    end if
                    solver_elastic_tti%pml = pml
                    solver_elastic_tti%free_surface = yn_free_surface
                    solver_elastic_tti%free_surface_dz_refine = free_surface_dz_refine
                    solver_elastic_tti%dz_max = dz_max
                    solver_elastic_tti%file_topo = tidy(file_topo)
                    solver_elastic_tti%topo_interp = tidy(topo_interp)
                    solver_elastic_tti%measure_source_depth_from_surface = measure_source_depth_from_surface
                    solver_elastic_tti%measure_receiver_depth_from_surface = measure_receiver_depth_from_surface
                    solver_elastic_tti%source_vertical_to_surface = source_vertical_to_surface
                    solver_elastic_tti%receiver_vertical_to_surface = receiver_vertical_to_surface
                    solver_elastic_tti%verbose = verbose
                    solver_elastic_tti%gmtr = gmtr(ishot)
                    solver_elastic_tti%reconstruct = yn_reconstruct
                    solver_elastic_tti%compx = yn_compx
                    solver_elastic_tti%compz = yn_compz
                    solver_elastic_tti%yn_update_medium = yn_update_medium
                    solver_elastic_tti%yn_update_source = yn_update_source
                    solver_elastic_tti%mt = slice(get_model('mt', 0.0), dim=2, index=ishot)

                    select case (aniso_param)

                        case ('iso')
                            solver_elastic_tti%anisotropy_type = 'iso'
                            solver_elastic_tti%vp = get_model('vp')
                            solver_elastic_tti%vs = get_model('vs')
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                        case ('thomsen')
                            solver_elastic_tti%anisotropy_type = 'thomsen'
                            solver_elastic_tti%vp = get_model('vp')
                            solver_elastic_tti%vs = get_model('vs')
                            solver_elastic_tti%tieps = get_model('epsilon', 0.0)
                            solver_elastic_tti%tidel = get_model('delta', 0.0)
                            solver_elastic_tti%tithe = get_model('theta', 0.0)
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                        case ('a-t')
                            solver_elastic_tti%anisotropy_type = 'a-t'
                            solver_elastic_tti%vp = get_model('vp')
                            solver_elastic_tti%vs = get_model('vs')
                            solver_elastic_tti%tieps = get_model('epsilon', 0.0)
                            solver_elastic_tti%tieta = get_model('eta', 0.0)
                            solver_elastic_tti%tithe = get_model('theta', 0.0)
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                        case ('cij')
                            solver_elastic_tti%anisotropy_type = 'cij'
                            solver_elastic_tti%c11 = get_model('c11')
                            solver_elastic_tti%c13 = get_model('c13')
                            solver_elastic_tti%c15 = get_model('c15')
                            solver_elastic_tti%c33 = get_model('c33')
                            solver_elastic_tti%c35 = get_model('c35')
                            solver_elastic_tti%c55 = get_model('c55')
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                    end select

                    call solver_elastic_tti%forward

            end select

            call warn(date_time_compact()//' Shot '//num2str(set_srcid(ishot))//' forward modeling completed. ')

            ! Process synthetic data
            call process_synthetic(ishot)

            if (sum(data_misfit) == 0) then
                ! At the first iteration, copy synthetic data to iteration 0 directory

                put_synthetic_in_scratch = .false.

                if (synthetic_processed) then
                    call make_directory(dir_iter_synthetic_processed(0))
                    dir_from = tidy(dir_iter_synthetic_processed(1))
                    dir_to = tidy(dir_iter_synthetic_processed(0))
                else
                    call make_directory(dir_iter_synthetic(0))
                    dir_from = tidy(dir_iter_synthetic(1))
                    dir_to = tidy(dir_iter_synthetic(0))
                end if

                do i = 1, ndata
                    call copy_file( &
                        tidy(dir_from)//'/shot_'//num2str(set_srcid(ishot))//'_seismogram_'//tidy(data_name(i))//'.su', &
                        tidy(dir_to)//'/shot_'//num2str(set_srcid(ishot))//'_seismogram_'//tidy(data_name(i))//'.su')
                end do

                ! Copy initial models to dir_working/iteration_0
                if (rankid == 0) then
                    call make_directory(tidy(dir_working)//'/iteration_0/model')
                    do i = 1, nmodel
                        call output_array(model_m(i)%array, tidy(dir_working)//'/iteration_0/model/' &
                            //tidy(model_name(i))//'.bin')
                    end do
                end if

            end if

            ! Compute gradient using adjoint-state method
            if (synthetic_processed) then
                ! if the synthetic is processed
                call compute_adjoint_source(ishot, step_misfit(ishot), dir_synthetic_processed, dir_field, dir_adjoint)
            else
                ! if the synthetic is intact
                call compute_adjoint_source(ishot, step_misfit(ishot), dir_synthetic, dir_field, dir_adjoint)
            end if
            call warn(date_time_compact()//' >> Shot '//num2str(set_srcid(ishot)) &
                //' misfit = '//num2str(step_misfit(ishot), '(es)'))

            if (yn_misfit_only) then
                ! If the step is for computing data misfit then
                ! delete temporary files and go to next shot
                ! call execute_command_line('rm -rf '//tidy(shot_prefix)//'_*')
                cycle
            end if

            call process_adjoint_source(ishot)

            select case (which_medium)

                case ('acoustic-iso')
                    solver_acoustic_iso%dir_adjoint = tidy(dir_adjoint)
                    solver_acoustic_iso%cc_step_interval = cc_step_interval
                    solver_acoustic_iso%energy_precond = yn_energy_precond
                    solver_acoustic_iso%kernel_v = kernel_v
                    solver_acoustic_iso%kernel_a = kernel_a

                    call solver_acoustic_iso%adjoint

                case ('elastic-iso', 'elastic-vhtiort')
                    solver_elastic_vhtiort%dir_adjoint = tidy(dir_adjoint)
                    solver_elastic_vhtiort%cc_step_interval = cc_step_interval
                    solver_elastic_vhtiort%energy_precond = yn_energy_precond
                    solver_elastic_vhtiort%kernel_v = kernel_v
                    solver_elastic_vhtiort%kernel_a = kernel_a

                    call solver_elastic_vhtiort%adjoint

                case ('elastic-tti')
                    solver_elastic_tti%dir_adjoint = tidy(dir_adjoint)
                    solver_elastic_tti%cc_step_interval = cc_step_interval
                    solver_elastic_tti%energy_precond = yn_energy_precond
                    solver_elastic_tti%kernel_v = kernel_v
                    solver_elastic_tti%kernel_a = kernel_a

                    call solver_elastic_tti%adjoint

            end select

            call warn(date_time_compact()//' >> Shot '//num2str(set_srcid(ishot)) &
                //' gradient computation completed. ')

            ! Process computed gradients
            do i = 1, nmodel

                select case (model_name(i))

                    case ('mt', 'stf')

                        call grd%input(tidy(shot_prefix)//'_grad_'//tidy(model_name(i))//'.grd')
                        model_grad(i)%array = model_grad(i)%array + grd%array
                        call warn(date_time_compact()//' Shot '//num2str(set_srcid(ishot))//' '//tidy(model_name(i))//' merged.')

                    case default

                        if(yn_shared_model_processing) then
                            call process_model_single_shot(ishot, model_grad(i)%array, 'grad', &
                                tidy(shot_prefix)//'_grad_'//tidy(model_name(i))//'.grd')
                        else
                            call process_model_single_shot(ishot, model_grad(i)%array, 'grad_'//tidy(model_name(i)), &
                                tidy(shot_prefix)//'_grad_'//tidy(model_name(i))//'.grd')
                        end if

                end select

            end do

            ! ! Remove temporary files
            ! call execute_command_line('rm -rf '//tidy(shot_prefix)//'_*')

        end do

        call mpibarrier

        ! collect misfit
        call allreduce_array(step_misfit)
        if (sum(data_misfit) == 0) then
            data_misfit(0) = sum(step_misfit)
            shot_misfit(:, 0) = step_misfit
        end if
        data_misfit(iter) = sum(step_misfit)

        if (yn_misfit_only) then
            return
        end if

        call mpibarrier

        do i = 1, nmodel
            call allreduce_array(model_grad(i)%array)
        end do

        call mpibarrier

    end subroutine compute_gradient_shots

    !
    !> Process gradient
    !
    subroutine process_gradient

        integer :: i

        do i = 1, nmodel

            select case (model_name(i))

                case ('mt', 'stf')

                    if (rankid == 0) then
                        call warn(' Model value range = '//num2str(minval(model_grad(i)%array), '(es)') &
                            //', '//num2str(maxval(model_grad(i)%array), '(es)'))
                    end if

                case default

                    if(yn_shared_model_processing) then
                        call process_model_single_parameter(model_grad(i)%array, 'grad', param_name=model_m(i)%name)
                    else
                        call process_model_single_parameter(model_grad(i)%array, 'grad_'//tidy(model_name(i)), param_name=model_m(i)%name)
                    end if

            end select

        end do

    end subroutine process_gradient

    !
    !> Output gradient
    !
    subroutine output_gradient

        integer :: i

        if (rankid == 0) then

            do i = 1, nmodel
                call output_array(model_grad(i)%array, dir_iter_model(iter)//'/grad_'//tidy(model_name(i))//'.bin')
            end do

            call warn(date_time_compact()//' >>>>>>>>>> Gradient saved. ')

        end if

    end subroutine output_gradient

    !
    !> Processing gradient associated with a source
    !
    subroutine process_model_single_shot(ishot, w, name, file_w)

        integer, intent(in) :: ishot
        real, dimension(:, :), intent(inout) :: w
        character(len=*), intent(in) :: name, file_w

        integer :: i, j
        character(len=1024) :: shot_prefix, dir_mask, file_mask
        type(grid2) :: grd
        real, allocatable, dimension(:, :) :: andf_aux, andf_coh
        real :: shot_w_moving_balancex, shot_w_moving_balancez
        real :: shot_w_median_filtx, shot_w_median_filtz
        real, allocatable, dimension(:) :: shot_w_taperx, shot_w_taperz
        type(andf_param) :: param
        real :: shot_w_smoothx, shot_w_smoothz
        real :: shot_w_adaptive_mutex, shot_w_adaptive_mutez
        real :: recmin, recmax, srcmin, srcmax
        integer :: lb, ub, px, l
        real, allocatable, dimension(:, :) :: wmask
        real, allocatable, dimension(:) :: st, tp
        real, allocatable, dimension(:) :: wavenums, wamps, fkdips, fkdipamps
        character(len=32), allocatable, dimension(:) :: process_shot_w
        real :: cone_mutex, cone_mutez, cone_mutepower, cone_mutetaper
        real :: srcx, srcz, depth
        real :: mdx, mdz
        character(len=1024) :: file_andf_aux, file_andf_coh

        call readpar_nstring(file_parameter, 'process_'//'shot_'//tidy(name), process_shot_w, [''])

        shot_prefix = tidy(dir_scratch)//'/shot_'//num2str(set_srcid(ishot))
        call grd%input(file_w)
        mdx = grd%d2
        mdz = grd%d1

        do i = 1, size(process_shot_w)

            select case (process_shot_w(i))

                case ('smooth')
                    ! Gaussian smooth
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smooth_x', shot_w_smoothx, 3.0*mdx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smooth_z', shot_w_smoothz, 3.0*mdz)
                    grd%array = gauss_filt(grd%array, [shot_w_smoothz/mdz, shot_w_smoothx/mdx])

                case ('max_balance')
                    if (maxval(grd%array) /= 0) then
                        grd%array = grd%array/maxval(grd%array)
                    end if

                case ('rms_balance')
                    ! Normalize with shot image energy
                    grd%array = grd%array/mean(grd%array, 2)

                case ('moving_balance')
                    ! Moving balance
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_moving_balance_x', shot_w_moving_balancex, 3*mdx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_moving_balance_z', shot_w_moving_balancez, 3*mdz)
                    grd%array = balance_filt(grd%array, nint([0.5*shot_w_moving_balancez/mdz, 0.5*shot_w_moving_balancex/mdx]), 0.01)

                case ('median_filt')
                    ! Median filtering
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_median_filt_x', shot_w_median_filtx, mdx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_median_filt_z', shot_w_median_filtz, mdz)
                    grd%array = median_filt(grd%array, nint([shot_w_median_filtz/mdz, shot_w_median_filtx/mdx]))

                case ('dip_filt')
                    ! Dip filtering
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dip_filt_zx', fkdips, [-100.0, 0.0, 100.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dip_filt_zx_coefs', fkdipamps, [0.0, 0.0, 0.0])
                    if (sum(abs(fkdipamps)) > 0) then
                        grd%array = dip_filt(grd%array, [1.0, dx/dz], fkdips, fkdipamps)
                    end if

                case ('remove_nan')
                    ! Remove NaN
                    grd%array = return_normal(grd%array)

                case ('andf_filt')
                    ! Structure-oriented nonlinear anisotropic diffusion filtering
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_smooth_x', param%smooth2, 2*mdx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_smooth_z', param%smooth1, 8*mdz)
                    param%smooth2 = param%smooth2/mdx
                    param%smooth1 = param%smooth1/mdz
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_powerm', param%powerm, 1.0)
                    call readpar_int(file_parameter, 'shot_'//tidy(name)//'_andf_t', param%niter, 5)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_sigma', param%sigma, 6*max(mdx, mdz))
                    param%sigma = param%sigma/max(mdx, mdz)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_alpha', param%lambda1, 0.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_beta', param%lambda2, 1.0)
                    call readpar_string(file_parameter, 'shot_'//tidy(name)//'_andf_aux', file_andf_aux, '')
                    call readpar_string(file_parameter, 'shot_'//tidy(name)//'_andf_coh', file_andf_coh, '')
                    if (file_andf_aux == '' .and. file_andf_coh == '') then
                        grd%array = andf_filt(grd%array, param)
                    else if (file_andf_aux /= '' .and. file_andf_coh == '') then
                        call prepare_model_single_parameter(andf_aux, 'andf_aux', file_andf_aux, update=.false.)
                        call alloc_array(andf_aux, [1, shot_nz, 1, shot_nx], &
                            source=andf_aux(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend))
                        grd%array = andf_filt(grd%array, param, aux=andf_aux)
                    else if (file_andf_aux == '' .and. file_andf_coh /= '') then
                        call prepare_model_single_parameter(andf_coh, 'andf_coh', file_andf_coh, update=.false.)
                        call alloc_array(andf_coh, [1, shot_nz, 1, shot_nx], &
                            source=andf_coh(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend))
                        grd%array = andf_filt(grd%array, param, acoh=andf_coh)
                    else if (file_andf_aux /= '' .and. file_andf_coh /= '') then
                        call prepare_model_single_parameter(andf_aux, 'andf_aux', file_andf_aux, update=.false.)
                        call alloc_array(andf_aux, [1, shot_nz, 1, shot_nx], &
                            source=andf_aux(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend))
                        call prepare_model_single_parameter(andf_coh, 'andf_coh', file_andf_coh, update=.false.)
                        call alloc_array(andf_coh, [1, shot_nz, 1, shot_nx], &
                            source=andf_coh(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend))
                        grd%array = andf_filt(grd%array, param, aux=andf_aux, acoh=andf_coh)
                    end if

                case ('wavenumber_filt')
                    ! Wavenumber-domain filtering in x-axis
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumber_filt_x', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumber_filt_x_coefs', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        !$omp parallel do private(j)
                        do j = 1, size(grd%array, 1)
                            grd%array(j, :) = fourier_filt(grd%array(j, :), dx, wavenums, wamps)
                        end do
                        !$omp end parallel do
                    end if

                    ! Wavenumber-domain filtering in z-axis
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumber_filt_z', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumber_filt_z_coefs', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        !$omp parallel do private(j)
                        do j = 1, size(grd%array, 2)
                            grd%array(:, j) = fourier_filt(grd%array(:, j), dz, wavenums, wamps)
                        end do
                        !$omp end parallel do
                    end if

                case ('taper')
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taper_x', shot_w_taperx, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taper_z', shot_w_taperz, [0.0, 0.0])
                    if (size(shot_w_taperx) == 1) then
                        call alloc_array(shot_w_taperx, [1, 2], source=[shot_w_taperx(1), shot_w_taperx(1)])
                    end if
                    if (size(shot_w_taperz) == 1) then
                        call alloc_array(shot_w_taperz, [1, 2], source=[shot_w_taperz(1), shot_w_taperz(1)])
                    end if
                    grd%array = taper(grd%array, nint([shot_w_taperz/mdz, shot_w_taperx/mdx]), &
                        ['blackman', 'blackman', 'blackman', 'blackman'])

                case ('mask')
                    ! Masking
                    call readpar_string(file_parameter, 'dir_shot_'//tidy(name)//'_mask', dir_mask, '')
                    if (dir_mask == '') then
                        call readpar_string(file_parameter, 'shot_'//tidy(name)//'_mask', file_mask, '')
                    else
                        file_mask = tidy(dir_mask)//'/'//tidy(shot_prefix)//'_mask.bin'
                    end if
                    call prepare_model_single_parameter(wmask, 'mask', file_mask, update=.false.)
                    call alloc_array(wmask, [1, shot_nz, 1, shot_nx], &
                        source=wmask(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend))
                    grd%array = mask(grd%array, wmask)

                case ('adaptive_mute')
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_adaptive_mute_x', shot_w_adaptive_mutex, -1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_adaptive_mute_z', shot_w_adaptive_mutez, -1.0)
                    if (shot_w_adaptive_mutex >= 0) then
                        ! find source-receiver widest possible range
                        recmin = minval(gmtr(ishot)%recr(:)%x - grd%o2)
                        recmax = maxval(gmtr(ishot)%recr(:)%x - grd%o2)
                        srcmin = minval(gmtr(ishot)%srcr(:)%x - grd%o2)
                        srcmax = maxval(gmtr(ishot)%srcr(:)%x - grd%o2)
                        ! ... and their integer grid point positions in the computed image
                        lb = nint(min(recmin, srcmin)/dx + 1)
                        ub = nint(max(recmax, srcmax)/dx + 1)
                        ! the taper length
                        px = nint(shot_w_adaptive_mutex/dx)
                        ! create the taper
                        call alloc_array(st, [lb, ub], pad=px)
                        st = 1.0
                        st = taper(st, [px, px], ['blackman', 'blackman'])
                        ! put the taper in the whole x range, which can be longer than the taper
                        call alloc_array(tp, [1, grd%n2])
                        do l = lb - px, ub + px
                            if (l >= 1 .and. l <= grd%n2) then
                                tp(l) = st(l)
                            end if
                        end do
                        ! now the taper has the same length with the image, and do the tapering
                        ! the resulting tapered image is now restricted to the region cropped by
                        ! the largest possible source-receiver offset
                        do l = 1, grd%n1
                            grd%array(l, :) = grd%array(l, :)*tp
                        end do
                    end if

                case ('cone_mute')
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_cone_mute_x', cone_mutex, -1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_cone_mute_z', cone_mutez, -1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_cone_mute_power', cone_mutepower, 2.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_cone_mute_taper', cone_mutetaper, 10*grd%d2)
                    if (cone_mutex >= 0 .and. cone_mutez >= 0) then

                        cone_mutepower = max(1.0, cone_mutepower)
                        srcx = mean(gmtr(ishot)%srcr(:)%x - grd%o2)
                        srcz = mean(gmtr(ishot)%srcr(:)%z)
                        call alloc_array(tp, [1, grd%n2])

                        do j = 1, grd%n1

                            depth = (j - 1)*grd%d1 + grd%o1

                            if (depth < srcz) then
                                ! When the depth is shallower than source then set zero
                                grd%array(j, :) = 0.0
                            else
                                ! Depth is deeper than the source
                                tp = 0.0
                                ! lower and upper spatial range
                                px = cone_mutex*(min(depth - srcz, cone_mutez)/cone_mutez)**(1.0/cone_mutepower)
                                lb = nint((srcx - px)/grd%d2 + 1)
                                ub = nint((srcx + px)/grd%d2 + 1)
                                ! the taper length
                                px = nint(cone_mutetaper/grd%d2)
                                ! create the taper
                                call alloc_array(st, [lb, ub], pad=px)
                                st = 1.0
                                st = taper(st, [px, px], ['blackman', 'blackman'])
                                ! put the taper in the whole x range, which can be longer than the taper
                                do l = lb - px, ub + px
                                    if (l >= 1 .and. l <= grd%n2) then
                                        tp(l) = st(l)
                                    end if
                                end do
                                grd%array(j, :) = grd%array(j, :)*tp
                            end if

                        end do

                    end if

            end select

            if (process_shot_w(i) /= '') then
                call warn(' Model value range = '//num2str(minval(grd%array), '(es)')//', '//num2str(maxval(grd%array), '(es)'))
                call warn(date_time_compact()//' Shot '//num2str(set_srcid(ishot))//' '//tidy(name) &
                    //' processing ('//tidy(process_shot_w(i))//') completed. ')
            end if

        end do

        ! Merge w
        w(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend) = &
            w(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend) + grd%array

        call warn(date_time_compact()//' Shot '//num2str(set_srcid(ishot))//' '//tidy(name)//' merged.')

    end subroutine

    !
    ! Processing model associated with a single parameter
    !
    !> Smoothing gradients can remove high-wavenumber noises
    !> generated during gradient calculations. This functiionality is
    !> not as fancy as anisotropic diffusion denoising, but
    !> is much more less computational expensive.
    !
    subroutine process_model_single_parameter(w, name, param_name)

        real, dimension(:, :), intent(inout) :: w
        character(len=*), intent(in) :: name
        character(len=*), intent(in), optional :: param_name

        integer :: i, j
        character(len=1024) :: file_mask
        real, allocatable, dimension(:) :: w_taperx, w_taperz
        real :: w_moving_balancex, w_moving_balancez
        real :: w_median_filtx, w_median_filtz
        real :: w_smoothx, w_smoothz
        type(andf_param) :: param
        real, allocatable, dimension(:, :) :: andf_aux, andf_coh, wmask
        character(len=32), allocatable, dimension(:) :: process_w
        real :: w_scalar
        real :: w_rms_balance_x, w_rms_balance_z
        integer :: wrx, wrz
        real, allocatable, dimension(:, :) :: wt
        integer, allocatable, dimension(:) :: update_iter
        character(len=1024) :: file_andf_aux, file_andf_coh

        call readpar_nstring(file_parameter, 'process_'//tidy(name), process_w, [''])

        if (present(param_name)) then
            call readpar_nint(file_parameter, tidy(param_name)//'_update_iter', update_iter, [1, niter_max])
            if (size(update_iter) == 1) then
                update_iter = [update_iter(1), niter_max]
            end if
            if (iter < update_iter(1) .or. iter > update_iter(2)) then
                w = 0
                return
            end if
        end if

        ! Process gradient
        do i = 1, size(process_w)

            if (rankid == 0) then
                call warn(' Model value range before processing = '//num2str(minval(w), '(es)')//', '//num2str(maxval(w), '(es)'))
            end if

            select case (process_w(i))

                case ('scale')
                    call readpar_xfloat(file_parameter, tidy(name)//'_scale', w_scalar, 1.0, iter*1.0)
                    w = w*w_scalar

                case ('max_balance')
                    if (maxval(w) /= 0) then
                        w = w/maxval(w)
                    end if

                case ('rms_balance')
                    if (mean(w, 2) /= 0) then
                        w = w/mean(w, 2)
                    end if

                case ('rms_balance_x')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rms_balance_x', w_rms_balance_x, 1.0*dx, iter*1.0)
                    wrx = nint(w_rms_balance_x/dx)
                    if (mod(wrx, 2) == 0) then
                        wrx = wrx - 1
                    end if
                    wt = w
                    call pad_array(wt, [0, 0, wrx, wrx])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 2)
                        w(:, j) = w(:, j)/norm2(wt(:, j - (wrx - 1)/2:j + (wrx - 1)/2))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('rms_balance_z')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rms_balance_z', w_rms_balance_z, 1.0*dz, iter*1.0)
                    wrz = nint(w_rms_balance_z/dz)
                    if (mod(wrz, 2) == 0) then
                        wrz = wrz - 1
                    end if
                    wt = w
                    call pad_array(wt, [wrz, wrz, 0, 0])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 1)
                        w(j, :) = w(j, :)/norm2(wt(j - (wrz - 1)/2:j + (wrz - 1)/2, :))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('moving_balance')
                    call readpar_xfloat(file_parameter, tidy(name)//'_moving_balance_x', w_moving_balancex, 6*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_moving_balance_z', w_moving_balancez, 6*dz, iter*1.0)
                    w = balance_filt(w, nint([0.5*w_moving_balancez/dz, 0.5*w_moving_balancex/dx]), 0.01)

                case ('taper')
                    call readpar_nfloat(file_parameter, tidy(name)//'_taper_x', w_taperx, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, tidy(name)//'_taper_z', w_taperz, [0.0, 0.0])
                    if (size(w_taperx) == 1) then
                        call alloc_array(w_taperx, [1, 2], source=[w_taperx(1), w_taperx(1)])
                    end if
                    if (size(w_taperz) == 1) then
                        call alloc_array(w_taperz, [1, 2], source=[w_taperz(1), w_taperz(1)])
                    end if
                    w = taper(w, nint([w_taperz/dz, w_taperx/dx]), ['blackman', 'blackman', 'blackman', 'blackman'])

                case ('smooth')
                    call readpar_xfloat(file_parameter, tidy(name)//'_smooth_x', w_smoothx, 3*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_smooth_z', w_smoothz, 3*dz, iter*1.0)
                    w = gauss_filt(w, [w_smoothz/dz, w_smoothx/dx])

                case ('andf_filt')
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smooth_x', param%smooth2, 2.0*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smooth_z', param%smooth1, 8.0*dz, iter*1.0)
                    param%smooth2 = param%smooth2/dx
                    param%smooth1 = param%smooth1/dz
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_powerm', param%powerm, 1.0, iter*1.0)
                    call readpar_xint(file_parameter, tidy(name)//'_andf_t', param%niter, 5, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_sigma', param%sigma, 6.0*max(dx, dz), iter*1.0)
                    param%sigma = param%sigma/max(dx, dz)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_alpha', param%lambda1, 1.0e-3, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_beta', param%lambda2, 1.0, iter*1.0)
                    call readpar_xstring(file_parameter, tidy(name)//'_andf_aux', file_andf_aux, '', iter*1.0)
                    call readpar_xstring(file_parameter, tidy(name)//'_andf_coh', file_andf_coh, '', iter*1.0)
                    if (file_andf_aux == '' .and. file_andf_coh == '') then
                        w = andf_filt(w, param)
                    else if (file_andf_aux /= '' .and. file_andf_coh == '') then
                        call prepare_model_single_parameter(andf_aux, 'andf_aux', file_andf_aux, update=.false.)
                        w = andf_filt(w, param, aux=andf_aux)
                    else if (file_andf_aux == '' .and. file_andf_coh /= '') then
                        call prepare_model_single_parameter(andf_coh, 'andf_coh', file_andf_coh, update=.false.)
                        w = andf_filt(w, param, acoh=andf_coh)
                    else if (file_andf_aux /= '' .and. file_andf_coh /= '') then
                        call prepare_model_single_parameter(andf_aux, 'andf_aux', file_andf_aux, update=.false.)
                        call prepare_model_single_parameter(andf_coh, 'andf_coh', file_andf_coh, update=.false.)
                        w = andf_filt(w, param, aux=andf_aux, acoh=andf_coh)
                    end if

                case ('median_filt')
                    call readpar_xfloat(file_parameter, tidy(name)//'_median_filt_x', w_median_filtx, dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_median_filt_z', w_median_filtz, dz, iter*1.0)
                    w = median_filt(w, nint([w_median_filtz/dz, w_median_filtx/dx]))

                case ('mask')
                    call readpar_xstring(file_parameter, tidy(name)//'_mask', file_mask, file_mask, iter*1.0)
                    call prepare_model_single_parameter(wmask, 'mask', file_mask, update=.false.)
                    w = mask(w, wmask)

            end select

            if (rankid == 0) then
                call warn(' Model value range after processing = '//num2str(minval(w), '(es)')//', '//num2str(maxval(w), '(es)'))
                if (process_w(i) /= '') then
                    call warn(date_time_compact()//' '//tidy(name)//' processing ('//tidy(process_w(i))//') completed. ')
                end if
            end if

        end do

    end subroutine

end module
