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
    use acoustic_iso_3d
    use elastic_vhtiort_3d
    use elastic_tti_3d

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
                    model_grad(i)%array = zeros(nz, ny, nx)
                case ('mt')
                    model_grad(i)%array = zeros(nc_mt, ns, 1)
            end select
        end do

    end subroutine zero_gradient

    !
    !> Compute gradients shot by shot and sum
    !
    subroutine compute_gradient_shots

        character(len=1024) :: dir_field, dir_from, dir_to
        integer :: i
        type(grid3) :: grd
        type(wave_solver_acoustic_iso_3d) :: solver_acoustic_iso
        type(wave_solver_elastic_vhtiort_3d) :: solver_elastic_vhtiort
        type(wave_solver_elastic_tti_3d) :: solver_elastic_tti

        ! temporary directory
        dir_scratch = tidy(dir_working)//'/scratch'
        dir_synthetic = dir_iter_synthetic(iter)
        dir_synthetic_processed = dir_iter_synthetic_processed(iter)
        dir_field = dir_iter_record(iter)
        dir_adjoint = dir_iter_adjoint_source(iter)
        if (rankid == 0) then
            call make_directory(dir_scratch)
        end if
        call mpibarrier

        ! step misfit set to zero
        step_misfit = 0.0d0

        ! Enforce Vp/Vs ratio in an appropriate range
        call clip_vpvsratio

        ! Show the medium parameter statistics
        if (rankid == 0) then
            do i = 1, nmodel
                call plot_histogram(model_m(i)%array, &
                    label=date_time_compact()//' '//tidy(model_name(i))//' distribution ')
            end do
            do i = 1, nmodel_aux
                call plot_histogram(model_aux(i)%array, &
                    label=date_time_compact()//' '//tidy(model_name_aux(i))//' distribution ')
            end do
        end if

        call mpibarrier

        ! FWI gradient computation, shot paralellization with hybrid OpenMP + MPI for each shot
        do ishot = shot_in_group(groupid, 1), shot_in_group(groupid, 2)

            shot_prefix = tidy(dir_scratch)//'/shot_'//num2str(set_srcid(ishot))

            call set_adaptive_range(gmtr(ishot))

            select case (which_medium)

                case ('acoustic-iso')

                    solver_acoustic_iso%nx = shot_nx
                    solver_acoustic_iso%ny = shot_ny
                    solver_acoustic_iso%nz = shot_nz
                    solver_acoustic_iso%dx = dx
                    solver_acoustic_iso%dy = dy
                    solver_acoustic_iso%dz = dz
                    solver_acoustic_iso%ox = shot_xbeg
                    solver_acoustic_iso%oy = shot_ybeg
                    solver_acoustic_iso%oz = shot_zbeg
                    solver_acoustic_iso%dt = dt
                    solver_acoustic_iso%tmax = tmax
                    solver_acoustic_iso%data_dt = data_dt
                    solver_acoustic_iso%data_tmax = data_tmax
                    solver_acoustic_iso%dir_synthetic = tidy(dir_synthetic)
                    solver_acoustic_iso%dir_working = tidy(dir_scratch)
                    if (sum(snaps) > 0) then
                        solver_acoustic_iso%dir_snapshot = tidy(dir_snapshot)
                        solver_acoustic_iso%snaps = regspace(snaps(1), snaps(2), snaps(3))
                    end if
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
                    solver_elastic_vhtiort%ny = shot_ny
                    solver_elastic_vhtiort%nz = shot_nz
                    solver_elastic_vhtiort%dx = dx
                    solver_elastic_vhtiort%dy = dy
                    solver_elastic_vhtiort%dz = dz
                    solver_elastic_vhtiort%ox = shot_xbeg
                    solver_elastic_vhtiort%oy = shot_ybeg
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
                    solver_elastic_vhtiort%compy = yn_compy
                    solver_elastic_vhtiort%compz = yn_compz
                    solver_elastic_vhtiort%yn_update_medium = yn_update_medium
                    solver_elastic_vhtiort%yn_update_source = yn_update_source
                    solver_elastic_vhtiort%mt = flatten(slice(get_model('mt', 0.0), dim=2, index=ishot))

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
                            solver_elastic_vhtiort%tigam = get_model('gamma', 0.0)
                            solver_elastic_vhtiort%tithe = get_model('theta', 0.0)
                            solver_elastic_vhtiort%tiphi = get_model('phi', 0.0)
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                        case ('a-t')
                            solver_elastic_vhtiort%anisotropy_type = 'a-t'
                            solver_elastic_vhtiort%vp = get_model('vp')
                            solver_elastic_vhtiort%vs = get_model('vs')
                            solver_elastic_vhtiort%tieps = get_model('epsilon', 0.0)
                            solver_elastic_vhtiort%tieta = get_model('eta', 0.0)
                            solver_elastic_vhtiort%tigam = get_model('gamma', 0.0)
                            solver_elastic_vhtiort%tithe = get_model('theta', 0.0)
                            solver_elastic_vhtiort%tiphi = get_model('phi', 0.0)
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                        case ('cij')
                            solver_elastic_vhtiort%anisotropy_type = 'cij'
                            solver_elastic_vhtiort%c11 = get_model('c11')
                            solver_elastic_vhtiort%c12 = get_model('c12')
                            solver_elastic_vhtiort%c13 = get_model('c13')
                            solver_elastic_vhtiort%c22 = get_model('c22')
                            solver_elastic_vhtiort%c23 = get_model('c23')
                            solver_elastic_vhtiort%c33 = get_model('c33')
                            solver_elastic_vhtiort%c44 = get_model('c44')
                            solver_elastic_vhtiort%c55 = get_model('c55')
                            solver_elastic_vhtiort%c66 = get_model('c66')
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                    end select

                    call solver_elastic_vhtiort%forward

                case ('elastic-tti')
                    solver_elastic_tti%nx = shot_nx
                    solver_elastic_tti%ny = shot_ny
                    solver_elastic_tti%nz = shot_nz
                    solver_elastic_tti%dx = dx
                    solver_elastic_tti%dy = dy
                    solver_elastic_tti%dz = dz
                    solver_elastic_tti%ox = shot_xbeg
                    solver_elastic_tti%oy = shot_ybeg
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
                    solver_elastic_tti%compy = yn_compy
                    solver_elastic_tti%compz = yn_compz
                    solver_elastic_tti%yn_update_medium = yn_update_medium
                    solver_elastic_tti%yn_update_source = yn_update_source
                    solver_elastic_tti%mt = flatten(slice(get_model('mt', 0.0), dim=2, index=ishot))

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
                            solver_elastic_tti%tigam = get_model('gamma', 0.0)
                            solver_elastic_tti%tithe = get_model('theta', 0.0)
                            solver_elastic_tti%tiphi = get_model('phi', 0.0)
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                        case ('a-t')
                            solver_elastic_tti%anisotropy_type = 'a-t'
                            solver_elastic_tti%vp = get_model('vp')
                            solver_elastic_tti%vs = get_model('vs')
                            solver_elastic_tti%tieps = get_model('epsilon', 0.0)
                            solver_elastic_tti%tieta = get_model('eta', 0.0)
                            solver_elastic_tti%tigam = get_model('gamma', 0.0)
                            solver_elastic_tti%tithe = get_model('theta', 0.0)
                            solver_elastic_tti%tiphi = get_model('phi', 0.0)
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                        case ('cij')
                            solver_elastic_tti%anisotropy_type = 'cij'
                            solver_elastic_tti%c11 = get_model('c11')
                            solver_elastic_tti%c12 = get_model('c12')
                            solver_elastic_tti%c13 = get_model('c13')
                            solver_elastic_tti%c14 = get_model('c14', 0.0)
                            solver_elastic_tti%c15 = get_model('c15', 0.0)
                            solver_elastic_tti%c16 = get_model('c16', 0.0)
                            solver_elastic_tti%c22 = get_model('c22')
                            solver_elastic_tti%c23 = get_model('c23')
                            solver_elastic_tti%c24 = get_model('c24', 0.0)
                            solver_elastic_tti%c25 = get_model('c25', 0.0)
                            solver_elastic_tti%c26 = get_model('c26', 0.0)
                            solver_elastic_tti%c33 = get_model('c33')
                            solver_elastic_tti%c34 = get_model('c34', 0.0)
                            solver_elastic_tti%c35 = get_model('c35', 0.0)
                            solver_elastic_tti%c36 = get_model('c36', 0.0)
                            solver_elastic_tti%c44 = get_model('c44')
                            solver_elastic_tti%c45 = get_model('c45', 0.0)
                            solver_elastic_tti%c46 = get_model('c46', 0.0)
                            solver_elastic_tti%c55 = get_model('c55')
                            solver_elastic_tti%c56 = get_model('c56', 0.0)
                            solver_elastic_tti%c66 = get_model('c66')
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                    end select

                    call solver_elastic_tti%forward

            end select

            if (rankid_group == 0) then

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
                    call make_directory(tidy(dir_working)//'/iteration_0/model')
                    do i = 1, nmodel
                        call output_array(model_m(i)%array, tidy(dir_working)//'/iteration_0/model/' &
                            //tidy(model_name(i))//'.bin')
                    end do

                end if

            end if
            ! Note that here synthetic_processed must be broadcasted as
            ! data processing only occurs in rankid_group = 0
            call bcast_group(synthetic_processed, source=0)
            call mpibarrier_group

            ! Compute adjoint source
            if (synthetic_processed) then
                ! If the synthetic is processed
                call compute_adjoint_source(ishot, step_misfit(ishot), dir_synthetic_processed, dir_field, dir_adjoint)
            else
                ! If the synthetic is intact
                call compute_adjoint_source(ishot, step_misfit(ishot), dir_synthetic, dir_field, dir_adjoint)
            end if

            if (rankid_group == 0) then

                call warn(date_time_compact()//' >> Shot '//num2str(set_srcid(ishot)) &
                    //' misfit = '//num2str(step_misfit(ishot), '(es)'))

                if (yn_misfit_only) then
                    ! call execute_command_line('rm -rf '//tidy(shot_prefix)//'_*')
                    cycle
                end if

                call process_adjoint_source(ishot)

            else

                if (yn_misfit_only) then
                    cycle
                end if

            end if
            call mpibarrier_group

            select case (which_medium)

                case ('acoustic-iso')
                    solver_acoustic_iso%dir_adjoint = tidy(dir_adjoint)
                    solver_acoustic_iso%cc_step_interval = cc_step_interval
                    solver_acoustic_iso%energy_precond = yn_energy_precond
                    call solver_acoustic_iso%adjoint

                case ('elastic-iso', 'elastic-vhtiort')
                    solver_elastic_vhtiort%dir_adjoint = tidy(dir_adjoint)
                    solver_elastic_vhtiort%cc_step_interval = cc_step_interval
                    solver_elastic_vhtiort%energy_precond = yn_energy_precond
                    call solver_elastic_vhtiort%adjoint

                case ('elastic-tti')
                    solver_elastic_tti%dir_adjoint = tidy(dir_adjoint)
                    solver_elastic_tti%cc_step_interval = cc_step_interval
                    solver_elastic_tti%energy_precond = yn_energy_precond
                    call solver_elastic_tti%adjoint

            end select

            if (rankid_group == 0) then
                call warn(date_time_compact()//' >> Shot '//num2str(set_srcid(ishot)) &
                    //' gradient computation completed. ')
            end if

            ! Process computed gradients
            do i = 1, nmodel

                select case (model_name(i))

                    case ('mt', 'stf')

                        call grd%input(tidy(shot_prefix)//'_grad_'//tidy(model_name(i))//'.grd')
                        model_grad(i)%array(:, ishot, 1) = model_grad(i)%array(:, ishot, 1) + grd%array(:, 1, 1)
                        if (rankid_group == 0) then
                            call warn(date_time_compact()//' Shot '//num2str(set_srcid(ishot))//' '//tidy(model_name(i))//' merged.')
                        end if

                    case default

                        if (yn_shared_model_processing) then
                            call process_model_single_shot(ishot, model_grad(i)%array, 'grad', &
                                tidy(shot_prefix)//'_grad_'//tidy(model_name(i))//'.grd')
                        else
                            call process_model_single_shot(ishot, model_grad(i)%array, 'grad_'//tidy(model_name(i)), &
                                tidy(shot_prefix)//'_grad_'//tidy(model_name(i))//'.grd')
                        end if

                end select

            end do

            !            if (rankid_group == 0) then
            !                call execute_command_line('rm -rf '//tidy(shot_prefix)//'_*')
            !            end if
            call mpibarrier_group

        end do

        call mpibarrier

        ! collect misfit
        if (rankid_group /= 0) then
            step_misfit = 0
        end if
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

                    if (yn_shared_model_processing) then
                        call process_model_single_parameter(model_grad(i)%array, 'grad')
                    else
                        call process_model_single_parameter(model_grad(i)%array, 'grad_'//tidy(model_name(i)))
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
    !> Process gradient of one shot
    !
    subroutine process_model_single_shot(ishot, w, name, file_w)

        integer, intent(in) :: ishot
        real, dimension(:, :, :), intent(inout) :: w
        character(len=*), intent(in) :: name, file_w

        integer :: i, j, k
        character(len=1024) :: shot_prefix, dir_mask, file_mask
        type(grid3) :: grd
        real, allocatable, dimension(:, :, :) :: andfaux, andfcoh
        real :: shot_w_movingbalx, shot_w_movingbaly, shot_w_movingbalz
        real, allocatable, dimension(:) :: shot_w_taperx, shot_w_tapery, shot_w_taperz
        real :: shot_w_medianfiltx, shot_w_medianfilty, shot_w_medianfiltz
        type(andf_param) :: param
        real :: shot_w_smoothx, shot_w_smoothy, shot_w_smoothz
        real, allocatable, dimension(:) :: wavenums, wamps, fkdips, fkdipamps
        real, allocatable, dimension(:, :, :) :: wmask
        character(len=32), allocatable, dimension(:) :: process_shot_w
        real :: mdx, mdy, mdz
        character(len=1024) :: file_andf_aux, file_andf_coh

        call readpar_nstring(file_parameter, 'process_'//'shot_'//tidy(name)//'', process_shot_w, [''])

        ! shot prefix for convenience
        shot_prefix = tidy(dir_scratch)//'/shot_'//num2str(set_srcid(ishot))
        call grd%input(file_w)
        mdx = grd%d3
        mdy = grd%d2
        mdz = grd%d1

        do i = 1, size(process_shot_w)

            select case (process_shot_w(i))

                case ('smooth')
                    ! Gaussian smooth
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smoothx', shot_w_smoothx, 3*mdx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smoothy', shot_w_smoothy, 3*mdy)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smoothz', shot_w_smoothz, 3*mdz)
                    grd%array = gauss_filt(grd%array, [shot_w_smoothz/mdz, shot_w_smoothy/mdy, shot_w_smoothx/mdx])

                case ('max_balance')
                    if (maxval(grd%array) /= 0) then
                        grd%array = grd%array/maxval(grd%array)
                    end if

                case ('rms_balance')
                    ! Normalize with shot image energy
                    grd%array = grd%array/mean(grd%array, 2)

                case ('moving_balance')
                    ! Moving balance
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbalx', shot_w_movingbalx, 6*mdx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbaly', shot_w_movingbaly, 6*mdy)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbalz', shot_w_movingbalz, 6*mdz)
                    grd%array = balance_filt(grd%array, nint([0.5*shot_w_movingbalz/mdz, 0.5*shot_w_movingbaly/mdy, 0.5*shot_w_movingbalx/mdx]), 0.01)

                case ('median_filt')
                    ! Median filtering
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfiltx', shot_w_medianfiltx, mdx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfilty', shot_w_medianfilty, mdy)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfiltz', shot_w_medianfiltz, mdz)
                    grd%array = median_filt(grd%array, nint([shot_w_medianfiltz/mdz, shot_w_medianfilty/mdy, shot_w_medianfiltx/mdx]))

                case ('dip_filt')
                    ! Dip filtering
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltzx', fkdips, [-100.0, 0.0, 100.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltzx_amps', fkdipamps, [0.0, 0.0, 0.0])
                    if (sum(abs(fkdipamps)) > 0) then
                        grd%array = dip_filt(grd%array, [1.0, dy/dz, dx/dz], fkdips, fkdipamps, 13)
                    end if

                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltzy', fkdips, [-100.0, 0.0, 100.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltzy_amps', fkdipamps, [0.0, 0.0, 0.0])
                    if (sum(abs(fkdipamps)) > 0) then
                        grd%array = dip_filt(grd%array, [1.0, dy/dz, dx/dz], fkdips, fkdipamps, 12)
                    end if

                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltyx', fkdips, [-100.0, 0.0, 100.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_dipfiltyx_amps', fkdipamps, [0.0, 0.0, 0.0])
                    if (sum(abs(fkdipamps)) > 0) then
                        grd%array = dip_filt(grd%array, [1.0, dy/dz, dx/dz], fkdips, fkdipamps, 23)
                    end if

                case ('remove_nan')
                    ! Remove NaN
                    grd%array = return_normal(grd%array)

                case ('laplace_filt')
                    grd%array = laplace_filt(grd%array)

                case ('andf_filt')
                    ! Structure-oriented nonlinear anisotropic diffusion filtering
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_smoothx', param%smooth3, 2*mdx)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_smoothy', param%smooth2, 2*mdy)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_smoothz', param%smooth1, 8*mdz)
                    param%smooth3 = param%smooth3/mdx
                    param%smooth2 = param%smooth2/mdy
                    param%smooth1 = param%smooth1/mdz
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_powerm', param%powerm, 1.0)
                    call readpar_int(file_parameter, 'shot_'//tidy(name)//'_andf_t', param%niter, 5)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_sigma', param%sigma, 6*max(mdx, mdy, mdz))
                    param%sigma = param%sigma/max(mdx, mdy, mdz)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_alpha', param%lambda1, 1.0e-3)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_beta', param%lambda2, 1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_andf_gamma', param%lambda3, 1.0)
                    call readpar_string(file_parameter, 'shot_'//tidy(name)//'_andf_aux', file_andf_aux, '')
                    call readpar_string(file_parameter, 'shot_'//tidy(name)//'_andf_coh', file_andf_coh, '')
                    if (file_andf_aux == '' .and. file_andf_coh == '') then
                        grd%array = andf_filt(grd%array, param)
                    else if (file_andf_aux /= '' .and. file_andf_coh == '') then
                        call prepare_model_single_parameter(andfaux, 'andf_aux', file_andf_aux, update=.false.)
                        call alloc_array(andfaux, [1, shot_nz, 1, shot_ny, 1, shot_nx], &
                            source=andfaux(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend))
                        grd%array = andf_filt(grd%array, param, aux=andfaux)
                    else if (file_andf_aux == '' .and. file_andf_coh /= '') then
                        call prepare_model_single_parameter(andfcoh, 'andf_coh', file_andf_coh, update=.false.)
                        call alloc_array(andfcoh, [1, shot_nz, 1, shot_ny, 1, shot_nx], &
                            source=andfcoh(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend))
                        grd%array = andf_filt(grd%array, param, acoh=andfcoh)
                    else if (file_andf_aux /= '' .and. file_andf_coh /= '') then
                        call prepare_model_single_parameter(andfaux, 'andf_aux', file_andf_aux, update=.false.)
                        call alloc_array(andfaux, [1, shot_nz, 1, shot_ny, 1, shot_nx], &
                            source=andfaux(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend))
                        call prepare_model_single_parameter(andfcoh, 'andf_coh', file_andf_coh, update=.false.)
                        call alloc_array(andfcoh, [1, shot_nz, 1, shot_ny, 1, shot_nx], &
                            source=andfcoh(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend))
                        grd%array = andf_filt(grd%array, param, aux=andfaux, acoh=andfcoh)
                    end if

                case ('wavenumber_filt')
                    ! Wavenumber-domain filtering in x-axis
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumx', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumx_amps', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        !$omp parallel do private(j, k)
                        do k = 1, size(grd%array, 1)
                            do j = 1, size(grd%array, 2)
                                grd%array(k, j, :) = fourier_filt(grd%array(k, j, :), dx, wavenums, wamps)
                            end do
                        end do
                        !$omp end parallel do
                    end if

                    ! Wavenumber-domain filtering in x-axis
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumy', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumy_amps', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        !$omp parallel do private(j, k)
                        do k = 1, size(grd%array, 1)
                            do j = 1, size(grd%array, 3)
                                grd%array(k, :, j) = fourier_filt(grd%array(k, :, j), dy, wavenums, wamps)
                            end do
                        end do
                        !$omp end parallel do
                    end if

                    ! Wavenumber-domain filtering in z-axis
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumz', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenumz_amps', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        !$omp parallel do private(j, k)
                        do k = 1, size(grd%array, 2)
                            do j = 1, size(grd%array, 3)
                                grd%array(:, k, j) = fourier_filt(grd%array(:, k, j), dz, wavenums, wamps)
                            end do
                        end do
                        !$omp end parallel do
                    end if

                case ('mask')
                    ! Masking
                    call readpar_string(file_parameter, 'dir_shot_'//tidy(name)//'_mask', dir_mask, '')
                    if (dir_mask == '') then
                        call readpar_string(file_parameter, 'shot_'//tidy(name)//'_mask', file_mask, '')
                    else
                        file_mask = tidy(dir_mask)//'/'//tidy(shot_prefix)//'_mask.bin'
                    end if
                    call prepare_model_single_parameter(wmask, 'mask', file_mask, update=.false.)
                    call alloc_array(wmask, [1, shot_nz, 1, shot_ny, 1, shot_nx], &
                        source=wmask(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend))
                    grd%array = mask(grd%array, wmask)

                case ('taper')
                    ! Tapering
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taperx', shot_w_taperx, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_tapery', shot_w_tapery, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taperz', shot_w_taperz, [0.0, 0.0])
                    if (size(shot_w_taperx) == 1) then
                        call alloc_array(shot_w_taperx, [1, 2], source=[shot_w_taperx(1), shot_w_taperx(1)])
                    end if
                    if (size(shot_w_tapery) == 1) then
                        call alloc_array(shot_w_tapery, [1, 2], source=[shot_w_tapery(1), shot_w_tapery(1)])
                    end if
                    if (size(shot_w_taperz) == 1) then
                        call alloc_array(shot_w_taperz, [1, 2], source=[shot_w_taperz(1), shot_w_taperz(1)])
                    end if
                    grd%array = taper(grd%array, nint([shot_w_taperz/mdz, shot_w_tapery/mdy, shot_w_taperx/mdx]), &
                        ['blackman', 'blackman', 'blackman', 'blackman', 'blackman', 'blackman'])

            end select

            if (process_shot_w(i) /= '') then
                call warn(' Model value range = '//num2str(minval(grd%array), '(es)')//', '//num2str(maxval(grd%array), '(es)'))
                call warn(date_time_compact()//' shot '//num2str(set_srcid(ishot))//' '//tidy(name) &
                    //' processing ('//tidy(process_shot_w(i))//') finished. ')
            end if

        end do

        ! Merge image
        w(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend) = &
            w(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend) + grd%array

        if (rankid_group == 0) then
            call warn(date_time_compact()//' Shot '//num2str(set_srcid(ishot))//' '//tidy(name)//' merged')
        end if

    end subroutine

    !
    !> Processing gradients for a single model parameter
    !
    !> @note Smoothing gradients can remove high-wavenumber noises
    !> generated during gradient calculations. This functiionality is
    !> not as fancy as anisotropic diffusion denoising, but
    !> is much more less computational expensive.
    !
    subroutine process_model_single_parameter(w, name, param_name)

        real, dimension(:, :, :), intent(inout) :: w
        character(len=*), intent(in) :: name
        character(len=*), intent(in), optional :: param_name

        integer :: i, j, k
        character(len=1024) :: file_mask
        real, allocatable, dimension(:) :: w_taperx, w_tapery, w_taperz
        real :: w_movingbalx, w_movingbaly, w_movingbalz
        real :: w_medianfiltx, w_medianfilty, w_medianfiltz
        real :: w_smoothx, w_smoothy, w_smoothz
        type(andf_param) :: param
        real, allocatable, dimension(:, :, :) :: andfaux, andfcoh, wmask
        character(len=32), allocatable, dimension(:) :: process_w
        real :: w_scalar
        real, allocatable, dimension(:, :, :) :: wt
        integer :: wrx, wry, wrz
        real :: w_rmsbalx, w_rmsbaly, w_rmsbalz
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

            select case (process_w(i))

                case ('scale')
                    call readpar_xfloat(file_parameter, tidy(name)//'_scale', w_scalar, 1.0, iter*1.0)
                    w = w*w_scalar

                case ('moving_balance')
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbalx', w_movingbalx, 3*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbaly', w_movingbaly, 3*dy, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbalz', w_movingbalz, 3*dz, iter*1.0)
                    w = balance_filt(w, nint([0.5*w_movingbalz/dz, 0.5*w_movingbaly/dy, 0.5*w_movingbalx/dx]), 0.01)

                case ('taper')
                    call readpar_nfloat(file_parameter, tidy(name)//'_taperx', w_taperx, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, tidy(name)//'_tapery', w_tapery, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, tidy(name)//'_taperz', w_taperz, [0.0, 0.0])
                    if (size(w_taperx) == 1) then
                        call alloc_array(w_taperx, [1, 2], source=[w_taperx(1), w_taperx(1)])
                    end if
                    if (size(w_tapery) == 1) then
                        call alloc_array(w_tapery, [1, 2], source=[w_tapery(1), w_tapery(1)])
                    end if
                    if (size(w_taperz) == 1) then
                        call alloc_array(w_taperz, [1, 2], source=[w_taperz(1), w_taperz(1)])
                    end if
                    w = taper(w, nint([w_taperz/dz, w_tapery/dy, w_taperx/dx]), &
                        ['blackman', 'blackman', 'blackman', 'blackman', 'blackman', 'blackman'])

                case ('smooth')
                    call readpar_xfloat(file_parameter, tidy(name)//'_smoothx', w_smoothx, 3*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_smoothy', w_smoothy, 3*dy, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_smoothz', w_smoothz, 3*dz, iter*1.0)
                    w = gauss_filt(w, [w_smoothz/dz, w_smoothy/dy, w_smoothx/dx])

                case ('andf_filt')
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smoothx', param%smooth3, 2*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smoothy', param%smooth2, 2*dy, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smoothz', param%smooth1, 8*dz, iter*1.0)
                    param%smooth3 = param%smooth3/dx
                    param%smooth2 = param%smooth2/dy
                    param%smooth1 = param%smooth1/dz
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_powerm', param%powerm, 1.0, iter*1.0)
                    call readpar_xint(file_parameter, tidy(name)//'_andf_t', param%niter, 5, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_sigma', param%sigma, 6*max(dx, dy, dz), iter*1.0)
                    param%sigma = param%sigma/max(dx, dy, dz)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_alpha', param%lambda1, 1.0e-3, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_beta', param%lambda2, 1.0, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_gamma', param%lambda3, 1.0, iter*1.0)
                    call readpar_xstring(file_parameter, tidy(name)//'_andf_aux', file_andf_aux, '', iter*1.0)
                    call readpar_xstring(file_parameter, tidy(name)//'_andf_coh', file_andf_coh, '', iter*1.0)
                    call readpar_int(file_parameter, tidy(name)//'_andf_rankx', rank3, 1)
                    call readpar_int(file_parameter, tidy(name)//'_andf_ranky', rank2, 1)
                    call readpar_int(file_parameter, tidy(name)//'_andf_rankz', rank1, 1)
                    if (file_andf_aux == '' .and. file_andf_coh == '') then
                        w = andf_filt_mpi(w, param)
                    else if (file_andf_aux /= '' .and. file_andf_coh == '') then
                        call prepare_model_single_parameter(andfaux, 'andf_aux', file_andf_aux, update=.false.)
                        w = andf_filt_mpi(w, param, aux=andfaux)
                    else if (file_andf_aux == '' .and. file_andf_coh /= '') then
                        call prepare_model_single_parameter(andfcoh, 'andf_coh', file_andf_coh, update=.false.)
                        w = andf_filt_mpi(w, param, acoh=andfcoh)
                    else if (file_andf_aux /= '' .and. file_andf_coh /= '') then
                        call prepare_model_single_parameter(andfaux, 'andf_aux', file_andf_aux, update=.false.)
                        call prepare_model_single_parameter(andfcoh, 'andf_coh', file_andf_coh, update=.false.)
                        w = andf_filt_mpi(w, param, aux=andfaux, acoh=andfcoh)
                    end if

                case ('median_filt')
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfiltx', w_medianfiltx, dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfilty', w_medianfilty, dy, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfiltz', w_medianfiltz, dz, iter*1.0)
                    w = median_filt(w, nint([w_medianfiltz/dz, w_medianfilty/dy, w_medianfiltx/dx]))

                case ('rms_balance')
                    if (mean(w, 2) /= 0) then
                        w = w/mean(w, 2)
                    end if

                case ('rms_balance_x')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbalx', w_rmsbalx, 1.0*dx, iter*1.0)
                    wrx = nint(w_rmsbalx/dx)
                    if (mod(wrx, 2) == 0) then
                        wrx = wrx - 1
                    end if
                    wt = w
                    call pad_array(wt, [0, 0, 0, 0, wrx, wrx])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 3)
                        w(:, :, j) = w(:, :, j)/norm2(wt(:, :, j - (wrx - 1)/2:j + (wrx - 1)/2))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('rms_balance_y')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbaly', w_rmsbaly, 1.0*dy, iter*1.0)
                    wry = nint(w_rmsbaly/dy)
                    if (mod(wry, 2) == 0) then
                        wry = wry - 1
                    end if
                    wt = w
                    call pad_array(wt, [0, 0, wry, wry, 0, 0])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 2)
                        w(:, j, :) = w(:, j, :)/norm2(wt(:, j - (wry - 1)/2:j + (wry - 1)/2, :))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('rms_balance_xy')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbalx', w_rmsbalx, 1.0*dx, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbaly', w_rmsbaly, 1.0*dy, iter*1.0)
                    wrx = nint(w_rmsbalx/dx)
                    if (mod(wrx, 2) == 0) then
                        wrx = wrx - 1
                    end if
                    wry = nint(w_rmsbaly/dy)
                    if (mod(wry, 2) == 0) then
                        wry = wry - 1
                    end if
                    wt = w
                    call pad_array(wt, [0, 0, wry, wry, wrx, wrx])
                    !$omp parallel do private(j, k)
                    do k = 1, size(w, 3)
                        do j = 1, size(w, 2)
                            w(:, j, k) = w(:, j, k)/norm2(wt(:, j - (wry - 1)/2:j + (wry - 1)/2, k - (wrx - 1)/2:k + (wrx - 1)/2))
                        end do
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('rms_balance_z')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbalz', w_rmsbalz, 1.0*dz, iter*1.0)
                    wrz = nint(w_rmsbalz/dz)
                    if (mod(wrz, 2) == 0) then
                        wrz = wrz - 1
                    end if
                    wt = w
                    call pad_array(wt, [wrz, wrz, 0, 0, 0, 0])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 1)
                        w(j, :, :) = w(j, :, :)/norm2(wt(j - (wrz - 1)/2:j + (wrz - 1)/2, :, :))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('mask')
                    call readpar_xstring(file_parameter, tidy(name)//'_mask', file_mask, file_mask, iter*1.0)
                    call prepare_model_single_parameter(wmask, 'mask', file_mask, update=.false.)
                    w = mask(w, wmask)

            end select

            if (rankid == 0) then
                if (process_w(i) /= '') then
                    call warn(' Model value range = '//num2str(minval(w), '(es)')//', '//num2str(maxval(w), '(es)'))
                    call warn(date_time_compact()//' '//tidy(name)//' processing ('//tidy(process_w(i))//') finished. ')
                end if
            end if

        end do

    end subroutine

end module gradient
