!
! © 2025-2026. Triad National Security, LLC. All rights reserved.
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


module inversion_step_size

    use vars
    use gradient
    use mod_parameters
    use mod_source_receiver
    use mod_utility
    use inversion_adjoint_source

    implicit none

#ifdef _dim2_
#define model_dimension dimension(:, :)
#endif

#ifdef _dim3_
#define model_dimension dimension(:, :, :)
#endif

    ! Maximum number of step size search
    ! In most cases, one (1) search is adequate
    integer :: nsearch_max = 6
    logical :: step_suitable

    public :: compute_step_size

contains

    !
    !> Check if a step size is suitable
    !
    subroutine check_step_range_single_model(step_model, &
            step_scale, srch, step_max)

        real, model_dimension, intent(in) :: srch
        real, intent(in) :: step_model, step_scale, step_max

        if (step_max == 0) then
            step_suitable = step_suitable .and. .true.
        else
            if (maxval(abs(step_model*srch*step_scale)) < step_max) then
                ! require all steps are within the prefined variation range
                step_suitable = step_suitable .and. .true.
            else
                ! require all steps are within the prefined variation range
                step_suitable = step_suitable .and. .false.
            end if
        end if

    end subroutine check_step_range_single_model

    !
    !> Compute the initial step size based on search direction values
    !
    subroutine initial_step_size(cij, srch, init_step, max_perturbation)

        ! arguments
        real, model_dimension, intent(in) :: cij, srch
        real, intent(inout) :: init_step
        real, intent(in) :: max_perturbation

        ! local variables
        real :: max_srch, max_cij

        max_cij = 0.0
        max_srch = 0.0

        ! Find out the maximum value of the search direciton
        max_srch = maxval(abs(srch))
        max_cij = maxval(abs(cij))

        ! Initial step size is maximum initial step size or 1% or maximum model value
        if (max_srch == 0) then
            init_step = 0.0
            return
        end if

        init_step = max_perturbation/max_srch

    end subroutine initial_step_size

    !
    !> Get min max of model purturbation
    !
    subroutine min_max_step_size(step, srch, step_scaling_factor, mname)

        real, model_dimension, intent(in) :: srch
        real, intent(in) :: step, step_scaling_factor
        character(len=*) :: mname

        real, allocatable, model_dimension :: deltam
        real :: min_deltam, max_deltam

        deltam = srch*step*step_scaling_factor

        ! Find out the maximum value of the search direciton
        min_deltam = minval(deltam)
        max_deltam = maxval(deltam)

        if (rankid == 0) then
            call warn(date_time_compact()//' Perturbation '//tidy(mname)//' range: ' &
                //num2str(min_deltam, '(es12.4)')//' ' &
                //num2str(max_deltam, '(es12.4)'))
        end if

    end subroutine min_max_step_size

    !
    !> Output updated model
    !
    subroutine output_updated_model

        integer :: i

        do i = 1, nmodel
            call output_array(model_m(i)%array, &
                dir_iter_model(iter)//'/updated_'//tidy(model_name(i))//'.bin')
        end do

    end subroutine output_updated_model

    !
    !> Check if search step size is within the desired range
    !
    subroutine check_step_range(s)

        real, intent(in) :: s
        integer :: i

        step_suitable = .true.

        do i = 1, nmodel
            call check_step_range_single_model(model_step(i), s, model_srch(i)%array, model_step_max(i))
        end do

    end subroutine check_step_range

    !
    !> Backup current model
    !
    subroutine backup_current_model

        integer :: i

        do i = 1, nmodel
            model_m_backup(i)%array = model_m(i)%array
        end do

    end subroutine backup_current_model

    !
    !> Restore current model
    !
    subroutine restore_current_model

        integer :: i

        do i = 1, nmodel
            model_m(i)%array = model_m_backup(i)%array
        end do

    end subroutine restore_current_model

    !
    !> Update model
    !
    subroutine update_model(step)

        real, intent(in) :: step
        integer :: i

        do i = 1, nmodel

            model_m(i)%array = model_m(i)%array + model_step(i)*model_srch(i)%array*step

            ! Box-clip model
            !            if (model_m(i)%name == 'vp' .or. model_m(i)%name == 'vs') then
            !                where (model_m(i)%array /= 0)
            !                    model_m(i)%array = clip(model_m(i)%array, model_min(i), model_max(i))
            !                end where
            !            else
            model_m(i)%array = clip(model_m(i)%array, model_min(i), model_max(i))
            !            end if

        end do

        ! Enforce Vp/Vs ratio in an appropriate range
        call clip_vpvsratio

    end subroutine update_model

    !
    !> Print perturbation information
    !
    subroutine print_perturb_info

        integer :: i

        do i = 1, nmodel
            call min_max_step_size(model_step(i), model_srch(i)%array, step_scaling_factor, model_name(i))
        end do

    end subroutine print_perturb_info

    !
    !> Calcualte initial step size for all parametes
    !
    !> Note that some parameters may be zero at the very beginning,
    !> then the initial step size is set to step_max_cij
    !> under such circumstances
    !
    subroutine compute_initial_step_size

        integer :: i
        real :: valstep

        do i = 1, nmodel

            select case (remove_string_after(model_name(i), ['c', 'C']))
                case ('vp', 'vs', 'rho')
                    valstep = 100.0
                case ('epsilon', 'delta', 'gamma', 'eta')
                    valstep = 0.1
                case ('theta', 'phi')
                    valstep = 0.1*const_pi_half
                case ('c', 'C')
                    valstep = 1.0e9
                case ('mt')
                    valstep = 1.0
            end select

            call readpar_xfloat(file_parameter, 'step_max_'//tidy(model_name(i)), model_step_max(i), valstep, iter*1.0)
            ! In most cases, step_max_scale_factor = 1
            model_step_max(i) = model_step_max(i)*step_max_scale_factor
            call initial_step_size(model_m(i)%array, model_srch(i)%array, model_step(i), model_step_max(i))

        end do

    end subroutine compute_initial_step_size

    !
    !> Calculate step size coefficients
    !
    subroutine compute_step_coef(srcindex, sum1, sum2)

        integer, intent(in) :: srcindex
        real, intent(inout) :: sum1, sum2

        integer :: i, nr, nn, nc, ic
        integer :: nt_base
        real :: dt_base
        real, allocatable, dimension(:) :: dobs, dsyn, dsyn_prev
        character(len=1024) :: filename, file_obs, file_syn, file_syn_prev
        type(su) :: seismo_obs, seismo_syn, seismo_syn_prev
        real, allocatable, dimension(:) :: err_trc_13, err_trc_23, weight
        real :: dd
        real, allocatable, dimension(:, :, :) :: data_obs, data_syn, data_syn_prev
        real, allocatable :: data1(:), data3(:, :, :)
        real, allocatable, dimension(:, :) :: dsyn2, dobs2, dsyn2_prev
        real, allocatable, dimension(:) :: norm_obs, norm_syn, norm_syn_prev
        real :: m_obs, m_syn, m_syn_prev
        real :: amp_obs, amp_syn, amp_syn_prev

        ! Get dimensions
        nc = size(data_name)
        nr = gmtr(set_gmtrid(srcindex))%nr

        ! Divide the traces into different group workers
        call alloc_array(trace_in_group_rank, [0, nrank_group - 1, 1, 2])
        call cut(1, nr, nrank_group, trace_in_group_rank)

        ! Set base nt, dt
        filename = '/shot_'//num2str(set_srcid(srcindex))//'_'//'seismogram_'//tidy(data_name(1))//'.su'

        if (synthetic_processed) then
            file_syn = dir_iter_synthetic_processed(iter)//tidy(filename)
        else
            file_syn = dir_iter_synthetic(iter)//tidy(filename)
        end if
        call seismo_syn%load(file_syn, nr=1)
        nt_base = seismo_syn%nt
        dt_base = seismo_syn%dt

        ! Get data
        data_obs = zeros(nt_base, nr, nc)
        data_syn = zeros(nt_base, nr, nc)
        data_syn_prev = zeros(nt_base, nr, nc)

        weight = ones(nr)
        do ic = 1, nc

            filename = '/shot_'//num2str(set_srcid(srcindex))//'_'//'seismogram_'//tidy(data_name(ic))//'.su'

            ! Observed data
            file_obs = dir_iter_record(iter)//tidy(filename)

            ! Synthetic data last iteration
            if (synthetic_processed) then
                file_syn_prev = dir_iter_synthetic_processed(iter - 1)//tidy(filename)
            else
                file_syn_prev = dir_iter_synthetic(iter - 1)//tidy(filename)
            end if

            ! Synthetic data current iteration
            if (synthetic_processed) then
                file_syn = dir_iter_synthetic_processed(iter)//tidy(filename)
            else
                file_syn = dir_iter_synthetic(iter)//tidy(filename)
            end if

            call seismo_syn_prev%load(file_syn_prev, nr=nr)
            call seismo_syn%load(file_syn, nr=nr)
            call seismo_obs%load(file_obs, nr=nr)

            norm_obs = zeros(nr)
            norm_syn = zeros(nr)
            norm_syn_prev = zeros(nr)
            !$omp parallel do private(i)
            do i = 1, nr
                norm_obs(i) = norm2(seismo_obs%trace(i)%data)
                norm_syn(i) = norm2(seismo_syn%trace(i)%data)
                norm_syn_prev(i) = norm2(seismo_syn_prev%trace(i)%data)
            end do
            !$omp end parallel do

            m_obs = median(norm_obs)
            m_syn = median(norm_syn)
            m_syn_prev = median(norm_syn_prev)

            !$omp parallel do private(i, amp_obs, amp_syn, amp_syn_prev) schedule(auto)
            do i = trace_in_group_rank(rankid_group, 1), trace_in_group_rank(rankid_group, 2)

                amp_obs = norm2(seismo_obs%trace(i)%data)
                amp_syn = norm2(seismo_syn%trace(i)%data)
                amp_syn_prev = norm2(seismo_syn_prev%trace(i)%data)

                if (gmtr(set_gmtrid(srcindex))%recr(i)%weight == 0 &
                        .or. amp_obs < trace_discard_threshold*m_obs &
                        .or. amp_syn < trace_discard_threshold*m_syn &
                        .or. amp_syn_prev < trace_discard_threshold*m_syn_prev &
                        .or. amp_obs*amp_syn*amp_syn_prev == 0) then

                    weight(i) = min(0.0, weight(i))

                else

                    ! Resample to make all data have consistent length
                    if (seismo_syn_prev%nt /= nt_base .or. seismo_syn_prev%dt /= dt_base) then
                        data_syn_prev(:, i, ic) = interp(seismo_syn_prev%trace(i)%data, seismo_syn_prev%nt, seismo_syn_prev%dt, 0.0, &
                            nt_base, dt_base, 0.0, 'cubic')
                    else
                        data_syn_prev(:, i, ic) = seismo_syn_prev%trace(i)%data
                    end if
                    if (seismo_syn%nt /= nt_base .or. seismo_syn%dt /= dt_base) then
                        data_syn(:, i, ic) = interp(seismo_syn%trace(i)%data, seismo_syn%nt, seismo_syn%dt, 0.0, &
                            nt_base, dt_base, 0.0, 'cubic')
                    else
                        data_syn(:, i, ic) = seismo_syn%trace(i)%data
                    end if
                    if (seismo_obs%nt /= nt_base .or. seismo_obs%dt /= dt_base) then
                        data_obs(:, i, ic) = interp(seismo_obs%trace(i)%data, seismo_obs%nt, seismo_obs%dt, 0.0, &
                            nt_base, dt_base, 0.0, 'cubic')
                    else
                        data_obs(:, i, ic) = seismo_obs%trace(i)%data
                    end if

                end if

            end do
            !$omp end parallel do

        end do

        call mpibarrier_group
        ! Each rank checked only its own traces; zero the other traces before the group sum
        weight(:trace_in_group_rank(rankid_group, 1) - 1) = 0.0
        weight(trace_in_group_rank(rankid_group, 2) + 1:) = 0.0
        call allreduce_array_group(weight)
        call allreduce_array_group(data_obs)
        call allreduce_array_group(data_syn)
        call allreduce_array_group(data_syn_prev)

        err_trc_13 = zeros(nr)
        err_trc_23 = zeros(nr)
        if (adj_nt > 0) then
            nn = adj_nt
            dd = (nt_base - 1.0)*dt_base/(nn - 1.0)
        else
            nn = nt_base
            dd = dt_base
        end if

        select case (misfit_type)

            case default
                ! Trace by trace

                data1 = zeros(1)

                do ic = 1, nc

                    err_trc_13 = 0.0
                    err_trc_23 = 0.0

                    !$omp parallel do private(i, data1, dobs, dsyn, dsyn_prev) schedule(dynamic)
                    do i = trace_in_group_rank(rankid_group, 1), trace_in_group_rank(rankid_group, 2)

                        if (weight(i) /= 0) then

                            ! Resample data to the desired length (adj_nt)
                            if (nt_base /= nn) then
                                dobs = interp_to(data_obs(:, i, ic), nn, method='cubic')
                                dsyn = interp_to(data_syn(:, i, ic), nn, method='cubic')
                                dsyn_prev = interp_to(data_syn_prev(:, i, ic), nn, method='cubic')
                            else
                                dobs = data_obs(:, i, ic)
                                dsyn = data_syn(:, i, ic)
                                dsyn_prev = data_syn_prev(:, i, ic)
                            end if

                            select case (misfit_type)

                                case ('waveform')
                                    call compute_adjsrc_waveform(dsyn, dsyn_prev, data1, err_trc_23(i))
                                    call compute_adjsrc_waveform(dobs, dsyn_prev, data1, err_trc_13(i))

                                case ('corr')
                                    call compute_adjsrc_corr(dsyn, dsyn_prev, data1, err_trc_23(i))
                                    call compute_adjsrc_corr(dobs, dsyn_prev, data1, err_trc_13(i))

                                case ('adaptive')
                                    call compute_adjsrc_adaptive(dsyn, dsyn_prev, dd, data1, err_trc_23(i))
                                    call compute_adjsrc_adaptive(dobs, dsyn_prev, dd, data1, err_trc_13(i))

                                case ('local-adaptive')
                                    call compute_adjsrc_local_adaptive(dsyn, dsyn_prev, dd, data1, err_trc_23(i))
                                    call compute_adjsrc_local_adaptive(dobs, dsyn_prev, dd, data1, err_trc_13(i))

                                case ('envelope')
                                    call compute_adjsrc_envelope(dsyn, dsyn_prev, data1, err_trc_23(i))
                                    call compute_adjsrc_envelope(dobs, dsyn_prev, data1, err_trc_13(i))

                                case ('phase')
                                    call compute_adjsrc_phase(dsyn, dsyn_prev, data1, err_trc_23(i))
                                    call compute_adjsrc_phase(dobs, dsyn_prev, data1, err_trc_13(i))

                                case default
                                    call compute_adjsrc_waveform(dsyn, dsyn_prev, data1, err_trc_23(i))
                                    call compute_adjsrc_waveform(dobs, dsyn_prev, data1, err_trc_13(i))

                            end select

                        end if

                    end do
                    !$omp end parallel do

                    call mpibarrier_group

                    call allreduce_array_group(err_trc_13)
                    call allreduce_array_group(err_trc_23)

                    sum1 = sum1 + sum(err_trc_23*err_trc_13)
                    sum2 = sum2 + sum(err_trc_23**2)

                end do

            case ('adaptive-spacetime')
                ! Group nearby traces for computing matching filter

                call readpar_xint(file_parameter, 'adaptive_half_window', adaptive_hw, 3, iter*1.0)

                do ic = 1, nc

                    data1 = zeros(1)

                    if (nt_base /= nn) then
                        dsyn2_prev = interp_to(data_syn_prev(:, :, ic), [nn, nr], method=['cubic', 'nearest'])
                        dsyn2 = interp_to(data_syn(:, :, ic), [nn, nr], method=['cubic', 'nearest'])
                        dobs2 = interp_to(data_obs(:, :, ic), [nn, nr], method=['cubic', 'nearest'])
                    else
                        dsyn2_prev = data_syn_prev(:, :, ic)
                        dsyn2 = data_syn(:, :, ic)
                        dobs2 = data_obs(:, :, ic)
                    end if

                    call alloc_array(dsyn2_prev, [1, nn, -adaptive_hw + 1, nr + adaptive_hw], &
                        source=pad(dsyn2_prev, [0, 0, adaptive_hw, adaptive_hw]))
                    call alloc_array(dsyn2, [1, nn, -adaptive_hw + 1, nr + adaptive_hw], &
                        source=pad(dsyn2, [0, 0, adaptive_hw, adaptive_hw]))
                    call alloc_array(dobs2, [1, nn, -adaptive_hw + 1, nr + adaptive_hw], &
                        source=pad(dobs2, [0, 0, adaptive_hw, adaptive_hw]))

                    !$omp parallel do private(i, dobs, dsyn, dsyn_prev) schedule(auto)
                    do i = trace_in_group_rank(rankid_group, 1), trace_in_group_rank(rankid_group, 2)

                        if (weight(i) /= 0) then

                            call compute_adjsrc_adaptive_spacetime(dsyn2(:, i - adaptive_hw:i + adaptive_hw), &
                                dsyn2_prev(:, i - adaptive_hw:i + adaptive_hw), dd, data1, err_trc_23(i))
                            call compute_adjsrc_adaptive_spacetime(dobs2(:, i - adaptive_hw:i + adaptive_hw), &
                                dsyn2_prev(:, i - adaptive_hw:i + adaptive_hw), dd, data1, err_trc_13(i))

                        end if

                    end do
                    !$omp end parallel do

                    call mpibarrier_group

                    call allreduce_array_group(err_trc_13)
                    call allreduce_array_group(err_trc_23)

                    sum1 = sum1 + sum(err_trc_23*err_trc_13)
                    sum2 = sum2 + sum(err_trc_23**2)

                end do

            case ('dtw')
                ! For DTW, the misfit is computed shot gather by shot gathter due to possible lateral smoothing

                data3 = zeros(1, nr, nc)

                do ic = 1, nc
                    call compute_adjsrc_dtw(data_syn(:, :, ic:ic), data_syn_prev(:, :, ic:ic), dt_base, weight, data3, err_trc_23)
                    call compute_adjsrc_dtw(data_obs(:, :, ic:ic), data_syn_prev(:, :, ic:ic), dt_base, weight, data3, err_trc_13)
                    sum1 = sum1 + sum(err_trc_23*err_trc_13)
                    sum2 = sum2 + sum(err_trc_23**2)
                end do

            case ('dtw-vector')
                ! For vector DTW, the misfit is computed shot gather by shot gather, but also all components together

                call compute_adjsrc_dtw(data_syn, data_syn_prev, dt_base, weight, data3, err_trc_23)
                call compute_adjsrc_dtw(data_obs, data_syn_prev, dt_base, weight, data3, err_trc_13)
                sum1 = sum1 + sum(err_trc_23*err_trc_13)
                sum2 = sum2 + sum(err_trc_23**2)

        end select

    end subroutine

    !
    !> Calculate the L2-norm error for a model
    !
    subroutine compute_step_misfit(step_scaling_factor, misfit)

        real, intent(in) :: step_scaling_factor
        real, intent(inout) :: misfit

        ! Update model
        call backup_current_model
        call update_model(step_scaling_factor)

        ! Forward modeling to get misfits
        yn_misfit_only = .true.
        yn_save_adjsrc = .false.
        yn_reconstruct = .false.
        ! yn_save_entire = .false.
        ! yn_encoded_synthetic_stage = .true.
        call compute_gradient_shots
        ! yn_encoded_synthetic_stage = .false.
        call mpibarrier

        misfit = sum(step_misfit)

        ! Backdate model
        call restore_current_model
        call mpibarrier

    end subroutine compute_step_misfit

    !
    !> Compute the optimal step size for a quasi-linear inversion
    !
    subroutine compute_step_size_linear

        ! Local variables
        real :: trial_misfit, sum1, sum2, data_misfit_current
        integer :: cnt, i
        character(len=1024) :: dir_from, dir_to

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Computing optimal step size')
        end if

        ! Calculate the initial step size
        call compute_initial_step_size

        ! default initial step scalar = 0.1
        step_scaling_factor = 0.2

        ! check if the initial step scalar is suitable
        step_suitable = .false.
        cnt = 1
        do while (.not. step_suitable .and. cnt < 20)
            call check_step_range(step_scaling_factor)
            step_scaling_factor = step_scaling_factor/2.0
            cnt = cnt + 1
        end do
        if (.not. step_suitable) then
            if (rankid == 0) then
                call warn(date_time_compact() &
                    //' Error: Cannot find suitable initial step size ')
            end if
            call mpibarrier
            call mpiend
        else
            if(rankid == 0) then
                call warn(date_time_compact() // ' Initial step scaling factor = ' &
                    // num2str(step_scaling_factor, '(es)'))
            end if
        end if

        ! Step size coefficients
        sum1 = 0.0
        sum2 = 0.0

        ! Update model
        call backup_current_model
        call update_model(step_scaling_factor)

        ! forward modeling and compute misfit
        yn_misfit_only = .true.
        yn_save_adjsrc = .false.
        yn_reconstruct = .false.
        ! yn_save_entire = .false.
        ! yn_encoded_synthetic_stage = .true.

        iter = iter + 1
        put_synthetic_in_scratch = .false.
        call make_iter_dir
        ! if (yn_source_encoding) then
        !    call encode_shot_gather
        ! end if
        call compute_gradient_shots
        ! yn_encoded_synthetic_stage = .false.

        call mpibarrier
        do i = shot_in_group(groupid, 1), shot_in_group(groupid, 2)
            call compute_step_coef(i, sum1, sum2)
        end do

        iter = iter - 1
        call mpibarrier
        if (rankid == 0) then
            call warn(date_time_compact()//' >> Step size coefficients computed ')
        end if

        ! Restore current iteration model
        call restore_current_model

        ! Calculate optimal step size; every rank of a shot group holds the same sums
        call mpibarrier
        if (rankid_group /= 0) then
            sum1 = 0.0
            sum2 = 0.0
        end if
        call mpi_allreduce(mpi_in_place, sum1, 1, mpi_real, mpi_sum, mpi_comm_world, mpi_ierr)
        call mpi_allreduce(mpi_in_place, sum2, 1, mpi_real, mpi_sum, mpi_comm_world, mpi_ierr)

        if (isnan(sum1) .or. isnan(sum2)) then
            call warn(date_time_compact()//' Step coefficient sum 1 = '//num2str(sum1, '(es)'))
            call warn(date_time_compact()//' Step coefficient sum 2 = '//num2str(sum2, '(es)'))
            call mpibarrier
            call mpiend
        else
            step_scaling_factor = sum1/(sum2 + float_tiny)
        end if
        if (rankid == 0) then
            call warn(date_time_compact()//' Step coefficient sum 1 = '//num2str(sum1, '(es)'))
            call warn(date_time_compact()//' Step coefficient sum 2 = '//num2str(sum2, '(es)'))
            call warn(date_time_compact()//' Step coefficient = '//num2str(step_scaling_factor, '(es)'))
        end if

        if (.not. (step_scaling_factor <= float_huge)) then
            if (rankid == 0) then
                call warn(date_time_compact()//' <compute_step_size_linear> Error: Step size NaN. Exiting. ')
                call mpi_abort(mpi_comm_world, mpi_err_other, mpi_ierr)
            end if
        end if

        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size scaling factor computed. ')
        end if

        ! Check step size
        step_scaling_factor = step_scaling_factor*2.0
        data_misfit_current = data_misfit(iter)
        trial_misfit = data_misfit_current*2.0
        cnt = 1

        ! Enforce update or not
        call readpar_xlogical(file_parameter, 'yn_enforce_update', yn_enforce_update, .false., iter*1.0)

        if (yn_enforce_update) then
            ! To enforce update regardless of misfit increase,
            ! use only one trial. This can be enabled after certain number of
            ! iterations, not necessarily from the begining

            do while (cnt < 2)

                step_scaling_factor = step_scaling_factor/2.0

                call check_step_range(step_scaling_factor)
                if (.not. step_suitable) then
                    cycle
                end if

                call compute_step_misfit(step_scaling_factor, trial_misfit)

                call print_step_info(cnt, step_scaling_factor, trial_misfit)

                cnt = cnt + 1

            end do

        else
            ! Otherwise, use conventional trial-error approach

            do while (trial_misfit > jumpout_factor*data_misfit_current .and. cnt < nsearch_max)

                step_scaling_factor = step_scaling_factor/2.0

                ! make sure step size and model values in suitable range
                ! if not in range, then half until suitable
                call check_step_range(step_scaling_factor)
                if (.not. step_suitable) then
                    cycle
                end if

                put_synthetic_in_scratch = .true.
                call compute_step_misfit(step_scaling_factor, trial_misfit)

                ! Print step size information
                call print_step_info(cnt, step_scaling_factor, trial_misfit)

                cnt = cnt + 1

            end do

            if (trial_misfit > jumpout_factor*data_misfit_current) then
                ! When searching the optimal step size, if the trial step size
                ! at the final trial cannot produce a smaller data misfit than
                ! that of the last iteration, then set step size to zero
                ! and the trial misfit to that in the last iteration
                step_scaling_factor = 0.0
                trial_misfit = data_misfit_current
            else
                ! Otherwise, copy the synthetic and synthetic processed data residing in scratch
                ! to the current synthetic and synthetic processed data
                ! so that they are the synthetic data of the current iteration
                if (rankid == 0) then
                    ! This only executed by rank0
                    put_synthetic_in_scratch = .true.
                    dir_from = dir_iter_synthetic(iter)
                    put_synthetic_in_scratch = .false.
                    dir_to = dir_iter_synthetic(iter)
                    call delete_directory(dir_to)
                    call move_directory(dir_from, dir_to)
                    if (synthetic_processed) then
                        put_synthetic_in_scratch = .true.
                        dir_from = dir_iter_synthetic_processed(iter)
                        put_synthetic_in_scratch = .false.
                        dir_to = dir_iter_synthetic_processed(iter)
                        call delete_directory(dir_to)
                        call move_directory(dir_from, dir_to)
                    end if
                end if
            end if

        end if

        ! Print perturbation infomation
        call print_perturb_info

        ! Current absolute data misfit
        data_misfit(iter) = trial_misfit

        ! Make all process ready
        call mpibarrier

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size computation completed.')
        end if

        ! Update model and output updated model for current iteration
        call update_model(step_scaling_factor)
        if (rankid == 0) then
            call output_updated_model
        end if
        put_synthetic_in_scratch = .false.

    end subroutine compute_step_size_linear

    !
    !> Calculate the optimal step size to update the model
    !>      with a quadratic + interval bounding+bisection line search method
    !
    subroutine compute_step_size_linesearch

        integer :: cnt
        real :: al, ar, ac, fl, fr, fc, a0
        real :: acr, arl, alc, bcr, brl, blc, num, den, an, fn, a2, f2
        real :: misfit0, step1, misfit1, step2, misfit2
        character(len=1024) :: dir_from, dir_to

        ! Print info
        if (rankid == 0) then
            call warn(' >>>>>>>>>> Computing step size ')
        end if

        ! Calculate the initial step size
        call compute_initial_step_size

        ! Default initial step scalar = 0.1
        step_scaling_factor = 1.0

        ! Check if the initial step scalar is suitable
        step_suitable = .false.
        cnt = 1
        do while (.not. step_suitable .and. cnt < 20)
            call check_step_range(step_scaling_factor)
            step_scaling_factor = step_scaling_factor/2.0
            cnt = cnt + 1
        end do
        if (.not. step_suitable) then
            if (rankid == 0) then
                call warn(date_time_compact()//' Error: Cannot find suitable initial step size ')
            end if
            call mpibarrier
            call mpiend
        else
            if (rankid == 0) then
                call warn(date_time_compact()//' Initial step scaling factor = ' &
                    //num2str(step_scaling_factor, '(es)'))
            end if
        end if

        a0 = step_scaling_factor

        ! The directory of synthetic
        put_synthetic_in_scratch = .true.

        ! Prepare the initial values
        al = 0.0
        fl = data_misfit(iter)
        misfit0 = fl
        call print_step_info(0, al, fl)
        call mpibarrier

        ar = a0
        call compute_step_misfit(ar, fr)
        call print_step_info(0, ar, fr)
        step1 = ar
        misfit1 = fr
        call mpibarrier

        ! If this step produces a smaller misfit,
        ! then skip trial
        if (fr < fl) then
            ac = ar
            fc = fr
            call mpibarrier
            goto 123
        end if

        ac = 0.5*(al + ar)
        call compute_step_misfit(ac, fc)
        call print_step_info(0, ac, fc)
        step2 = ac
        misfit2 = fc
        call mpibarrier

        ! If this step produces a smaller misfit,
        ! then skip trial
        if (fc < fl) then
            call mpibarrier
            goto 123
        end if

        cnt = 3

        call mpibarrier

        ! Quadratic fit + bisection line search
        do while (cnt < nsearch_max)

            if ((fc < fl) .and. (fc < fr)) then

                acr = ac - ar
                bcr = ac**2 - ar**2

                arl = ar - al
                brl = ar**2 - al**2

                alc = al - ac
                blc = al**2 - ac**2

                num = bcr*fl + brl*fc + blc*fr
                den = acr*fl + arl*fc + alc*fr

                if (den == 0.0) then
                    exit
                end if

                an = 0.5*num/den
                call compute_step_misfit(an, fn)
                call print_step_info(cnt - 2, an, fn)
                call mpibarrier

                cnt = cnt + 1

                if (an > ac) then
                    if (fn >= fc) then
                        ar = an
                        fr = fn
                    else
                        al = ac
                        fl = fc
                        ac = an
                        fc = fn
                    end if
                else
                    if (fn >= fc) then
                        al = an
                        fl = fn
                    else
                        ar = ac
                        fr = fc
                        ac = an
                        fc = fn
                    end if
                end if

            else
                ! Bisection safe switchover

                a2 = ac + 1.0e-1*(ar - al)
                call compute_step_misfit(a2, f2)
                call print_step_info(cnt - 2, a2, f2)
                call mpibarrier

                cnt = cnt + 1

                if (fc < f2) then
                    ar = a2
                    fr = f2
                else
                    al = ac
                    fl = fc
                end if

                ac = 0.5*(al + ar)
                call compute_step_misfit(ac, fc)
                call print_step_info(cnt - 2, ac, fc)
                call mpibarrier

                cnt = cnt + 1

            end if

            ! Stop search if range too small
            if (abs(ar - al) <= 0.05) then
                if (fc >= misfit1) then
                    ac = step1
                    fc = misfit1
                end if
                if (fc >= misfit2) then
                    ac = step2
                    fc = misfit2
                end if
                exit
            end if

        end do

        ! Current relative data misfit
        123 continue

        if (fc >= jumpout_factor*misfit0) then
            ! When searching optimal step size, if the trial step size
            ! at the final trial cannot produce a smaller data misfit than
            ! that of the last iteration, then set step size to zero
            ! and the trial misfit to that in the last iteration
            step_scaling_factor = 0.0
            data_misfit(iter) = misfit0
        else
            ! Otherwise, copy the synthetic and synthetic processed data residing in scratch
            ! to the current synthetic and synthetic processed data
            ! so that they are the synthetic data of the current iteration
            step_scaling_factor = ac
            data_misfit(iter) = fc
            if (rankid == 0) then
                ! This only executed by rank0
                put_synthetic_in_scratch = .true.
                dir_from = dir_iter_synthetic(iter)
                put_synthetic_in_scratch = .false.
                dir_to = dir_iter_synthetic(iter)
                call delete_directory(dir_to)
                call move_directory(dir_from, dir_to)
                if (synthetic_processed) then
                    put_synthetic_in_scratch = .true.
                    dir_from = dir_iter_synthetic_processed(iter)
                    put_synthetic_in_scratch = .false.
                    dir_to = dir_iter_synthetic_processed(iter)
                    call delete_directory(dir_to)
                    call move_directory(dir_from, dir_to)
                end if
            end if
        end if

        ! Make all process ready
        call mpibarrier

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size calculation completed. ')
        end if

        ! Update model and output
        call update_model(step_scaling_factor)
        if (rankid == 0) then
            call output_updated_model
        end if

        put_synthetic_in_scratch = .false.

    end subroutine compute_step_size_linesearch

    subroutine compute_step_size_quadratic

        integer :: cnt
        real :: step0, misfit0, step1, misfit1, step2, misfit2
        character(len=1024) :: dir_from, dir_to

        ! Print info
        if (rankid == 0) then
            call warn(' >>>>>>>>>> Computing step size ')
        end if

        ! Calculate the initial step size
        call compute_initial_step_size

        ! default initial step scalar = 0.1
        step_scaling_factor = 1.0

        ! check if the initial step scalar is suitable
        step_suitable = .false.
        cnt = 1
        do while (.not. step_suitable .and. cnt < 20)
            call check_step_range(step_scaling_factor)
            step_scaling_factor = step_scaling_factor/2.0
            cnt = cnt + 1
        end do
        if (.not. step_suitable) then
            if (rankid == 0) then
                call warn(date_time_compact() &
                    //' <compute_step_size_quadratic> Error: Cannot find a suitable initial step size.')
            end if
            call mpibarrier
            call mpiend
        else
            if (rankid == 0) then
                call warn(date_time_compact()//' Initial step scaling factor = ' &
                    //num2str(step_scaling_factor, '(es)'))
            end if
        end if

        ! The directory of synthetic
        put_synthetic_in_scratch = .true.

        step0 = 0.0
        misfit0 = data_misfit(iter)
        call print_step_info(0, step0, misfit0)
        call mpibarrier

        step1 = 0.1*step_scaling_factor
        call compute_step_misfit(step1, misfit1)
        call print_step_info(1, step1, misfit1)
        call mpibarrier

        step2 = step_scaling_factor
        call compute_step_misfit(step2, misfit2)
        call print_step_info(2, step2, misfit2)
        call mpibarrier

        step_scaling_factor = &
            0.5*((misfit0 - misfit2)*step1**2 + (-misfit0 + misfit1)*step2**2) &
            /(-(misfit2*step1) + misfit0*(step1 - step2) + misfit1*step2)
        call compute_step_misfit(step_scaling_factor, data_misfit(iter))
        call print_step_info(3, step_scaling_factor, data_misfit(iter))
        call mpibarrier

        if (rankid == 0) then
            ! This only executed by rank0
            put_synthetic_in_scratch = .true.
            dir_from = dir_iter_synthetic(iter)
            put_synthetic_in_scratch = .false.
            dir_to = dir_iter_synthetic(iter)
            call delete_directory(dir_to)
            call move_directory(dir_from, dir_to)
            if (synthetic_processed) then
                put_synthetic_in_scratch = .true.
                dir_from = dir_iter_synthetic_processed(iter)
                put_synthetic_in_scratch = .false.
                dir_to = dir_iter_synthetic_processed(iter)
                call delete_directory(dir_to)
                call move_directory(dir_from, dir_to)
            end if
        end if

        ! Make all process ready
        call mpibarrier

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size calculation completed. ')
        end if

        ! update model and output
        call update_model(step_scaling_factor)
        if (rankid == 0) then
            call output_updated_model
        end if

        put_synthetic_in_scratch = .false.

    end subroutine compute_step_size_quadratic

    !
    !> Compute step size using perturbation method or the conventional line search method
    !
    subroutine compute_step_size

        ! Define a relaxation factor that allows later-iteration misfit > current iteration
        if (trigger_jumpout) then
            ! When there are three same misfits, jump out by default
            call readpar_xfloat(file_parameter, 'jumpout_factor', jumpout_factor, 1.05, 1.0*iter)
        else
            call readpar_xfloat(file_parameter, 'jumpout_factor', jumpout_factor, 1.0, 1.0*iter)
        end if

        select case (step_size_method)
            case ('linear')
                call compute_step_size_linear
            case ('quadratic')
                call compute_step_size_quadratic
            case ('line-search')
                call compute_step_size_linesearch
        end select

        call mpibarrier

    end subroutine compute_step_size

end module inversion_step_size
