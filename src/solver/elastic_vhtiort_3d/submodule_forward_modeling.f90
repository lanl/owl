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


submodule(elastic_vhtiort_3d) elastic_vhtiort_3d_forward_modeling

    use libflit
    use elastic_vhtiort_3d_vars
    use elastic_vhtiort_3d_boundary_saving
    use elastic_vhtiort_3d_wavefield
    use elastic_vhtiort_3d_cfspml

    implicit none

contains

    module subroutine compute_forward(this)

        class(wave_solver_elastic_vhtiort_3d), intent(inout) :: this

        integer :: l, ir, irx, iry, irz, rgx, rgy, rgz, kr
        integer :: i, j, k, t
        real :: wmin1, wmin2, wmin3
        real :: wmax1, wmax2, wmax3
        logical :: wnan
        real, allocatable, dimension(:, :) :: seis_vx, seis_vy, seis_vz

        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield

        seis_vx = zeros(nt, sgmtr%nr)
        seis_vy = zeros(nt, sgmtr%nr)
        seis_vz = zeros(nt, sgmtr%nr)

        if (yn_reconstruct) then
            call prepare_boundary_saving
            call open_boundary_saving
        end if

        l = 1
        do t = 1, nt

            ! Save boundary wavefield for reconstruction
            if (yn_reconstruct) then
                call save_boundary_wavefield(t)
                if (t == nt) then
                    call output_final_step_wavefield
                end if
            end if

            ! Record seismogram. With a free surface, the interpolation points above it take the
            ! velocities mirrored about it, vx(1 - k) = vx(1 + k), vy likewise, vz(2 - k) = vz(1 + k),
            ! read from the rows below the surface (those above it are the halo of the top block)
            !$omp parallel do private(ir, irx, iry, irz, rgx, rgy, rgz, kr) schedule(auto)
            do ir = 1, sgmtr%nr
                if (sgmtr%recr(ir)%weight /= 0) then

                    if (yn_compx) then
                        rgx = sgmtr%recr(ir)%hx
                        rgy = sgmtr%recr(ir)%gy
                        rgz = sgmtr%recr(ir)%gz
                        do irz = -nkw, nkw
                            kr = rgz + irz
                            if (yn_free_surface .and. kr < 1) then
                                kr = 2 - kr
                            end if
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(rgx + irx, rgy + iry, kr)) then
                                        seis_vx(t, ir) = seis_vx(t, ir) + vx(rgx + irx, rgy + iry, kr) &
                                            *sgmtr%recr(ir)%interp_hx(irx) &
                                            *sgmtr%recr(ir)%interp_iy(iry) &
                                            *sgmtr%recr(ir)%interp_iz(irz) &
                                            *sgmtr%recr(ir)%weight
                                    end if
                                end do
                            end do
                        end do
                    end if

                    if (yn_compy) then
                        rgx = sgmtr%recr(ir)%gx
                        rgy = sgmtr%recr(ir)%hy
                        rgz = sgmtr%recr(ir)%gz
                        do irz = -nkw, nkw
                            kr = rgz + irz
                            if (yn_free_surface .and. kr < 1) then
                                kr = 2 - kr
                            end if
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(rgx + irx, rgy + iry, kr)) then
                                        seis_vy(t, ir) = seis_vy(t, ir) + vy(rgx + irx, rgy + iry, kr) &
                                            *sgmtr%recr(ir)%interp_ix(irx) &
                                            *sgmtr%recr(ir)%interp_hy(iry) &
                                            *sgmtr%recr(ir)%interp_iz(irz) &
                                            *sgmtr%recr(ir)%weight
                                    end if
                                end do
                            end do
                        end do
                    end if

                    if (yn_compz) then
                        rgx = sgmtr%recr(ir)%gx
                        rgy = sgmtr%recr(ir)%gy
                        rgz = sgmtr%recr(ir)%hz
                        do irz = -nkw, nkw
                            kr = rgz + irz
                            if (yn_free_surface .and. kr < 2) then
                                kr = 3 - kr
                            end if
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(rgx + irx, rgy + iry, kr)) then
                                        seis_vz(t, ir) = seis_vz(t, ir) + vz(rgx + irx, rgy + iry, kr) &
                                            *sgmtr%recr(ir)%interp_ix(irx) &
                                            *sgmtr%recr(ir)%interp_iy(iry) &
                                            *sgmtr%recr(ir)%interp_hz(irz) &
                                            *sgmtr%recr(ir)%weight
                                    end if
                                end do
                            end do
                        end do
                    end if

                end if
            end do
            !$omp end parallel do

            ! Record wavefield snapshot if necessary
            if (np /= 0 .and. l <= np) then
                if (t - 1 == nint(snaps(l)/dt)) then

                    call alloc_array(snapvx, [1, nx, 1, ny, 1, nz], pad=pml)
                    call alloc_array(snapvy, [1, nx, 1, ny, 1, nz], pad=pml)
                    call alloc_array(snapvz, [1, nx, 1, ny, 1, nz], pad=pml)

                    ! The velocities changed after their halo exchange; refresh the
                    ! halos that the averages below read across block faces
                    call commute_array_group(vx, fdhalf)
                    call commute_array_group(vy, fdhalf)
                    call commute_array_group(vz, fdhalf)

                    if (yn_free_surface) then

                        !$omp parallel do private(i, j, k) collapse(3)
                        do k = nz1, nz2
                            do j = ny1, ny2
                                do i = nx1, nx2
                                    snapvx(i, j, k) = 0.5*sum(vx(i:i + 1, j, k))
                                    snapvy(i, j, k) = 0.5*sum(vy(i, j:j + 1, k))
                                    snapvz(i, j, k) = 0.5*sum(vz(i, j, k:k + 1))
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        call reduce_array_group(snapvx)
                        call reduce_array_group(snapvy)
                        call reduce_array_group(snapvz)

                        if (rankid_group == 0) then

                            call map_irregular_to_regular(snapvx, this, [1, this%nx, 1, this%ny, 1, this%nz])
                            call map_irregular_to_regular(snapvy, this, [1, this%nx, 1, this%ny, 1, this%nz])
                            call map_irregular_to_regular(snapvz, this, [1, this%nx, 1, this%ny, 1, this%nz])

                            call output_array(snapvx, tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_forward_wavefield_x_' &
                                //num2str(l)//'.bin', store=321)
                            call output_array(snapvy, tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_forward_wavefield_y_' &
                                //num2str(l)//'.bin', store=321)
                            call output_array(snapvz, tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_forward_wavefield_z_' &
                                //num2str(l)//'.bin', store=321)

                        end if

                    else

                        !$omp parallel do private(i, j, k) collapse(3)
                        do k = nz1, nz2
                            do j = ny1, ny2
                                do i = nx1, nx2
                                    snapvx(i, j, k) = 0.5*sum(vx(i:i + 1, j, k))
                                    snapvy(i, j, k) = 0.5*sum(vy(i, j:j + 1, k))
                                    snapvz(i, j, k) = 0.5*sum(vz(i, j, k:k + 1))
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        call reduce_array_group(snapvx)
                        call reduce_array_group(snapvy)
                        call reduce_array_group(snapvz)

                        ! Output
                        if (rankid_group == 0) then
                            call output_array(snapvx(1:nx, 1:ny, 1:nz), tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_forward_wavefield_x_' &
                                //num2str(l)//'.bin', store=321)
                            call output_array(snapvy(1:nx, 1:ny, 1:nz), tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_forward_wavefield_y_' &
                                //num2str(l)//'.bin', store=321)
                            call output_array(snapvz(1:nx, 1:ny, 1:nz), tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_forward_wavefield_z_' &
                                //num2str(l)//'.bin', store=321)
                        end if

                    end if

                    l = l + 1
                end if
            end if

            ! Wavefield update
            if (yn_free_surface) then
                call update_wavefield_free_surface(dt, &
                    stressxx, stressyy, stresszz, &
                    stressyz, stressxz, stressxy, &
                    vx, vy, vz, &
                    memory_pdxxx, memory_pdyxy, memory_pdzxz, &
                    memory_pdxxy, memory_pdyyy, memory_pdzyz, &
                    memory_pdxxz, memory_pdyyz, memory_pdzzz, &
                    memory_pdxvx, memory_pdyvx, memory_pdzvx, &
                    memory_pdxvy, memory_pdyvy, memory_pdzvy, &
                    memory_pdxvz, memory_pdyvz, memory_pdzvz)
            else
                call update_wavefield(dt, &
                    stressxx, stressyy, stresszz, &
                    stressyz, stressxz, stressxy, &
                    vx, vy, vz, &
                    memory_pdxxx, memory_pdyxy, memory_pdzxz, &
                    memory_pdxxy, memory_pdyyy, memory_pdzyz, &
                    memory_pdxxz, memory_pdyyz, memory_pdzzz, &
                    memory_pdxvx, memory_pdyvx, memory_pdzvx, &
                    memory_pdxvy, memory_pdyvy, memory_pdzvy, &
                    memory_pdxvz, memory_pdyvz, memory_pdzvz)
            end if

            ! Source term
            call add_source(t)

            if (verbose .and. (mod(t, max(nint(nt/10.0), 1)) == 0 .or. t == 1 .or. t == nt)) then

                wnan = group_and(any(isnan(vx)) .or. any(isnan(vy)) .or. any(isnan(vz)))

                wmin1 = group_min(minval(vx))
                wmax1 = group_max(maxval(vx))

                wmin2 = group_min(minval(vy))
                wmax2 = group_max(maxval(vy))

                wmin3 = group_min(minval(vz))
                wmax3 = group_max(maxval(vz))

                if (rankid_group == 0) then
                    call warn(date_time_compact()//' >> Shot '//num2str(sgmtr%id) &
                        //' forward modeling step '//num2str(t)//' of '//num2str(nt))
                    if (wnan) then
                        call warn(date_time_compact()//' >> Vx, Vy, Vz contain NaN!')
                        stop
                    else
                        call warn(date_time_compact()//' >> Vx, Vy, Vz value range = ')
                        call warn(date_time_compact()//'      '//num2str(wmin1, '(es)')//' ~ '//num2str(wmax1, '(es)'))
                        call warn(date_time_compact()//'      '//num2str(wmin2, '(es)')//' ~ '//num2str(wmax2, '(es)'))
                        call warn(date_time_compact()//'      '//num2str(wmin3, '(es)')//' ~ '//num2str(wmax3, '(es)'))
                    end if
                end if

            end if

        end do

        ! close boundary saving
        if (yn_reconstruct) then
            call close_boundary_saving
        end if

        ! resample if necessary
        if (yn_compx) then
            call allreduce_array_group(seis_vx)
            do i = 1, sgmtr%nr
                if (.not. (i >= trace_range(rankid_group, 1) .and. i <= trace_range(rankid_group, 2))) then
                    seis_vx(:, i) = 0.0
                end if
            end do
            call this%seis_vx%init(nt=nt, dt=dt, nr=sgmtr%nr)
            call this%seis_vx%from_array(seis_vx)
            call this%seis_vx%resamp(nnt=data_nt, ddt=data_dt)
            call this%seis_vx%collect_group
            if (rankid_group == 0) then
                call this%seis_vx%output(tidy(dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_seismogram_x.su')
            end if
        end if

        if (yn_compy) then
            call allreduce_array_group(seis_vy)
            do i = 1, sgmtr%nr
                if (.not. (i >= trace_range(rankid_group, 1) .and. i <= trace_range(rankid_group, 2))) then
                    seis_vy(:, i) = 0.0
                end if
            end do
            call this%seis_vy%init(nt=nt, dt=dt, nr=sgmtr%nr)
            call this%seis_vy%from_array(seis_vy)
            call this%seis_vy%resamp(nnt=data_nt, ddt=data_dt)
            call this%seis_vy%collect_group
            if (rankid_group == 0) then
                call this%seis_vy%output(tidy(dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_seismogram_y.su')
            end if
        end if

        if (yn_compz) then
            call allreduce_array_group(seis_vz)
            do i = 1, sgmtr%nr
                if (.not. (i >= trace_range(rankid_group, 1) .and. i <= trace_range(rankid_group, 2))) then
                    seis_vz(:, i) = 0.0
                end if
            end do
            call this%seis_vz%init(nt=nt, dt=dt, nr=sgmtr%nr)
            call this%seis_vz%from_array(seis_vz)
            call this%seis_vz%resamp(nnt=data_nt, ddt=data_dt)
            call this%seis_vz%collect_group
            if (rankid_group == 0) then
                call this%seis_vz%output(tidy(dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_seismogram_z.su')
            end if
        end if

        call mpibarrier_group

    end subroutine compute_forward

end submodule
