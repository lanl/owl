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


submodule(acoustic_iso_3d) acoustic_iso_3d_forward_modeling

    use libflit
    use acoustic_iso_3d_vars
    use acoustic_iso_3d_boundary_saving
    use acoustic_iso_3d_wavefield

    implicit none

contains

    module subroutine compute_forward(this)

        class(wave_solver_acoustic_iso_3d), intent(inout) :: this

        integer :: i, j, k, l, t, ir, irx, iry, irz, rgx, rgy, rgz, kr
        real :: wmin, wmax, sr
        logical :: wnan

        ! Prepare modeling
        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield

        ! Seismic data
        call this%seis_p%init(nt=nt, dt=dt, nr=sgmtr%nr)

        ! Prepare boundary saving
        if (yn_reconstruct) then
            call prepare_boundary_saving
            call open_boundary_saving
        end if

        call alloc_array(snapp, [1, nx, 1, ny, 1, nz], pad=pml)

        ! Modeling
        l = 1
        do t = 1, nt

            ! Save boundary wavefield for FWI
            if (yn_reconstruct) then
                call save_boundary_wavefield(t)
                if (t == nt) then
                    call output_final_step_wavefield
                end if
            end if

            ! Compute seismic data; each block sums only the stencil points it holds, and
            ! collect_group sums the partial traces over the shot group. With a free surface
            ! (p = 0 on row 1), the interpolation points above it take the pressure mirrored
            ! about it, p(1 - k) = -p(1 + k), read from the rows below the surface (the rows
            ! above it are the halo of the top block)
            !$omp parallel do private(ir, irx, iry, irz, rgx, rgy, rgz, kr, sr)
            do ir = 1, sgmtr%nr
                rgx = sgmtr%recr(ir)%gx
                rgy = sgmtr%recr(ir)%gy
                rgz = sgmtr%recr(ir)%gz
                if (sgmtr%recr(ir)%weight /= 0) then
                    do irz = -nkw, nkw
                        kr = rgz + irz
                        sr = 1.0
                        if (yn_free_surface .and. kr < 1) then
                            kr = 2 - kr
                            sr = -1.0
                        end if
                        do iry = -nkw, nkw
                            do irx = -nkw, nkw
                                if (is_in_block(rgx + irx, rgy + iry, kr)) then
                                    this%seis_p%trace(ir)%data(t) = this%seis_p%trace(ir)%data(t) &
                                        + sr*p(rgx + irx, rgy + iry, kr) &
                                        *sgmtr%recr(ir)%interp_ix(irx) &
                                        *sgmtr%recr(ir)%interp_iy(iry) &
                                        *sgmtr%recr(ir)%interp_iz(irz) &
                                        *sgmtr%recr(ir)%weight
                                end if
                            end do
                        end do
                    end do
                end if
            end do
            !$omp end parallel do

            ! Record wavefield snapshot if necessary
            if (np /= 0 .and. l <= np) then
                if (t - 1 == nint(snaps(l)/dt)) then

                    snapp = 0.0

                    !$omp parallel do private(i, j, k) collapse(3)
                    do k = nz1, nz2
                        do j = ny1, ny2
                            do i = nx1, nx2
                                snapp(i, j, k) = p(i, j, k)
                            end do
                        end do
                    end do
                    !$omp end parallel do

                    call reduce_array_group(snapp)

                    if (rankid_group == 0) then
                        call output_array(snapp(1:nx, 1:ny, 1:nz), tidy(dir_snapshot)//'/shot_' &
                            //num2str(sgmtr%id) &
                            //'_forward_wavefield_p_' &
                            //num2str(l)//'.bin', store=321)
                    end if

                    l = l + 1

                end if
            end if

            ! Update wavefield
            if (yn_free_surface) then
                call update_wavefield_free_surface(dt, &
                    p, vx, vy, vz, &
                    memory_pdxp_xmin, memory_pdxp_xmax, &
                    memory_pdyp_ymin, memory_pdyp_ymax, &
                    memory_pdzp_zmax, &
                    memory_pdxvx_xmin, memory_pdxvx_xmax, &
                    memory_pdyvy_ymin, memory_pdyvy_ymax, &
                    memory_pdzvz_zmax)
            else
                call update_wavefield(dt, &
                    p, vx, vy, vz, &
                    memory_pdxp_xmin, memory_pdxp_xmax, &
                    memory_pdyp_ymin, memory_pdyp_ymax, &
                    memory_pdzp_zmin, memory_pdzp_zmax, &
                    memory_pdxvx_xmin, memory_pdxvx_xmax, &
                    memory_pdyvy_ymin, memory_pdyvy_ymax, &
                    memory_pdzvz_zmin, memory_pdzvz_zmax)
            end if

            ! Source term
            call add_source(t)

            if (verbose .and. (mod(t, max(nint(nt/10.0), 1)) == 0 .or. t == 1 .or. t == nt)) then

                wnan = group_and(any(isnan(p)))

                wmin = group_min(minval(p))
                wmax = group_max(maxval(p))

                if (rankid_group == 0) then
                    call warn(date_time_compact()//' >> Shot '//num2str(sgmtr%id) &
                        //' forward modeling step '//num2str(t)//' of '//num2str(nt))
                    if (wnan) then
                        call warn(date_time_compact()//' >> P contain NaN!')
                        stop
                    else
                        call warn(date_time_compact()//' >> P value range = ')
                        call warn(date_time_compact()//'      '//num2str(wmin, '(es)')//' ~ '//num2str(wmax, '(es)'))
                    end if
                end if

            end if

        end do

        ! Close boundary saving files
        if (yn_reconstruct) then
            call close_boundary_saving
        end if

        ! Resample data if necessary
        call this%seis_p%resamp(nnt=data_nt, ddt=data_dt)
        call this%seis_p%collect_group

        ! Output
        if (rankid_group == 0) then
            call this%seis_p%output(tidy(dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_seismogram_p.su')
        end if

        call mpibarrier_group

    end subroutine

end submodule
