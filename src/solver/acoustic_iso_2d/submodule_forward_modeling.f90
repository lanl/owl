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


submodule(acoustic_iso_2d) acoustic_iso_2d_forward_modeling

    use libflit
    use acoustic_iso_2d_vars
    use acoustic_iso_2d_boundary_saving
    use acoustic_iso_2d_wavefield

    implicit none

contains

    module subroutine compute_forward(this)

        class(wave_solver_acoustic_iso_2d), intent(inout) :: this

        integer :: l, ir, irx, irz, rgx, rgz, kr, t
        real :: sr

        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield

        call this%seis_p%init(nt=nt, dt=dt, nr=sgmtr%nr)

        if (yn_reconstruct) then
            call prepare_boundary_saving
            call open_boundary_saving
        end if

        l = 1
        do t = 1, nt

            if (yn_reconstruct) then
                call save_boundary_wavefield(t)
                if (t == nt) then
                    call output_final_step_wavefield
                end if
            end if

            ! Record seismogram. With a free surface (p = 0 on row 1), the interpolation points
            ! above it take the pressure mirrored about it, p(1 - k) = -p(1 + k), read from the
            ! rows below the surface (the rows above it are refreshed only once per time step)
            !$omp parallel do private(ir, irx, irz, rgx, rgz, kr, sr)
            do ir = 1, sgmtr%nr

                if (sgmtr%recr(ir)%weight /= 0) then

                    rgx = sgmtr%recr(ir)%gx
                    rgz = sgmtr%recr(ir)%gz

                    do irz = -nkw, nkw
                        kr = rgz + irz
                        sr = 1.0
                        if (yn_free_surface .and. kr < 1) then
                            kr = 2 - kr
                            sr = -1.0
                        end if
                        do irx = -nkw, nkw
                            this%seis_p%trace(ir)%data(t) = &
                                this%seis_p%trace(ir)%data(t) + &
                                sr*p(rgx + irx, kr) &
                                *sgmtr%recr(ir)%interp_ix(irx) &
                                *sgmtr%recr(ir)%interp_iz(irz) &
                                *sgmtr%recr(ir)%weight
                        end do
                    end do

                end if
            end do
            !$omp end parallel do

            ! Record wavefield snapshot if necessary
            if (np /= 0 .and. l <= np) then
                if (t - 1 == nint(snaps(l)/dt)) then

                    call output_array(p(1:nx, 1:nz), tidy(dir_snapshot)//'/shot_' &
                        //num2str(sgmtr%id) &
                        //'_forward_wavefield_p_' &
                        //num2str(l)//'.bin', transp=.true.)

                    l = l + 1
                end if
            end if

            ! Update wavefields
            if (yn_free_surface) then
                call update_wavefield_free_surface(dt, p, vx, vz, memory_pdxvx, memory_pdzvz, memory_pdxp, memory_pdzp)
            else
                call update_wavefield(dt, p, vx, vz, memory_pdxvx, memory_pdzvz, memory_pdxp, memory_pdzp)
            end if

            ! Source term
            call add_source(t)

            if (verbose .and. (mod(t, max(nint(nt/10.0), 1)) == 0 .or. t == 1 .or. t == nt)) then
                call warn(date_time_compact()//' >> Shot '//num2str(sgmtr%id) &
                    //' forward modeling step '//num2str(t)//' of '//num2str(nt))
                if (any(isnan(p))) then
                    call warn(date_time_compact()//' >> P contains NaN!')
                    stop
                else
                    call warn(date_time_compact()//' >> P value range = ')
                    call warn(date_time_compact()//'      '//num2str(minval(p), '(es)') &
                        //' ~ '//num2str(maxval(p), '(es)'))
                end if
            end if

        end do

        ! close boundary saving
        if (yn_reconstruct) then
            call close_boundary_saving
        end if

        ! resample if necessary
        call this%seis_p%resamp(nnt=data_nt, ddt=data_dt)

        ! output
        call this%seis_p%output(tidy(dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_seismogram_p.su')

        call mpibarrier_group

    end subroutine

end submodule
