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


submodule(acoustic_iso_2d) acoustic_iso_2d_fwi_gradient

    use libflit
    use acoustic_iso_2d_vars
    use acoustic_iso_2d_boundary_saving
    use acoustic_iso_2d_wavefield

    implicit none

contains

    module subroutine compute_adjoint(this)

        class(wave_solver_acoustic_iso_2d), intent(inout) :: this

        integer :: l, ir, irx, irz, rgx, rgz, t
        type(grid2) :: grd

        yn_energy_precond = this%energy_precond

        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield
        call alloc_adjoint_wavefield

        call this%seis_pr%load(tidy(dir_adjoint)//'/shot_'//tidy(num2str(sgmtr%id))//'_seismogram_p.su')
        call this%seis_pr%resamp(nnt=nt, ddt=dt)

        call alloc_array(prev_p, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_vx, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(prev_vz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(p_lrsh, [1, nx], pad=pml)
        call alloc_array(p_lrrh, [1, nx], pad=pml)
        call alloc_array(p_udsh, [1, nz], pad=pml)
        call alloc_array(p_udrh, [1, nz], pad=pml)

        call alloc_array(src_p, [1, nx, 1, nz], pad=pml)
        call alloc_array(rec_p, [1, nx, 1, nz], pad=pml)
        call alloc_array(src_vx, [1, nx, 1, nz], pad=pml)
        call alloc_array(rec_vx, [1, nx, 1, nz], pad=pml)
        call alloc_array(src_vz, [1, nx, 1, nz], pad=pml)
        call alloc_array(rec_vz, [1, nx, 1, nz], pad=pml)

        grad_vp = zeros(nx, nz)
        grad_rho = zeros(nx, nz)

        energy_src_v = zeros(nx, nz)
        energy_rec_v = zeros(nx, nz)
        energy_src_a = zeros(nx, nz)
        energy_rec_a = zeros(nx, nz)

        call prepare_boundary_saving
        call open_boundary_saving

        l = np
        do t = nt, sgmtr%srcr(1)%hnt, -1

            ! Save previous step wavefield
            prev_p = p
            prev_vx = vx
            prev_vz = vz

            ! ===========================================================================
            !  Source wavefield reconstruction

            ! Read final-step wavefield
            if (t == nt) then
                call input_final_step_wavefield
            end if

            ! Inject boundary wavefields
            call inject_boundary_wavefield(t)

            ! Reconstruct
            if (yn_free_surface) then
                call update_wavefield_free_surface(-dt, p, vx, vz, memory_pdxvx, memory_pdzvz, memory_pdxp, memory_pdzp)
            else
                call update_wavefield(-dt, p, vx, vz, memory_pdxvx, memory_pdzvz, memory_pdxp, memory_pdzp)
            end if

            ! ===========================================================================
            ! Adjoint wavefield reverse-time propagation

            if (yn_free_surface) then
                call update_wavefield_free_surface(-dt, pr, vxr, vzr, memory_pdxvxr, memory_pdzvzr, memory_pdxpr, memory_pdzpr)
            else
                call update_wavefield(-dt, pr, vxr, vzr, memory_pdxvxr, memory_pdzvzr, memory_pdxpr, memory_pdzpr)
            end if

            ! Record wavefield snapshot if necessary
            if (np /= 0 .and. l >= 1) then
                if (t - 1 == nint(snaps(l)/dt)) then

                    call output_array(p(1:nx, 1:nz), tidy(dir_working)//'/shot_' &
                        //num2str(sgmtr%id) &
                        //'_reconstructed_wavefield_p_' &
                        //num2str(l)//'.bin', transp=.true.)

                    l = l - 1
                end if
            end if

            ! Read and add adjoint source -- the receivers must be independent as when
            ! receivers are not at grid points, adjoint source may add to a same grid point
            ! at different/adjacent receivers.
            ! The following openmp does not parallize receivers, but only parallize within each receiver.
            !$omp parallel private(ir, irx, irz, rgx, rgz)
            do ir = 1, sgmtr%nr
                if (sgmtr%recr(ir)%weight /= 0) then

                    rgx = sgmtr%recr(ir)%gx
                    rgz = sgmtr%recr(ir)%gz

                    !$omp do collapse(2)
                    do irz = -nkw, nkw
                        do irx = -nkw, nkw
                            pr(rgx + irx, rgz + irz) = pr(rgx + irx, rgz + irz) &
                                + this%seis_pr%trace(ir)%data(t) &
                                *sgmtr%recr(ir)%interp_ix(irx) &
                                *sgmtr%recr(ir)%interp_iz(irz) &
                                *sgmtr%recr(ir)%weight
                        end do
                    end do
                    !$omp end do

                end if
            end do
            !$omp end parallel

            ! Gradient
            if (mod(t, cc_step_interval) == 0) then
                call compute_gradient
            end if

            if (verbose .and. (mod(t, max(nint(nt/10.0), 1)) == 0 .or. t == 1 .or. t == nt)) then
                call warn(date_time_compact()//' >> Shot '//num2str(sgmtr%id) &
                    //' FWI gradient computation step '//num2str(t)//' of '//num2str(nt))
                if (any(isnan(pr))) then
                    call warn(date_time_compact()//' >> Pr contains NaN!')
                    stop
                else
                    call warn(date_time_compact()//' >> Pr value range = ')
                    call warn(date_time_compact()//'      '//num2str(minval(pr), '(es)') &
                        //' ~ '//num2str(maxval(pr), '(es)'))
                end if
            end if

        end do

        ! Delete boundary saving files
        call close_boundary_saving(delete=.true.)

        ! Process and output gradients
        if (yn_energy_precond) then

            energy_src_v = energy_src_v + 1.0e-3*maxval(energy_src_v)
            energy_rec_v = energy_rec_v + 1.0e-3*maxval(energy_rec_v)
            energy_src_v = sqrt(energy_src_v*energy_rec_v)
            grad_vp = grad_vp/energy_src_v

            energy_src_a = energy_src_a + 1.0e-3*maxval(energy_src_a)
            energy_rec_a = energy_rec_a + 1.0e-3*maxval(energy_rec_a)
            energy_src_a = sqrt(energy_src_a*energy_rec_a)
            grad_rho = grad_rho/energy_src_a

        end if

        call grd%init(n=[nz, nx], d=[dz, dx], o=[oz, ox])

        grd%array = transpose(return_normal(grad_vp))
        grd%array = return_normal(grd%array)
        call grd%output(tidy(dir_working)//'/shot_'//tidy(num2str(sgmtr%id, '(i)'))//'_grad_vp.grd')

        grd%array = transpose(return_normal(grad_rho))
        grd%array = return_normal(grd%array)
        call grd%output(tidy(dir_working)//'/shot_'//tidy(num2str(sgmtr%id, '(i)'))//'_grad_rho.grd')

        call mpibarrier_group

    end subroutine

    subroutine compute_gradient

        integer :: i, j

        ! Vp
        !$omp parallel do private(i, j) collapse(2)
        do j = -pml + 1, nz + pml
            do i = -pml + 1, nx + pml
                src_p(i, j) = p(i, j) - prev_p(i, j)
                rec_p(i, j) = pr(i, j)
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j) collapse(2)
        do j = 1, nz
            do i = 1, nx
                grad_vp(i, j) = grad_vp(i, j) - src_p(i, j)*rec_p(i, j)/(rho(i, j)*vp(i, j)**3)
            end do
        end do
        !$omp end parallel do

        if (yn_energy_precond) then

            !$omp parallel do private(i, j) collapse(2)
            do j = 1, nz
                do i = 1, nx
                    energy_src_v(i, j) = energy_src_v(i, j) + src_p(i, j)**2/(rho(i, j)*vp(i, j)**3)
                    energy_rec_v(i, j) = energy_rec_v(i, j) + rec_p(i, j)**2/(rho(i, j)*vp(i, j)**3)
                end do
            end do
            !$omp end parallel do

        end if

        ! Density
        !$omp parallel do private(i, j) collapse(2)
        do j = -pml + 1, nz + pml
            do i = -pml + 1, nx + pml
                src_vx(i, j) = sum(vx(i:i + 1, j)) - sum(prev_vx(i:i + 1, j))
                rec_vx(i, j) = sum(vxr(i:i + 1, j))
                src_vz(i, j) = sum(vz(i, j:j + 1)) - sum(prev_vz(i, j:j + 1))
                rec_vz(i, j) = sum(vzr(i, j:j + 1))
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j) collapse(2)
        do j = 1, nz
            do i = 1, nx
                grad_rho(i, j) = grad_rho(i, j) + src_vx(i, j)*rec_vx(i, j) + src_vz(i, j)*rec_vz(i, j)
            end do
        end do
        !$omp end parallel do

        if (yn_energy_precond) then

            !$omp parallel do private(i, j) collapse(2)
            do j = 1, nz
                do i = 1, nx
                    energy_src_a(i, j) = energy_src_a(i, j) + src_vx(i, j)**2 + src_vz(i, j)**2
                    energy_rec_a(i, j) = energy_rec_a(i, j) + rec_vx(i, j)**2 + rec_vz(i, j)**2
                end do
            end do
            !$omp end parallel do

        end if

    end subroutine

end submodule
