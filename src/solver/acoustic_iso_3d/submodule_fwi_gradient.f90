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


submodule(acoustic_iso_3d) acoustic_iso_3d_fwi_gradient

    use libflit
    use acoustic_iso_3d_vars
    use acoustic_iso_3d_boundary_saving
    use acoustic_iso_3d_wavefield

    implicit none

#define interior_region nx1_interior:nx2_interior, ny1_interior:ny2_interior, nz1_interior:nz2_interior

contains

    module subroutine compute_adjoint(this)

        class(wave_solver_acoustic_iso_3d), intent(inout) :: this

        integer :: l, i, j, k, ir, irx, iry, irz, rgx, rgy, rgz, t
        type(grid3) :: grd
        logical :: wnan
        real :: wmin, wmax

        yn_energy_precond = this%energy_precond

        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield
        call alloc_adjoint_wavefield

        call this%seis_pr%load(tidy(dir_adjoint)//'/shot_'//tidy(num2str(sgmtr%id))//'_seismogram_p.su')
        call this%seis_pr%zero_foreign_rank_traces_group
        call this%seis_pr%resamp(nnt=nt, ddt=dt)
        call this%seis_pr%collect_group

        call alloc_array(prev_p, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_vx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_vy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_vz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(src_p, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(rec_p, [nx1, nx2, ny1, ny2, nz1, nz2])

        call alloc_array(src_vx, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(rec_vx, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(src_vy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(rec_vy, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(src_vz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(rec_vz, [nx1, nx2, ny1, ny2, nz1, nz2])

        grad_vp = zeros(nx, ny, nz)
        grad_rho = zeros(nx, ny, nz)

        energy_src_v = zeros(nx, ny, nz)
        energy_rec_v = zeros(nx, ny, nz)
        energy_src_a = zeros(nx, ny, nz)
        energy_rec_a = zeros(nx, ny, nz)

        call alloc_array(snapp, [1, nx, 1, ny, 1, nz], pad=pml)

        call prepare_boundary_saving
        call open_boundary_saving

        l = np
        do t = nt, sgmtr%srcr(1)%hnt, -1

            prev_p = p
            prev_vx = vx
            prev_vy = vy
            prev_vz = vz

            ! ===========================================================================
            !  Source wavefield reconstruction

            ! Read final step wavefield
            if (t == nt) then
                call input_final_step_wavefield
            end if

            ! Inject boundary wavefield
            call inject_boundary_wavefield(t)

            if (yn_free_surface) then
                call update_wavefield_free_surface(-dt, &
                    p, vx, vy, vz, &
                    memory_pdxp_xmin, memory_pdxp_xmax, &
                    memory_pdyp_ymin, memory_pdyp_ymax, &
                    memory_pdzp_zmax, &
                    memory_pdxvx_xmin, memory_pdxvx_xmax, &
                    memory_pdyvy_ymin, memory_pdyvy_ymax, &
                    memory_pdzvz_zmax)
            else
                call update_wavefield(-dt, &
                    p, vx, vy, vz, &
                    memory_pdxp_xmin, memory_pdxp_xmax, &
                    memory_pdyp_ymin, memory_pdyp_ymax, &
                    memory_pdzp_zmin, memory_pdzp_zmax, &
                    memory_pdxvx_xmin, memory_pdxvx_xmax, &
                    memory_pdyvy_ymin, memory_pdyvy_ymax, &
                    memory_pdzvz_zmin, memory_pdzvz_zmax)
            end if

            ! ---------------------- Adjoint wavefield reverse-time propagation ----------------------
            if (yn_free_surface) then
                call update_wavefield_free_surface(-dt, &
                    pr, vxr, vyr, vzr, &
                    memory_pdxpr_xmin, memory_pdxpr_xmax, &
                    memory_pdypr_ymin, memory_pdypr_ymax, &
                    memory_pdzpr_zmax, &
                    memory_pdxvxr_xmin, memory_pdxvxr_xmax, &
                    memory_pdyvyr_ymin, memory_pdyvyr_ymax, &
                    memory_pdzvzr_zmax)
            else
                call update_wavefield(-dt, &
                    pr, vxr, vyr, vzr, &
                    memory_pdxpr_xmin, memory_pdxpr_xmax, &
                    memory_pdypr_ymin, memory_pdypr_ymax, &
                    memory_pdzpr_zmin, memory_pdzpr_zmax, &
                    memory_pdxvxr_xmin, memory_pdxvxr_xmax, &
                    memory_pdyvyr_ymin, memory_pdyvyr_ymax, &
                    memory_pdzvzr_zmin, memory_pdzvzr_zmax)
            end if

            ! Save reconstructed source wavefield if necessary
            if (np /= 0 .and. l >= 1) then
                if (t - 1 == nint(snaps(l)/dt)) then

                    snapp = 0

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
                        call output_array(snapp(1:nx, 1:ny, 1:nz), tidy(dir_working)//'/shot_' &
                            //num2str(sgmtr%id) &
                            //'_reconstructed_wavefield_p_' &
                            //num2str(l)//'.bin', store=321)
                    end if

                    l = l - 1
                end if
            end if

            ! Read and add adjoint source
            !$omp parallel private(ir, irx, iry, irz, rgx, rgy, rgz)
            do ir = 1, sgmtr%nr
                if (sgmtr%recr(ir)%weight /= 0) then

                    rgx = sgmtr%recr(ir)%gx
                    rgy = sgmtr%recr(ir)%gy
                    rgz = sgmtr%recr(ir)%gz

                    if (is_in_block(rgx, rgy, rgz)) then

                        !$omp do collapse(3)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    pr(rgx + irx, rgy + iry, rgz + irz) = &
                                        pr(rgx + irx, rgy + iry, rgz + irz) &
                                        + this%seis_pr%trace(ir)%data(t) &
                                        *sgmtr%recr(ir)%interp_ix(irx) &
                                        *sgmtr%recr(ir)%interp_iy(iry) &
                                        *sgmtr%recr(ir)%interp_iz(irz) &
                                        *sgmtr%recr(ir)%weight
                                end do
                            end do
                        end do
                        !$omp end do

                    end if

                end if
            end do
            !$omp end parallel

            ! Compute gradients
            if (mod(t, cc_step_interval) == 0) then
                call compute_gradient
            end if

            ! Print progress
            if (verbose .and. (mod(t, max(nint(nt/10.0), 1)) == 0 .or. t == 1 .or. t == nt)) then

                wnan = group_and(any(isnan(pr)))

                wmin = group_min(minval(pr))
                wmax = group_max(maxval(pr))

                if (rankid_group == 0) then
                    call warn(date_time_compact()//' >> Shot '//num2str(sgmtr%id) &
                        //' FWI gradient computation step '//num2str(t)//' of '//num2str(nt))
                    if (wnan) then
                        call warn(date_time_compact()//' >> Pr contains NaN!')
                        stop
                    else
                        call warn(date_time_compact()//' >> Pr value range = ')
                        call warn(date_time_compact()//'      '//num2str(wmin, '(es)')//' ~ '//num2str(wmax, '(es)'))
                    end if
                end if

            end if

        end do

        ! Delete boundary saving files
        call close_boundary_saving(delete=.true.)

        ! Process and output gradients
        if (yn_energy_precond) then

            call allreduce_array_group(grad_vp)
            call allreduce_array_group(grad_rho)
            call allreduce_array_group(energy_src_v)
            call allreduce_array_group(energy_rec_v)
            call allreduce_array_group(energy_src_a)
            call allreduce_array_group(energy_rec_a)

            energy_src_v = energy_src_v + 1.0e-3*maxval(energy_src_v)
            energy_rec_v = energy_rec_v + 1.0e-3*maxval(energy_rec_v)
            energy_src_v = sqrt(energy_src_v*energy_rec_v)
            grad_vp = grad_vp/energy_src_v

            energy_src_a = energy_src_a + 1.0e-3*maxval(energy_src_a)
            energy_rec_a = energy_rec_a + 1.0e-3*maxval(energy_rec_a)
            energy_src_a = sqrt(energy_src_a*energy_rec_a)
            grad_rho = grad_rho/energy_src_a

        end if

        if (rankid_group == 0) then

            call grd%init(n=[nz, ny, nx], d=[dz, dy, dx], o=[oz, oy, ox])

            grd%array = permute(return_normal(grad_vp), 321)
            grd%array = return_normal(grd%array)
            call grd%output(tidy(dir_working)//'/shot_'//tidy(num2str(sgmtr%id, '(i)'))//'_grad_vp.grd')

            grd%array = permute(return_normal(grad_rho), 321)
            grd%array = return_normal(grd%array)
            call grd%output(tidy(dir_working)//'/shot_'//tidy(num2str(sgmtr%id, '(i)'))//'_grad_rho.grd')

        end if

        call mpibarrier_group

    end subroutine

    subroutine compute_gradient

        integer :: i, j, k

        ! Vp
        !$omp parallel do private(i, j, k) collapse(3)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1, nx2
                    src_p(i, j, k) = p(i, j, k) - prev_p(i, j, k)
                    rec_p(i, j, k) = pr(i, j, k)
                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k) collapse(3)
        do k = nz1_interior, nz2_interior
            do j = ny1_interior, ny2_interior
                do i = nx1_interior, nx2_interior
                    grad_vp(i, j, k) = grad_vp(i, j, k) - src_p(i, j, k)*rec_p(i, j, k)/(rho(i, j, k)*vp(i, j, k)**3)
                end do
            end do
        end do
        !$omp end parallel do

        if (yn_energy_precond) then

            !$omp parallel do private(i, j, k) collapse(3)
            do k = nz1_interior, nz2_interior
                do j = ny1_interior, ny2_interior
                    do i = nx1_interior, nx2_interior
                        energy_src_v(i, j, k) = energy_src_v(i, j, k) + src_p(i, j, k)**2/(rho(i, j, k)*vp(i, j, k)**3)
                        energy_rec_v(i, j, k) = energy_rec_v(i, j, k) + rec_p(i, j, k)**2/(rho(i, j, k)*vp(i, j, k)**3)
                    end do
                end do
            end do
            !$omp end parallel do

        end if

        ! Density
        !$omp parallel do private(i, j, k) collapse(3)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1, nx2
                    src_vx(i, j, k) = sum(vx(i:i + 1, j, k)) - sum(prev_vx(i:i + 1, j, k))
                    rec_vx(i, j, k) = sum(vxr(i:i + 1, j, k))
                    src_vy(i, j, k) = sum(vy(i, j:j + 1, k)) - sum(prev_vy(i, j:j + 1, k))
                    rec_vy(i, j, k) = sum(vyr(i, j:j + 1, k))
                    src_vz(i, j, k) = sum(vz(i, j, k:k + 1)) - sum(prev_vz(i, j, k:k + 1))
                    rec_vz(i, j, k) = sum(vzr(i, j, k:k + 1))
                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k) collapse(3)
        do k = nz1_interior, nz2_interior
            do j = ny1_interior, ny2_interior
                do i = nx1_interior, nx2_interior
                    grad_rho(i, j, k) = grad_rho(i, j, k) + src_vx(i, j, k)*rec_vx(i, j, k) &
                        + src_vy(i, j, k)*rec_vy(i, j, k) + src_vz(i, j, k)*rec_vz(i, j, k)
                end do
            end do
        end do
        !$omp end parallel do

        if (yn_energy_precond) then

            !$omp parallel do private(i, j, k) collapse(3)
            do k = nz1_interior, nz2_interior
                do j = ny1_interior, ny2_interior
                    do i = nx1_interior, nx2_interior
                        energy_src_a(i, j, k) = energy_src_a(i, j, k) &
                            + src_vx(i, j, k)**2 + src_vy(i, j, k)**2 + src_vz(i, j, k)**2
                        energy_rec_a(i, j, k) = energy_rec_a(i, j, k) &
                            + rec_vx(i, j, k)**2 + src_vy(i, j, k)**2 + rec_vz(i, j, k)**2
                    end do
                end do
            end do
            !$omp end parallel do

        end if

    end subroutine

end submodule
