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


submodule(elastic_tti_3d) elastic_tti_3d_forward_modeling

    use libflit
    use elastic_tti_3d_vars
    use elastic_tti_3d_boundary_saving
    use elastic_tti_3d_wavefield
    use elastic_tti_3d_cfspml

    implicit none

contains

    module subroutine compute_forward(this)

        class(wave_solver_elastic_tti_3d), intent(inout) :: this

        integer :: l, ir, irx, iry, irz, rgx, rgy, rgz, kr
        integer :: i, j, k, t
        real :: wmin1, wmin2, wmin3
        real :: wmax1, wmax2, wmax3
        logical :: wnan
        real :: amp1, amp2, amp3, amp4
        real, allocatable, dimension(:, :) :: seis_vx, seis_vy, seis_vz
        real, allocatable, dimension(:) :: pt, rec_slopex, rec_slopey

        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield

        call alloc_array(snapvx, [1, nx, 1, ny, 1, nz], pad=pml)
        call alloc_array(snapvy, [1, nx, 1, ny, 1, nz], pad=pml)
        call alloc_array(snapvz, [1, nx, 1, ny, 1, nz], pad=pml)

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
            ! velocities mirrored about it, v(1 - k) = v(1 + k) on the integer-z sets and
            ! v(2 - k) = v(1 + k) on the half-z sets, read from the rows below the surface (those
            ! above it are the halo of the top block)
            !$omp parallel do private(ir, irx, iry, irz, rgx, rgy, rgz, kr, amp1, amp2, amp3, amp4) schedule(auto)
            do ir = 1, sgmtr%nr
                if (sgmtr%recr(ir)%weight /= 0) then

                    ! Set a
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

                                    amp1 = sgmtr%recr(ir)%interp_hx(irx) &
                                        *sgmtr%recr(ir)%interp_iy(iry) &
                                        *sgmtr%recr(ir)%interp_iz(irz) &
                                        *sgmtr%recr(ir)%weight*0.25

                                    seis_vx(t, ir) = &
                                        seis_vx(t, ir) &
                                        + vx_hxiyiz(rgx + irx, rgy + iry, kr)*amp1
                                    seis_vy(t, ir) = &
                                        seis_vy(t, ir) &
                                        + vy_hxiyiz(rgx + irx, rgy + iry, kr)*amp1
                                    seis_vz(t, ir) = &
                                        seis_vz(t, ir) &
                                        + vz_hxiyiz(rgx + irx, rgy + iry, kr)*amp1

                                end if

                            end do
                        end do
                    end do

                    ! Set b
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

                                    amp2 = sgmtr%recr(ir)%interp_ix(irx) &
                                        *sgmtr%recr(ir)%interp_hy(iry) &
                                        *sgmtr%recr(ir)%interp_iz(irz) &
                                        *sgmtr%recr(ir)%weight*0.25

                                    seis_vx(t, ir) = &
                                        seis_vx(t, ir) &
                                        + vx_ixhyiz(rgx + irx, rgy + iry, kr)*amp2
                                    seis_vy(t, ir) = &
                                        seis_vy(t, ir) &
                                        + vy_ixhyiz(rgx + irx, rgy + iry, kr)*amp2
                                    seis_vz(t, ir) = &
                                        seis_vz(t, ir) &
                                        + vz_ixhyiz(rgx + irx, rgy + iry, kr)*amp2
                                end if

                            end do
                        end do
                    end do

                    ! Set c
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

                                    amp3 = sgmtr%recr(ir)%interp_ix(irx) &
                                        *sgmtr%recr(ir)%interp_iy(iry) &
                                        *sgmtr%recr(ir)%interp_hz(irz) &
                                        *sgmtr%recr(ir)%weight*0.25

                                    seis_vx(t, ir) = &
                                        seis_vx(t, ir) &
                                        + vx_ixiyhz(rgx + irx, rgy + iry, kr)*amp3
                                    seis_vy(t, ir) = &
                                        seis_vy(t, ir) &
                                        + vy_ixiyhz(rgx + irx, rgy + iry, kr)*amp3
                                    seis_vz(t, ir) = &
                                        seis_vz(t, ir) &
                                        + vz_ixiyhz(rgx + irx, rgy + iry, kr)*amp3
                                end if

                            end do
                        end do
                    end do

                    ! Set d
                    rgx = sgmtr%recr(ir)%hx
                    rgy = sgmtr%recr(ir)%hy
                    rgz = sgmtr%recr(ir)%hz
                    do irz = -nkw, nkw
                        kr = rgz + irz
                        if (yn_free_surface .and. kr < 2) then
                            kr = 3 - kr
                        end if
                        do iry = -nkw, nkw
                            do irx = -nkw, nkw

                                if (is_in_block(rgx + irx, rgy + iry, kr)) then

                                    amp4 = sgmtr%recr(ir)%interp_hx(irx) &
                                        *sgmtr%recr(ir)%interp_hy(iry) &
                                        *sgmtr%recr(ir)%interp_hz(irz) &
                                        *sgmtr%recr(ir)%weight*0.25

                                    seis_vx(t, ir) = &
                                        seis_vx(t, ir) &
                                        + vx_hxhyhz(rgx + irx, rgy + iry, kr)*amp4
                                    seis_vy(t, ir) = &
                                        seis_vy(t, ir) &
                                        + vy_hxhyhz(rgx + irx, rgy + iry, kr)*amp4
                                    seis_vz(t, ir) = &
                                        seis_vz(t, ir) &
                                        + vz_hxhyhz(rgx + irx, rgy + iry, kr)*amp4

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

                    ! The velocities changed after their halo exchange; refresh the
                    ! halos that the averages below read across block faces
                    call commute_array_group(vx_hxiyiz, fdhalf)
                    call commute_array_group(vy_hxiyiz, fdhalf)
                    call commute_array_group(vz_hxiyiz, fdhalf)
                    call commute_array_group(vx_ixhyiz, fdhalf)
                    call commute_array_group(vy_ixhyiz, fdhalf)
                    call commute_array_group(vz_ixhyiz, fdhalf)
                    call commute_array_group(vx_ixiyhz, fdhalf)
                    call commute_array_group(vy_ixiyhz, fdhalf)
                    call commute_array_group(vz_ixiyhz, fdhalf)
                    call commute_array_group(vx_hxhyhz, fdhalf)
                    call commute_array_group(vy_hxhyhz, fdhalf)
                    call commute_array_group(vz_hxhyhz, fdhalf)

                    if (yn_free_surface) then

                        call alloc_array(snapvx, [1, nx, 1, ny, 1, nz], pad=pml)
                        call alloc_array(snapvy, [1, nx, 1, ny, 1, nz], pad=pml)
                        call alloc_array(snapvz, [1, nx, 1, ny, 1, nz], pad=pml)

                        !$omp parallel do private(i, j, k) collapse(3)
                        do k = nz1, nz2
                            do j = ny1, ny2
                                do i = nx1, nx2
                                    snapvx(i, j, k) = &
                                        0.5*sum(vx_hxiyiz(i:i + 1, j, k)) &
                                        + 0.5*sum(vx_ixhyiz(i, j:j + 1, k)) &
                                        + 0.5*sum(vx_ixiyhz(i, j, k:k + 1)) &
                                        + 0.125*sum(vx_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                    snapvy(i, j, k) = &
                                        0.5*sum(vy_hxiyiz(i:i + 1, j, k)) &
                                        + 0.5*sum(vy_ixhyiz(i, j:j + 1, k)) &
                                        + 0.5*sum(vy_ixiyhz(i, j, k:k + 1)) &
                                        + 0.125*sum(vy_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                    snapvz(i, j, k) = &
                                        0.5*sum(vz_hxiyiz(i:i + 1, j, k)) &
                                        + 0.5*sum(vz_ixhyiz(i, j:j + 1, k)) &
                                        + 0.5*sum(vz_ixiyhz(i, j, k:k + 1)) &
                                        + 0.125*sum(vz_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        call reduce_array_group(snapvx)
                        call reduce_array_group(snapvy)
                        call reduce_array_group(snapvz)

                        snapvx = 0.25*snapvx
                        snapvy = 0.25*snapvy
                        snapvz = 0.25*snapvz

                        if (rankid_group == 0) then

                            open(3, file=tidy(dir_snapshot)//'/shot_' &
                                //num2str(sgmtr%id) &
                                //'_forward_wavefield_' &
                                //num2str(l)//'_'//num2str((nx + 2*pml)*(ny + 2*pml)*(nz + pml))//'x6.bin', &
                                form='unformatted', access='stream')
                            write(3) &
                                meshgrid([nz + pml, ny + 2*pml, nx + 2*pml], [1.0, dy, dx], [0.0, oy - pml*dy, ox - pml*dx], dim=3), &
                                meshgrid([nz + pml, ny + 2*pml, nx + 2*pml], [1.0, dy, dx], [0.0, oy - pml*dy, ox - pml*dx], dim=2), &
                                zz_i(-pml + 1:nx + pml, -pml + 1:ny + pml, 1:nz + pml), &
                                snapvx(:, :, 1:nz + pml), &
                                snapvy(:, :, 1:nz + pml), &
                                snapvz(:, :, 1:nz + pml)
                            close(3)

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

                        snapvx = 0.0
                        snapvy = 0.0
                        snapvz = 0.0

                        !$omp parallel do private(i, j, k) collapse(3)
                        do k = nz1, nz2
                            do j = ny1, ny2
                                do i = nx1, nx2
                                    snapvx(i, j, k) = &
                                        0.5*sum(vx_hxiyiz(i:i + 1, j, k)) &
                                        + 0.5*sum(vx_ixhyiz(i, j:j + 1, k)) &
                                        + 0.5*sum(vx_ixiyhz(i, j, k:k + 1)) &
                                        + 0.125*sum(vx_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                    snapvy(i, j, k) = &
                                        0.5*sum(vy_hxiyiz(i:i + 1, j, k)) &
                                        + 0.5*sum(vy_ixhyiz(i, j:j + 1, k)) &
                                        + 0.5*sum(vy_ixiyhz(i, j, k:k + 1)) &
                                        + 0.125*sum(vy_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                    snapvz(i, j, k) = &
                                        0.5*sum(vz_hxiyiz(i:i + 1, j, k)) &
                                        + 0.5*sum(vz_ixhyiz(i, j:j + 1, k)) &
                                        + 0.5*sum(vz_ixiyhz(i, j, k:k + 1)) &
                                        + 0.125*sum(vz_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        call reduce_array_group(snapvx)
                        call reduce_array_group(snapvy)
                        call reduce_array_group(snapvz)

                        snapvx = 0.25*snapvx
                        snapvy = 0.25*snapvy
                        snapvz = 0.25*snapvz

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
                    stressxx_ixiyiz, stressyy_ixiyiz, stresszz_ixiyiz, &
                    stressxy_ixiyiz, stressxz_ixiyiz, stressyz_ixiyiz, &
                    memory_pdxvx_ixiyiz, &
                    memory_pdxvy_ixiyiz, &
                    memory_pdxvz_ixiyiz, &
                    memory_pdyvx_ixiyiz, &
                    memory_pdyvy_ixiyiz, &
                    memory_pdyvz_ixiyiz, &
                    memory_pdzvx_ixiyiz, &
                    memory_pdzvy_ixiyiz, &
                    memory_pdzvz_ixiyiz, &
                    stressxx_hxhyiz, stressyy_hxhyiz, stresszz_hxhyiz, &
                    stressxy_hxhyiz, stressxz_hxhyiz, stressyz_hxhyiz, &
                    memory_pdxvx_hxhyiz, &
                    memory_pdxvy_hxhyiz, &
                    memory_pdxvz_hxhyiz, &
                    memory_pdyvx_hxhyiz, &
                    memory_pdyvy_hxhyiz, &
                    memory_pdyvz_hxhyiz, &
                    memory_pdzvx_hxhyiz, &
                    memory_pdzvy_hxhyiz, &
                    memory_pdzvz_hxhyiz, &
                    stressxx_hxiyhz, stressyy_hxiyhz, stresszz_hxiyhz, &
                    stressxy_hxiyhz, stressxz_hxiyhz, stressyz_hxiyhz, &
                    memory_pdxvx_hxiyhz, &
                    memory_pdxvy_hxiyhz, &
                    memory_pdxvz_hxiyhz, &
                    memory_pdyvx_hxiyhz, &
                    memory_pdyvy_hxiyhz, &
                    memory_pdyvz_hxiyhz, &
                    memory_pdzvx_hxiyhz, &
                    memory_pdzvy_hxiyhz, &
                    memory_pdzvz_hxiyhz, &
                    stressxx_ixhyhz, stressyy_ixhyhz, stresszz_ixhyhz, &
                    stressxy_ixhyhz, stressxz_ixhyhz, stressyz_ixhyhz, &
                    memory_pdxvx_ixhyhz, &
                    memory_pdxvy_ixhyhz, &
                    memory_pdxvz_ixhyhz, &
                    memory_pdyvx_ixhyhz, &
                    memory_pdyvy_ixhyhz, &
                    memory_pdyvz_ixhyhz, &
                    memory_pdzvx_ixhyhz, &
                    memory_pdzvy_ixhyhz, &
                    memory_pdzvz_ixhyhz, &
                    vx_hxiyiz, vy_hxiyiz, vz_hxiyiz, &
                    memory_pdxxx_hxiyiz, &
                    memory_pdxxy_hxiyiz, &
                    memory_pdxxz_hxiyiz, &
                    memory_pdyxy_hxiyiz, &
                    memory_pdyyy_hxiyiz, &
                    memory_pdyyz_hxiyiz, &
                    memory_pdzxx_hxiyiz, memory_pdzxy_hxiyiz, memory_pdzxz_hxiyiz, &
                    memory_pdzyy_hxiyiz, memory_pdzyz_hxiyiz, &
                    memory_pdzzz_hxiyiz, &
                    vx_ixhyiz, vy_ixhyiz, vz_ixhyiz, &
                    memory_pdxxx_ixhyiz, &
                    memory_pdxxy_ixhyiz, &
                    memory_pdxxz_ixhyiz, &
                    memory_pdyxy_ixhyiz, &
                    memory_pdyyy_ixhyiz, &
                    memory_pdyyz_ixhyiz, &
                    memory_pdzxx_ixhyiz, memory_pdzxy_ixhyiz, memory_pdzxz_ixhyiz, &
                    memory_pdzyy_ixhyiz, memory_pdzyz_ixhyiz, &
                    memory_pdzzz_ixhyiz, &
                    vx_ixiyhz, vy_ixiyhz, vz_ixiyhz, &
                    memory_pdxxx_ixiyhz, &
                    memory_pdxxy_ixiyhz, &
                    memory_pdxxz_ixiyhz, &
                    memory_pdyxy_ixiyhz, &
                    memory_pdyyy_ixiyhz, &
                    memory_pdyyz_ixiyhz, &
                    memory_pdzxx_ixiyhz, memory_pdzxy_ixiyhz, memory_pdzxz_ixiyhz, &
                    memory_pdzyy_ixiyhz, memory_pdzyz_ixiyhz, &
                    memory_pdzzz_ixiyhz, &
                    vx_hxhyhz, vy_hxhyhz, vz_hxhyhz, &
                    memory_pdxxx_hxhyhz, &
                    memory_pdxxy_hxhyhz, &
                    memory_pdxxz_hxhyhz, &
                    memory_pdyxy_hxhyhz, &
                    memory_pdyyy_hxhyhz, &
                    memory_pdyyz_hxhyhz, &
                    memory_pdzxx_hxhyhz, memory_pdzxy_hxhyhz, memory_pdzxz_hxhyhz, &
                    memory_pdzyy_hxhyhz, memory_pdzyz_hxhyhz, &
                    memory_pdzzz_hxhyhz)
            else
                call update_wavefield(dt, &
                    stressxx_ixiyiz, stressyy_ixiyiz, stresszz_ixiyiz, &
                    stressxy_ixiyiz, stressxz_ixiyiz, stressyz_ixiyiz, &
                    memory_pdxvx_ixiyiz, &
                    memory_pdxvy_ixiyiz, &
                    memory_pdxvz_ixiyiz, &
                    memory_pdyvx_ixiyiz, &
                    memory_pdyvy_ixiyiz, &
                    memory_pdyvz_ixiyiz, &
                    memory_pdzvx_ixiyiz, &
                    memory_pdzvy_ixiyiz, &
                    memory_pdzvz_ixiyiz, &
                    stressxx_hxhyiz, stressyy_hxhyiz, stresszz_hxhyiz, &
                    stressxy_hxhyiz, stressxz_hxhyiz, stressyz_hxhyiz, &
                    memory_pdxvx_hxhyiz, &
                    memory_pdxvy_hxhyiz, &
                    memory_pdxvz_hxhyiz, &
                    memory_pdyvx_hxhyiz, &
                    memory_pdyvy_hxhyiz, &
                    memory_pdyvz_hxhyiz, &
                    memory_pdzvx_hxhyiz, &
                    memory_pdzvy_hxhyiz, &
                    memory_pdzvz_hxhyiz, &
                    stressxx_hxiyhz, stressyy_hxiyhz, stresszz_hxiyhz, &
                    stressxy_hxiyhz, stressxz_hxiyhz, stressyz_hxiyhz, &
                    memory_pdxvx_hxiyhz, &
                    memory_pdxvy_hxiyhz, &
                    memory_pdxvz_hxiyhz, &
                    memory_pdyvx_hxiyhz, &
                    memory_pdyvy_hxiyhz, &
                    memory_pdyvz_hxiyhz, &
                    memory_pdzvx_hxiyhz, &
                    memory_pdzvy_hxiyhz, &
                    memory_pdzvz_hxiyhz, &
                    stressxx_ixhyhz, stressyy_ixhyhz, stresszz_ixhyhz, &
                    stressxy_ixhyhz, stressxz_ixhyhz, stressyz_ixhyhz, &
                    memory_pdxvx_ixhyhz, &
                    memory_pdxvy_ixhyhz, &
                    memory_pdxvz_ixhyhz, &
                    memory_pdyvx_ixhyhz, &
                    memory_pdyvy_ixhyhz, &
                    memory_pdyvz_ixhyhz, &
                    memory_pdzvx_ixhyhz, &
                    memory_pdzvy_ixhyhz, &
                    memory_pdzvz_ixhyhz, &
                    vx_hxiyiz, vy_hxiyiz, vz_hxiyiz, &
                    memory_pdxxx_hxiyiz, &
                    memory_pdxxy_hxiyiz, &
                    memory_pdxxz_hxiyiz, &
                    memory_pdyxy_hxiyiz, &
                    memory_pdyyy_hxiyiz, &
                    memory_pdyyz_hxiyiz, &
                    memory_pdzxz_hxiyiz, &
                    memory_pdzyz_hxiyiz, &
                    memory_pdzzz_hxiyiz, &
                    vx_ixhyiz, vy_ixhyiz, vz_ixhyiz, &
                    memory_pdxxx_ixhyiz, &
                    memory_pdxxy_ixhyiz, &
                    memory_pdxxz_ixhyiz, &
                    memory_pdyxy_ixhyiz, &
                    memory_pdyyy_ixhyiz, &
                    memory_pdyyz_ixhyiz, &
                    memory_pdzxz_ixhyiz, &
                    memory_pdzyz_ixhyiz, &
                    memory_pdzzz_ixhyiz, &
                    vx_ixiyhz, vy_ixiyhz, vz_ixiyhz, &
                    memory_pdxxx_ixiyhz, &
                    memory_pdxxy_ixiyhz, &
                    memory_pdxxz_ixiyhz, &
                    memory_pdyxy_ixiyhz, &
                    memory_pdyyy_ixiyhz, &
                    memory_pdyyz_ixiyhz, &
                    memory_pdzxz_ixiyhz, &
                    memory_pdzyz_ixiyhz, &
                    memory_pdzzz_ixiyhz, &
                    vx_hxhyhz, vy_hxhyhz, vz_hxhyhz, &
                    memory_pdxxx_hxhyhz, &
                    memory_pdxxy_hxhyhz, &
                    memory_pdxxz_hxhyhz, &
                    memory_pdyxy_hxhyhz, &
                    memory_pdyyy_hxhyhz, &
                    memory_pdyyz_hxhyhz, &
                    memory_pdzxz_hxhyhz, &
                    memory_pdzyz_hxhyhz, &
                    memory_pdzzz_hxhyhz)
            end if

            ! Source term
            call add_source(t)

            if (verbose .and. (mod(t, max(nint(nt/10.0), 1)) == 0 .or. t == 1 .or. t == nt)) then

                wnan = group_and(any(isnan(vx_hxiyiz)) .or. any(isnan(vy_hxiyiz)) .or. any(isnan(vz_hxiyiz)))

                wmin1 = group_min(min(minval(vx_hxiyiz), minval(vx_ixhyiz), minval(vx_ixiyhz), minval(vx_hxhyhz)))
                wmax1 = group_max(max(maxval(vx_hxiyiz), maxval(vx_ixhyiz), maxval(vx_ixiyhz), maxval(vx_hxhyhz)))

                wmin2 = group_min(min(minval(vy_hxiyiz), minval(vy_ixhyiz), minval(vy_ixiyhz), minval(vy_hxhyhz)))
                wmax2 = group_max(max(maxval(vy_hxiyiz), maxval(vy_ixhyiz), maxval(vy_ixiyhz), maxval(vy_hxhyhz)))

                wmin3 = group_min(min(minval(vz_hxiyiz), minval(vz_ixhyiz), minval(vz_ixiyhz), minval(vz_hxhyhz)))
                wmax3 = group_max(max(maxval(vz_hxiyiz), maxval(vz_ixhyiz), maxval(vz_ixiyhz), maxval(vz_hxhyhz)))

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

        ! Close boundary saving files
        if (yn_reconstruct) then
            call close_boundary_saving
        end if

        call allreduce_array_group(seis_vx)
        call allreduce_array_group(seis_vy)
        call allreduce_array_group(seis_vz)
        do i = 1, sgmtr%nr
            if (.not. (i >= trace_range(rankid_group, 1) .and. i <= trace_range(rankid_group, 2))) then
                seis_vx(:, i) = 0.0
                seis_vy(:, i) = 0.0
                seis_vz(:, i) = 0.0
            end if
        end do

        ! Rotate data to local coordinates if necessary.
        ! The rotation follows the formula developed by
        !   Hestholm and Ruud, 1998, Geophysics
        !   3-D finite-difference elastic wave modeling including surface topography
        !   Equations C-7 -- C-9
        if (this%receiver_vertical_to_surface) then

            rec_slopex = ginterp(topox, topoy, flatten(slopex_ixhy), sgmtr%recr(:)%x, sgmtr%recr(:)%y)
            rec_slopex = atan(rec_slopex)
            rec_slopey = ginterp(topox, topoy, flatten(slopey_ixhy), sgmtr%recr(:)%x, sgmtr%recr(:)%y)
            rec_slopey = atan(rec_slopey*cos(rec_slopex))

            ! z positive direction points downward
            rec_slopex = -rec_slopex
            rec_slopey = -rec_slopey

            pt = zeros(3)
            !$omp parallel do private(i, t, pt) collapse(2)
            do t = 1, nt
                do i = trace_range(rankid_group, 1), trace_range(rankid_group, 2)

                    pt = matx(reshape([ &
                        cos(rec_slopex(i)), &
                        -sin(rec_slopex(i))*sin(rec_slopey(i)), &
                        -sin(rec_slopex(i))*cos(rec_slopey(i)), &
                        0.0, &
                        cos(rec_slopey(i)), &
                        -sin(rec_slopey(i)), &
                        sin(rec_slopex(i)), &
                        cos(rec_slopex(i))*sin(rec_slopey), &
                        cos(rec_slopex(i))*cos(rec_slopey)], [3, 3]), &
                        [seis_vx(t, ir), seis_vy(t, ir), seis_vz(t, ir)])

                    seis_vx(t, ir) = pt(1)
                    seis_vy(t, ir) = pt(2)
                    seis_vz(t, ir) = pt(3)

                end do
            end do
            !$omp end parallel do

        end if

        ! resample if necessary
        if (yn_compx) then
            call this%seis_vx%init(nt=nt, dt=dt, nr=sgmtr%nr)
            call this%seis_vx%from_array(seis_vx)
            call this%seis_vx%resamp(nnt=data_nt, ddt=data_dt)
            call this%seis_vx%collect_group
            if (rankid_group == 0) then
                call this%seis_vx%output(tidy(dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_seismogram_x.su')
            end if
        end if

        if (yn_compy) then
            call this%seis_vy%init(nt=nt, dt=dt, nr=sgmtr%nr)
            call this%seis_vy%from_array(seis_vy)
            call this%seis_vy%resamp(nnt=data_nt, ddt=data_dt)
            call this%seis_vy%collect_group
            if (rankid_group == 0) then
                call this%seis_vy%output(tidy(dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_seismogram_y.su')
            end if
        end if

        if (yn_compz) then
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
