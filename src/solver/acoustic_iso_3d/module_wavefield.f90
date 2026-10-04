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

module acoustic_iso_3d_wavefield

    use libflit
    use acoustic_iso_3d_vars
    use acoustic_iso_3d_cfspml

    implicit none

    ! FD stencil
#include 'macro_fd_stencil.f90'

    ! Average rho on a half_x, half_y or half_z node
#define rho_eff_x (0.5*(rho(i + 1, j, k) + rho(i, j, k)))
#define rho_eff_y (0.5*(rho(i, j + 1, k) + rho(i, j, k)))
#define rho_eff_z (0.5*(rho(i, j, k + 1) + rho(i, j, k)))

contains

    !
    !> Update wavefield
    !
    subroutine update_wavefield(dt, &
            p, vx, vy, vz, &
            memory_pdxp_xmin, memory_pdxp_xmax, &
            memory_pdyp_ymin, memory_pdyp_ymax, &
            memory_pdzp_zmin, memory_pdzp_zmax, &
            memory_pdxvx_xmin, memory_pdxvx_xmax, &
            memory_pdyvy_ymin, memory_pdyvy_ymax, &
            memory_pdzvz_zmin, memory_pdzvz_zmax)

        real, intent(in) :: dt
        real, allocatable, dimension(:, :, :), intent(inout) :: p, vx, vy, vz
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxp_xmin, memory_pdxp_xmax
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdyp_ymin, memory_pdyp_ymax
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdzp_zmin, memory_pdzp_zmax
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxvx_xmin, memory_pdxvx_xmax
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdyvy_ymin, memory_pdyvy_ymax
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdzvz_zmin, memory_pdzvz_zmax

        integer :: i, j, k
        real :: pdxvx, pdyvy, pdzvz
        real :: pdxp, pdyp, pdzp

        call commute_array_group(vx, fdhalf)
        call commute_array_group(vy, fdhalf)
        call commute_array_group(vz, fdhalf)

        !$omp parallel do private(i, j, k, pdxvx, pdyvy, pdzvz) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1, nx2

                    pdxvx = idx*pdxvx_stencil
                    pdyvy = idy*pdyvy_stencil
                    pdzvz = idz*pdzvz_stencil

                    if (i <= 0) then
                        memory_pdxvx_xmin(i, j, k) = axi(i)*pdxvx + bxi(i)*memory_pdxvx_xmin(i, j, k)
                        pdxvx = (pdxvx + memory_pdxvx_xmin(i, j, k))/kxi(i)
                    else if (i >= nx + 1) then
                        memory_pdxvx_xmax(i, j, k) = axi(i)*pdxvx + bxi(i)*memory_pdxvx_xmax(i, j, k)
                        pdxvx = (pdxvx + memory_pdxvx_xmax(i, j, k))/kxi(i)
                    end if
                    if (j <= 0) then
                        memory_pdyvy_ymin(i, j, k) = ayi(j)*pdyvy + byi(j)*memory_pdyvy_ymin(i, j, k)
                        pdyvy = (pdyvy + memory_pdyvy_ymin(i, j, k))/kyi(j)
                    else if (j >= ny + 1) then
                        memory_pdyvy_ymax(i, j, k) = ayi(j)*pdyvy + byi(j)*memory_pdyvy_ymax(i, j, k)
                        pdyvy = (pdyvy + memory_pdyvy_ymax(i, j, k))/kyi(j)
                    end if
                    if (k <= 0) then
                        memory_pdzvz_zmin(i, j, k) = azi(k)*pdzvz + bzi(k)*memory_pdzvz_zmin(i, j, k)
                        pdzvz = (pdzvz + memory_pdzvz_zmin(i, j, k))/kzi(k)
                    else if (k >= nz + 1) then
                        memory_pdzvz_zmax(i, j, k) = azi(k)*pdzvz + bzi(k)*memory_pdzvz_zmax(i, j, k)
                        pdzvz = (pdzvz + memory_pdzvz_zmax(i, j, k))/kzi(k)
                    end if

                    p(i, j, k) = p(i, j, k) - dt*bk(i, j, k)*(pdxvx + pdyvy + pdzvz)

                end do
            end do
        end do
        !$omp end parallel do

        call commute_array_group(p, fdhalf)

        !$omp parallel do private(i, j, k, pdxp, pdyp, pdzp) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1 - 1, nx2 - 1

                    pdxp = idx*pdxp_stencil

                    if (i + 1 <= 1) then
                        memory_pdxp_xmin(i + 1, j, k) = axh(i + 1)*pdxp + bxh(i + 1)*memory_pdxp_xmin(i + 1, j, k)
                        pdxp = (pdxp + memory_pdxp_xmin(i + 1, j, k))/kxh(i + 1)
                    else if (i + 1 >= nx + 1) then
                        memory_pdxp_xmax(i + 1, j, k) = axh(i + 1)*pdxp + bxh(i + 1)*memory_pdxp_xmax(i + 1, j, k)
                        pdxp = (pdxp + memory_pdxp_xmax(i + 1, j, k))/kxh(i + 1)
                    end if

                    vx(i + 1, j, k) = vx(i + 1, j, k) - dt/rho_eff_x*pdxp

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, pdxp, pdyp, pdzp) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1 - 1, ny2 - 1
                do i = nx1, nx2

                    pdyp = idy*pdyp_stencil

                    if (j + 1 <= 1) then
                        memory_pdyp_ymin(i, j + 1, k) = ayh(j + 1)*pdyp + byh(j + 1)*memory_pdyp_ymin(i, j + 1, k)
                        pdyp = (pdyp + memory_pdyp_ymin(i, j + 1, k))/kyh(j + 1)
                    else if (j + 1 >= ny + 1) then
                        memory_pdyp_ymax(i, j + 1, k) = ayh(j + 1)*pdyp + byh(j + 1)*memory_pdyp_ymax(i, j + 1, k)
                        pdyp = (pdyp + memory_pdyp_ymax(i, j + 1, k))/kyh(j + 1)
                    end if

                    vy(i, j + 1, k) = vy(i, j + 1, k) - dt/rho_eff_y*pdyp

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, pdxp, pdyp, pdzp) collapse(3) schedule(auto)
        do k = nz1 - 1, nz2 - 1
            do j = ny1, ny2
                do i = nx1, nx2

                    pdzp = idz*pdzp_stencil

                    if (k + 1 <= 1) then
                        memory_pdzp_zmin(i, j, k + 1) = azh(k + 1)*pdzp + bzh(k + 1)*memory_pdzp_zmin(i, j, k + 1)
                        pdzp = (pdzp + memory_pdzp_zmin(i, j, k + 1))/kzh(k + 1)
                    else if (k + 1 >= nz + 1) then
                        memory_pdzp_zmax(i, j, k + 1) = azh(k + 1)*pdzp + bzh(k + 1)*memory_pdzp_zmax(i, j, k + 1)
                        pdzp = (pdzp + memory_pdzp_zmax(i, j, k + 1))/kzh(k + 1)
                    end if

                    vz(i, j, k + 1) = vz(i, j, k + 1) - dt/rho_eff_z*pdzp

                end do
            end do
        end do
        !$omp end parallel do

    end subroutine update_wavefield

    !
    !> Update wavefield
    !
    subroutine update_wavefield_free_surface(dt, &
            p, vx, vy, vz, &
            memory_pdxp_xmin, memory_pdxp_xmax, &
            memory_pdyp_ymin, memory_pdyp_ymax, &
            memory_pdzp_zmax, &
            memory_pdxvx_xmin, memory_pdxvx_xmax, &
            memory_pdyvy_ymin, memory_pdyvy_ymax, &
            memory_pdzvz_zmax)

        real, intent(in) :: dt
        real, allocatable, dimension(:, :, :), intent(inout) :: p, vx, vy, vz
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxp_xmin, memory_pdxp_xmax
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdyp_ymin, memory_pdyp_ymax
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdzp_zmax
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxvx_xmin, memory_pdxvx_xmax
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdyvy_ymin, memory_pdyvy_ymax
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdzvz_zmax

        integer :: i, j, k
        real :: pdxvx, pdyvy, pdzvz
        real :: pdxp, pdyp, pdzp

        call commute_array_group(vx, fdhalf)
        call commute_array_group(vy, fdhalf)
        call commute_array_group(vz, fdhalf)

        ! Update p
        !$omp parallel do private(i, j, k, pdxvx, pdyvy, pdzvz) collapse(3) schedule(auto)
        do k = max(2, nz1), nz2
            do j = ny1, ny2
                do i = nx1, nx2

                    pdxvx = idx*pdxvx_stencil
                    pdyvy = idy*pdyvy_stencil
                    pdzvz = idz*pdzvz_stencil

                    if (i <= 0) then
                        memory_pdxvx_xmin(i, j, k) = axi(i)*pdxvx + bxi(i)*memory_pdxvx_xmin(i, j, k)
                        pdxvx = (pdxvx + memory_pdxvx_xmin(i, j, k))/kxi(i)
                    else if (i >= nx + 1) then
                        memory_pdxvx_xmax(i, j, k) = axi(i)*pdxvx + bxi(i)*memory_pdxvx_xmax(i, j, k)
                        pdxvx = (pdxvx + memory_pdxvx_xmax(i, j, k))/kxi(i)
                    end if
                    if (j <= 0) then
                        memory_pdyvy_ymin(i, j, k) = ayi(j)*pdyvy + byi(j)*memory_pdyvy_ymin(i, j, k)
                        pdyvy = (pdyvy + memory_pdyvy_ymin(i, j, k))/kyi(j)
                    else if (j >= ny + 1) then
                        memory_pdyvy_ymax(i, j, k) = ayi(j)*pdyvy + byi(j)*memory_pdyvy_ymax(i, j, k)
                        pdyvy = (pdyvy + memory_pdyvy_ymax(i, j, k))/kyi(j)
                    end if
                    if (k >= nz + 1) then
                        memory_pdzvz_zmax(i, j, k) = azi(k)*pdzvz + bzi(k)*memory_pdzvz_zmax(i, j, k)
                        pdzvz = (pdzvz + memory_pdzvz_zmax(i, j, k))/kzi(k)
                    end if

                    p(i, j, k) = p(i, j, k) - dt*bk(i, j, k)*(pdxvx + pdyvy + pdzvz)

                end do
            end do
        end do
        !$omp end parallel do

        ! Mirroring p to above-free-surface region
        if (nz1 <= 0) then
            !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
            do k = max(-pml + 1, nz1), min(0, nz2)
                do j = ny1, ny2
                    do i = nx1, nx2
                        p(i, j, k) = -p(i, j, 2 - k)
                    end do
                end do
            end do
            !$omp end parallel do
        end if

        call commute_array_group(p, fdhalf)

        !$omp parallel do private(i, j, k, pdxp, pdyp, pdzp) collapse(3) schedule(auto)
        do k = max(1, nz1), nz2
            do j = ny1, ny2
                do i = nx1 - 1, nx2 - 1

                    pdxp = idx*pdxp_stencil

                    if (i + 1 <= 1) then
                        memory_pdxp_xmin(i + 1, j, k) = axh(i + 1)*pdxp + bxh(i + 1)*memory_pdxp_xmin(i + 1, j, k)
                        pdxp = (pdxp + memory_pdxp_xmin(i + 1, j, k))/kxh(i + 1)
                    else if (i + 1 >= nx + 1) then
                        memory_pdxp_xmax(i + 1, j, k) = axh(i + 1)*pdxp + bxh(i + 1)*memory_pdxp_xmax(i + 1, j, k)
                        pdxp = (pdxp + memory_pdxp_xmax(i + 1, j, k))/kxh(i + 1)
                    end if

                    vx(i + 1, j, k) = vx(i + 1, j, k) - dt/rho_eff_x*pdxp

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, pdxp, pdyp, pdzp) collapse(3) schedule(auto)
        do k = max(1, nz1), nz2
            do j = ny1 - 1, ny2 - 1
                do i = nx1, nx2

                    pdyp = idy*pdyp_stencil

                    if (j + 1 <= 1) then
                        memory_pdyp_ymin(i, j + 1, k) = ayh(j + 1)*pdyp + byh(j + 1)*memory_pdyp_ymin(i, j + 1, k)
                        pdyp = (pdyp + memory_pdyp_ymin(i, j + 1, k))/kyh(j + 1)
                    else if (j + 1 >= ny + 1) then
                        memory_pdyp_ymax(i, j + 1, k) = ayh(j + 1)*pdyp + byh(j + 1)*memory_pdyp_ymax(i, j + 1, k)
                        pdyp = (pdyp + memory_pdyp_ymax(i, j + 1, k))/kyh(j + 1)
                    end if

                    vy(i, j + 1, k) = vy(i, j + 1, k) - dt/rho_eff_y*pdyp

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, pdxp, pdyp, pdzp) collapse(3) schedule(auto)
        do k = max(1, nz1 - 1), nz2 - 1
            do j = ny1, ny2
                do i = nx1, nx2

                    pdzp = idz*pdzp_stencil

                    if (k + 1 >= nz + 1) then
                        memory_pdzp_zmax(i, j, k + 1) = azh(k + 1)*pdzp + bzh(k + 1)*memory_pdzp_zmax(i, j, k + 1)
                        pdzp = (pdzp + memory_pdzp_zmax(i, j, k + 1))/kzh(k + 1)
                    end if

                    vz(i, j, k + 1) = vz(i, j, k + 1) - dt/rho_eff_z*pdzp

                end do
            end do
        end do
        !$omp end parallel do

    end subroutine update_wavefield_free_surface

    !
    !> Add source
    !
    subroutine add_source(t)

        integer, intent(in) :: t

        integer :: k
        integer :: sgx, sgy, sgz, nbeg, nend
        real :: polar, azimuth, amp
        integer :: irx, iry, irz
        real :: rho_s(1:1)

        do k = 1, sgmtr%ns

            nbeg = nint(sgmtr%srcr(k)%t0/dt) + 1
            nend = nbeg + sgmtr%srcr(k)%nt - 1

            if (t >= nbeg .and. t <= nend) then

                amp = sgmtr%srcr(k)%stf(t - nbeg + 1)*sgmtr%srcr(k)%amp*dt

                select case (sgmtr%srcr(k)%mechanism)

                    case ('force')
                        ! Force vector
                        polar = sgmtr%srcr(k)%polar
                        azimuth = sgmtr%srcr(k)%azimuth

                        sgx = sgmtr%srcr(k)%gx
                        sgy = sgmtr%srcr(k)%gy
                        sgz = sgmtr%srcr(k)%gz

                        rho_s = 0
                        if (is_in_block(sgx, sgy, sgz)) then
                            rho_s = rho(sgx, sgy, sgz)
                        end if
                        call allreduce_array_group(rho_s)
                        amp = amp/rho_s(1)

                        sgx = sgmtr%srcr(k)%hx
                        sgy = sgmtr%srcr(k)%gy
                        sgz = sgmtr%srcr(k)%gz
                        !$omp parallel do private(irx, iry, irz) collapse(3) schedule(auto)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        vx(sgx + irx, sgy + iry, sgz + irz) = &
                                            vx(sgx + irx, sgy + iry, sgz + irz) + sin(polar)*cos(azimuth)*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        sgx = sgmtr%srcr(k)%gx
                        sgy = sgmtr%srcr(k)%hy
                        sgz = sgmtr%srcr(k)%gz
                        !$omp parallel do private(irx, iry, irz) collapse(3) schedule(auto)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        vy(sgx + irx, sgy + iry, sgz + irz) = &
                                            vy(sgx + irx, sgy + iry, sgz + irz) + sin(polar)*sin(azimuth)*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        sgx = sgmtr%srcr(k)%gx
                        sgy = sgmtr%srcr(k)%gy
                        sgz = sgmtr%srcr(k)%hz
                        !$omp parallel do private(irx, iry, irz) collapse(3) schedule(auto)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        vz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vz(sgx + irx, sgy + iry, sgz + irz) + cos(polar)*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                    case ('explosion')
                        ! Explosive source

                        sgx = sgmtr%srcr(k)%gx
                        sgy = sgmtr%srcr(k)%gy
                        sgz = sgmtr%srcr(k)%gz
                        !$omp parallel do private(irx, iry, irz) collapse(3) schedule(auto)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        p(sgx + irx, sgy + iry, sgz + irz) = &
                                            p(sgx + irx, sgy + iry, sgz + irz) + amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                end select

            end if

        end do

    end subroutine add_source

end module
