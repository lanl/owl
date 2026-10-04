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

        ! The pressure is odd about the free surface (k = 1), so vz is even: mirror vz above the
        ! surface, vz(-z) = vz(z). The z-derivative stencil of the pressure update reads up to fdhalf
        ! rows above the surface; leaving zeros there makes vz discontinuous at the surface, an error
        ! that does not vanish with grid refinement. With a free surface the blocks start at k = 1,
        ! so these rows are the halo above the top block; the mirror follows the exchange so that it
        ! reads valid rows.
        !$omp parallel do private(k) schedule(auto)
        do k = 1, fdhalf
            if (2 - k >= nz1 - fdhalf .and. 1 + k <= nz2 + fdhalf) then
                vz(:, :, 2 - k) = vz(:, :, 1 + k)
            end if
        end do
        !$omp end parallel do

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

        call commute_array_group(p, fdhalf)

        ! Mirror p above the free surface, p(-z) = -p(z). With a free surface the blocks start at
        ! k = 1, so these rows are the halo above the top block; the mirror follows the exchange so
        ! that it reads valid rows also in a top block thinner than fdhalf + 1
        !$omp parallel do private(k) schedule(auto)
        do k = 1, fdhalf
            if (1 - k >= nz1 - fdhalf .and. 1 + k <= nz2 + fdhalf) then
                p(:, :, 1 - k) = -p(:, :, 1 + k)
            end if
        end do
        !$omp end parallel do

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
    !> With a free surface (p = 0), every source comes with its mirror image (see
    !> add_source_value_3d): the pressure and the horizontal particle velocities are odd about the
    !> surface (image opposite to the source), and the vertical particle velocity is even (image equal
    !> to the source). A pressure source on the surface therefore radiates nothing.
    !
    subroutine add_source(t)

        integer, intent(in) :: t

        integer :: k, nbeg, nend
        real :: polar, azimuth, amp, a
        real :: rho_s(1:1)
        integer, dimension(6) :: block

        block = [nx1_interior, nx2_interior, ny1_interior, ny2_interior, nz1_interior, nz2_interior]

        do k = 1, sgmtr%ns

            nbeg = nint(sgmtr%srcr(k)%t0/dt) + 1
            nend = nbeg + sgmtr%srcr(k)%nt - 1

            if (t >= nbeg .and. t <= nend) then

                amp = sgmtr%srcr(k)%stf(t - nbeg + 1)*sgmtr%srcr(k)%amp*dt

                associate (s => sgmtr%srcr(k))

                select case (s%mechanism)

                    case ('force')
                        ! Force vector
                        polar = s%polar
                        azimuth = s%azimuth

                        rho_s = 0
                        if (is_in_block(s%gx, s%gy, s%gz)) then
                            rho_s = rho(s%gx, s%gy, s%gz)
                        end if
                        call allreduce_array_group(rho_s)
                        amp = amp/rho_s(1)

                        a = sin(polar)*cos(azimuth)*amp
                        call add_source_value_3d(vx, s%hx, s%gy, s%gz, s%interp_hx, s%interp_iy, s%interp_iz, &
                            a, -a, .false., yn_free_surface, block)
                        a = sin(polar)*sin(azimuth)*amp
                        call add_source_value_3d(vy, s%gx, s%hy, s%gz, s%interp_ix, s%interp_hy, s%interp_iz, &
                            a, -a, .false., yn_free_surface, block)
                        a = cos(polar)*amp
                        call add_source_value_3d(vz, s%gx, s%gy, s%hz, s%interp_ix, s%interp_iy, s%interp_hz, &
                            a, a, .true., yn_free_surface, block)

                    case ('explosion')
                        ! Explosive source
                        call add_source_value_3d(p, s%gx, s%gy, s%gz, s%interp_ix, s%interp_iy, s%interp_iz, &
                            amp, -amp, .false., yn_free_surface, block)

                end select

                end associate

            end if

        end do

    end subroutine add_source

end module
