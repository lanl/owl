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


module acoustic_iso_2d_wavefield

    use libflit
    use acoustic_iso_2d_vars
    use acoustic_iso_2d_cfspml

    implicit none

    ! FD stencils
#include 'macro_fd_stencil.f90'

    ! Average rho on a half_x or half_z node
#define rho_eff_x (0.5*(rho(i + 1, j) + rho(i, j)))
#define rho_eff_z (0.5*(rho(i, j + 1) + rho(i, j)))

contains

    !
    !> Update model in model with CFS-PMLs on four boundaries
    !
    subroutine update_wavefield(dt, &
            p, vx, vz, &
            memory_pdxvx, memory_pdzvz, memory_pdxp, memory_pdzp)

        real, intent(in) :: dt
        real, allocatable, dimension(:, :), intent(inout) :: &
            p, vx, vz, memory_pdxvx, memory_pdzvz, memory_pdxp, memory_pdzp

        integer :: i, j
        real :: pdxvx, pdzvz, pdxp, pdzp

        !$omp parallel do private(i, j, pdxvx, pdzvz) collapse(2) schedule(auto)
        do j = -pml + 1, nz + pml
            do i = -pml + 1, nx + pml

                pdxvx = idx*pdxvx_stencil
                pdzvz = idz*pdzvz_stencil

                if (i <= 0 .or. i >= nx + 1) then
                    memory_pdxvx(i, j) = axi(i)*pdxvx + bxi(i)*memory_pdxvx(i, j)
                    pdxvx = (pdxvx + memory_pdxvx(i, j))/kxi(i)
                end if
                if (j <= 0 .or. j >= nz + 1) then
                    memory_pdzvz(i, j) = azi(j)*pdzvz + bzi(j)*memory_pdzvz(i, j)
                    pdzvz = (pdzvz + memory_pdzvz(i, j))/kzi(j)
                end if

                p(i, j) = p(i, j) - dt*bk(i, j)*(pdxvx + pdzvz)

            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, pdxp, pdzp) collapse(2) schedule(auto)
        do j = -pml + 1, nz + pml
            do i = -pml, nx + pml - 1

                pdxp = idx*pdxp_stencil

                if (i + 1 <= 1 .or. i + 1 >= nx + 1) then
                    memory_pdxp(i + 1, j) = axh(i + 1)*pdxp + bxh(i + 1)*memory_pdxp(i + 1, j)
                    pdxp = (pdxp + memory_pdxp(i + 1, j))/kxh(i + 1)
                end if

                vx(i + 1, j) = vx(i + 1, j) - dt/rho_eff_x*pdxp

            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, pdxp, pdzp) collapse(2) schedule(auto)
        do j = -pml, nz + pml - 1
            do i = -pml + 1, nx + pml

                pdzp = idz*pdzp_stencil

                if (j + 1 <= 1 .or. j + 1 >= nz + 1) then
                    memory_pdzp(i, j + 1) = azh(j + 1)*pdzp + bzh(j + 1)*memory_pdzp(i, j + 1)
                    pdzp = (pdzp + memory_pdzp(i, j + 1))/kzh(j + 1)
                end if

                vz(i, j + 1) = vz(i, j + 1) - dt/rho_eff_z*pdzp

            end do
        end do
        !$omp end parallel do

    end subroutine update_wavefield

    !
    !> Update wavefield in model with free surface
    !
    subroutine update_wavefield_free_surface(dt, &
            p, vx, vz, &
            memory_pdxvx, memory_pdzvz, memory_pdxp, memory_pdzp)

        real, intent(in) :: dt
        real, allocatable, dimension(:, :), intent(inout) :: &
            p, vx, vz, memory_pdxvx, memory_pdzvz, memory_pdxp, memory_pdzp

        integer :: i, j
        real :: pdxvx, pdzvz, pdxp, pdzp

        ! The pressure is odd about the free surface (row 1), so vz is even: mirror vz above the
        ! surface, vz(-z) = vz(z). The z-derivative stencil below reads up to fdhalf rows above the
        ! surface; leaving zeros there makes vz discontinuous at the surface, an error that does not
        ! vanish with grid refinement.
        !$omp parallel do private(j) schedule(auto)
        do j = 1, fdhalf
            vz(:, 2 - j) = vz(:, 1 + j)
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, pdxvx, pdzvz) collapse(2) schedule(auto)
        do j = 2, nz + pml
            do i = -pml + 1, nx + pml

                pdxvx = idx*pdxvx_stencil
                pdzvz = idz*pdzvz_stencil

                if (i <= 0 .or. i >= nx + 1) then
                    memory_pdxvx(i, j) = axi(i)*pdxvx + bxi(i)*memory_pdxvx(i, j)
                    pdxvx = (pdxvx + memory_pdxvx(i, j))/kxi(i)
                end if
                if (j >= nz + 1) then
                    memory_pdzvz(i, j) = azi(j)*pdzvz + bzi(j)*memory_pdzvz(i, j)
                    pdzvz = (pdzvz + memory_pdzvz(i, j))/kzi(j)
                end if

                p(i, j) = p(i, j) - dt*bk(i, j)*(pdxvx + pdzvz)

            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j) collapse(2) schedule(auto)
        do j = -pml + 1, 0
            do i = -pml + 1, nx + pml
                p(i, j) = -p(i, 2 - j)
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, pdxp, pdzp) collapse(2) schedule(auto)
        do j = 1, nz + pml
            do i = -pml, nx + pml - 1

                pdxp = idx*pdxp_stencil

                if (i + 1 <= 1 .or. i + 1 >= nx + 1) then
                    memory_pdxp(i + 1, j) = axh(i + 1)*pdxp + bxh(i + 1)*memory_pdxp(i + 1, j)
                    pdxp = (pdxp + memory_pdxp(i + 1, j))/kxh(i + 1)
                end if

                vx(i + 1, j) = vx(i + 1, j) - dt/rho_eff_x*pdxp

            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, pdxp, pdzp) collapse(2) schedule(auto)
        do j = 1, nz + pml - 1
            do i = -pml + 1, nx + pml

                pdzp = idz*pdzp_stencil

                if (j + 1 >= nz + 1) then
                    memory_pdzp(i, j + 1) = azh(j + 1)*pdzp + bzh(j + 1)*memory_pdzp(i, j + 1)
                    pdzp = (pdzp + memory_pdzp(i, j + 1))/kzh(j + 1)
                end if

                vz(i, j + 1) = vz(i, j + 1) - dt/rho_eff_z*pdzp

            end do
        end do
        !$omp end parallel do

    end subroutine update_wavefield_free_surface

    !
    !> Add source term
    !
    subroutine add_source(t)

        integer, intent(in) :: t

        integer :: k
        integer :: sgx, sgz, nbeg, nend, shx, shz
        real :: amp, polar, a

        do k = 1, sgmtr%ns

            nbeg = nint(sgmtr%srcr(k)%t0/dt) + 1
            nend = nbeg + sgmtr%srcr(k)%nt - 1

            sgx = sgmtr%srcr(k)%gx
            sgz = sgmtr%srcr(k)%gz
            shx = sgmtr%srcr(k)%hx
            shz = sgmtr%srcr(k)%hz

            if (t >= nbeg .and. t <= nend) then

                amp = sgmtr%srcr(k)%stf(t - nbeg + 1)*sgmtr%srcr(k)%amp*dt

                select case (sgmtr%srcr(k)%mechanism)

                    case ('force')
                        ! Force vector. Under a free surface (p = 0) the image of a horizontal force is
                        ! opposite to it (vx is odd about the surface) and the image of a vertical
                        ! force equals it (vz is even about the surface)
                        polar = sgmtr%srcr(k)%polar
                        a = sin(polar)*amp/rho(shx, sgz)
                        call add_source_value_2d(vx, shx, sgz, sgmtr%srcr(k)%interp_hx, sgmtr%srcr(k)%interp_iz, &
                            a, -a, .false., yn_free_surface)
                        a = cos(polar)*amp/rho(sgx, shz)
                        call add_source_value_2d(vz, sgx, shz, sgmtr%srcr(k)%interp_ix, sgmtr%srcr(k)%interp_hz, &
                            a, a, .true., yn_free_surface)

                    case ('explosion')
                        ! Explosive source. Under a free surface its image is opposite to it (p is odd
                        ! about the surface), so a source on the surface radiates nothing
                        call add_source_value_2d(p, sgx, sgz, sgmtr%srcr(k)%interp_ix, sgmtr%srcr(k)%interp_iz, &
                            amp, -amp, .false., yn_free_surface)

                end select

            end if

        end do

    end subroutine add_source

end module
