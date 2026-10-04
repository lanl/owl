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

module elastic_vhtiort_2d_wavefield

    use libflit
    use elastic_vhtiort_2d_vars
    use elastic_vhtiort_2d_cfspml

    implicit none

    ! Finite-difference stencils
#include 'macro_fd_stencil.f90'

    ! Average rho at (1/2, 0) and (0, 1/2) nodes
#define rho_eff_x (0.5*sum(rho(i:i + 1, j)))
#define rho_eff_z (0.5*sum(rho(i, j:j + 1)))

    ! Average C55 on (1/2, 1/2) node
#define c55_eff_xz (4.0/sum(1.0/c55(i:i + 1, j:j + 1)))

contains

    !
    !> Update wavefields
    !
    subroutine update_wavefield(dt, &
            stressxx, stresszz, stressxz, vx, vz, &
            memory_pdxvx, memory_pdzvx, memory_pdxvz, memory_pdzvz, &
            memory_pdxxx, memory_pdzxz, memory_pdxxz, memory_pdzzz)

        real, intent(in) :: dt
        real, allocatable, dimension(:, :), intent(inout) :: &
            stressxx, stresszz, stressxz, vx, vz, &
            memory_pdxvx, memory_pdzvx, memory_pdxvz, memory_pdzvz, &
            memory_pdxxx, memory_pdzxz, memory_pdxxz, memory_pdzzz

        integer :: i, j
        real :: pdxvx, pdzvx, pdxvz, pdzvz
        real :: pdxxx, pdzxz, pdxxz, pdzzz

        !$omp parallel do private(i, j, pdxvx, pdzvx, pdxvz, pdzvz) collapse(2) schedule(auto)
        do j = -pml + 1, nz + pml
            do i = -pml + 1, nx + pml

                pdxvx = idx*pdxvx_stencil
                pdzvz = idz*pdzvz_stencil

                if (i <= 0 .or. i >= nx + 1 .or. j <= 0 .or. j >= nz + 1) then
                    memory_pdxvx(i, j) = axii(i, j)*pdxvx + bxii(i, j)*memory_pdxvx(i, j)
                    pdxvx = (pdxvx + memory_pdxvx(i, j))/kxii(i, j)
                    memory_pdzvz(i, j) = azii(i, j)*pdzvz + bzii(i, j)*memory_pdzvz(i, j)
                    pdzvz = (pdzvz + memory_pdzvz(i, j))/kzii(i, j)
                end if

                stressxx(i, j) = stressxx(i, j) + dt*(c11(i, j)*pdxvx + c13(i, j)*pdzvz)
                stresszz(i, j) = stresszz(i, j) + dt*(c13(i, j)*pdxvx + c33(i, j)*pdzvz)

            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, pdxvx, pdzvx, pdxvz, pdzvz) collapse(2) schedule(auto)
        do j = -pml, nz + pml - 1
            do i = -pml, nx + pml - 1

                pdzvx = idz*pdzvx_stencil
                pdxvz = idx*pdxvz_stencil

                if (i + 1 <= 1 .or. i + 1 >= nx + 1 .or. j + 1 <= 1 .or. j + 1 >= nz + 1) then
                    memory_pdxvz(i + 1, j + 1) = axhh(i + 1, j + 1)*pdxvz + bxhh(i + 1, j + 1)*memory_pdxvz(i + 1, j + 1)
                    pdxvz = (pdxvz + memory_pdxvz(i + 1, j + 1))/kxhh(i + 1, j + 1)
                    memory_pdzvx(i + 1, j + 1) = azhh(i + 1, j + 1)*pdzvx + bzhh(i + 1, j + 1)*memory_pdzvx(i + 1, j + 1)
                    pdzvx = (pdzvx + memory_pdzvx(i + 1, j + 1))/kzhh(i + 1, j + 1)
                end if

                stressxz(i + 1, j + 1) = stressxz(i + 1, j + 1) + dt*c55_eff_xz*(pdzvx + pdxvz)

            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, pdxxx, pdzxz, pdxxz, pdzzz) collapse(2) schedule(auto)
        do j = -pml + 1, nz + pml
            do i = -pml, nx + pml - 1

                pdxxx = idx*pdxxx_stencil
                pdzxz = idz*pdzxz_stencil

                if (i + 1 <= 1 .or. i + 1 >= nx + 1 .or. j <= 0 .or. j >= nz + 1) then
                    memory_pdxxx(i + 1, j) = axhi(i + 1, j)*pdxxx + bxhi(i + 1, j)*memory_pdxxx(i + 1, j)
                    pdxxx = (pdxxx + memory_pdxxx(i + 1, j))/kxhi(i + 1, j)
                    memory_pdzxz(i + 1, j) = azhi(i + 1, j)*pdzxz + bzhi(i + 1, j)*memory_pdzxz(i + 1, j)
                    pdzxz = (pdzxz + memory_pdzxz(i + 1, j))/kzhi(i + 1, j)
                end if

                vx(i + 1, j) = vx(i + 1, j) + dt/rho_eff_x*(pdxxx + pdzxz)

            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, pdxxx, pdzxz, pdxxz, pdzzz) collapse(2) schedule(auto)
        do j = -pml, nz + pml - 1
            do i = -pml + 1, nx + pml

                pdxxz = idx*pdxxz_stencil
                pdzzz = idz*pdzzz_stencil

                if (i <= 0 .or. i >= nx + 1 .or. j + 1 <= 1 .or. j + 1 >= nz + 1) then
                    memory_pdxxz(i, j + 1) = axih(i, j + 1)*pdxxz + bxih(i, j + 1)*memory_pdxxz(i, j + 1)
                    pdxxz = (pdxxz + memory_pdxxz(i, j + 1))/kxih(i, j + 1)
                    memory_pdzzz(i, j + 1) = azih(i, j + 1)*pdzzz + bzih(i, j + 1)*memory_pdzzz(i, j + 1)
                    pdzzz = (pdzzz + memory_pdzzz(i, j + 1))/kzih(i, j + 1)
                end if

                vz(i, j + 1) = vz(i, j + 1) + dt/rho_eff_z*(pdxxz + pdzzz)

            end do
        end do
        !$omp end parallel do

    end subroutine update_wavefield

    !
    !> Update wavefields with planar free surface
    !
    subroutine update_wavefield_free_surface(dt, &
            stressxx, stresszz, stressxz, vx, vz, &
            memory_pdxvx, memory_pdzvx, memory_pdxvz, memory_pdzvz, &
            memory_pdxxx, memory_pdzxz, memory_pdxxz, memory_pdzzz)

        real, intent(in) :: dt
        real, allocatable, dimension(:, :), intent(inout) :: &
            stressxx, stresszz, stressxz, vx, vz, &
            memory_pdxvx, memory_pdzvx, memory_pdxvz, memory_pdzvz, &
            memory_pdxxx, memory_pdzxz, memory_pdxxz, memory_pdzzz

        integer :: i, j
        real :: pdxvx, pdzvx, pdxvz, pdzvz
        real :: pdxxx, pdzxz, pdxxz, pdzzz

        ! Mirror the particle velocities above the free surface, v(-z) = v(z). The z-derivative
        ! stencils of the stress update read up to fdhalf rows above the surface; zeros there make
        ! the velocities discontinuous at the surface, an error that does not vanish with grid
        ! refinement (surface amplitudes 6-8% too low), while the mirrored values keep them continuous.
        !$omp parallel do private(j) schedule(auto)
        do j = 1, fdhalf
            vx(:, 1 - j) = vx(:, 1 + j)
            vz(:, 2 - j) = vz(:, 1 + j)
        end do
        !$omp end parallel do

        ! Update stress xx, zz
        !$omp parallel do private(i, j, pdxvx, pdzvx, pdxvz, pdzvz) collapse(2) schedule(auto)
        do j = 1, nz + pml
            do i = -pml + 1, nx + pml

                pdxvx = idx*pdxvx_stencil

                if (j == 1) then
                    ! At the free surface, pdzvz is computed from pdxvx to satisfy free-surface condition

                    if (i <= 0 .or. i >= nx + 1) then
                        memory_pdxvx(i, j) = axii(i, j)*pdxvx + bxii(i, j)*memory_pdxvx(i, j)
                        pdxvx = (pdxvx + memory_pdxvx(i, j))/kxii(i, j)
                    end if

                    pdzvz = -c13(i, j)/c33(i, j)*pdxvx

                else
                    ! Below the free surface, pdzvz is computed normally
                    pdzvz = idz*pdzvz_stencil*dz_scaling_i(j)

                    if (i <= 0 .or. i >= nx + 1 .or. j >= nz + 1) then
                        memory_pdxvx(i, j) = axii(i, j)*pdxvx + bxii(i, j)*memory_pdxvx(i, j)
                        pdxvx = (pdxvx + memory_pdxvx(i, j))/kxii(i, j)
                        memory_pdzvz(i, j) = azii(i, j)*pdzvz + bzii(i, j)*memory_pdzvz(i, j)
                        pdzvz = (pdzvz + memory_pdzvz(i, j))/kzii(i, j)
                    end if

                end if

                stressxx(i, j) = stressxx(i, j) + dt*(c11(i, j)*pdxvx + c13(i, j)*pdzvz)
                stresszz(i, j) = stresszz(i, j) + dt*(c13(i, j)*pdxvx + c33(i, j)*pdzvz)

            end do
        end do
        !$omp end parallel do

        ! Update stress xz
        !$omp parallel do private(i, j, pdxvx, pdzvx, pdxvz, pdzvz) collapse(2) schedule(auto)
        do j = 1, nz + pml - 1
            do i = -pml, nx + pml - 1

                pdxvz = idx*pdxvz_stencil
                pdzvx = idz*pdzvx_stencil*dz_scaling_h(j + 1)

                if (i + 1 <= 1 .or. i + 1 >= nx + 1 .or. j + 1 >= nz + 1) then
                    memory_pdxvz(i + 1, j + 1) = axhh(i + 1, j + 1)*pdxvz + bxhh(i + 1, j + 1)*memory_pdxvz(i + 1, j + 1)
                    pdxvz = (pdxvz + memory_pdxvz(i + 1, j + 1))/kxhh(i + 1, j + 1)
                    memory_pdzvx(i + 1, j + 1) = azhh(i + 1, j + 1)*pdzvx + bzhh(i + 1, j + 1)*memory_pdzvx(i + 1, j + 1)
                    pdzvx = (pdzvx + memory_pdzvx(i + 1, j + 1))/kzhh(i + 1, j + 1)
                end if

                stressxz(i + 1, j + 1) = stressxz(i + 1, j + 1) + dt*c55_eff_xz*(pdzvx + pdxvz)

            end do
        end do
        !$omp end parallel do

        ! Stress mirror
        !$omp parallel do private(j) schedule(auto)
        do j = 1, fdhalf

            ! The normal stress is mirror-symmetric w.r.t. the free surface, i.e., j = 1
            stresszz(:, 1 - j) = -stresszz(:, 1 + j)

            ! At the free surface, normal stress is strictly zero
            stresszz(:, 1) = 0.0

            ! The shear stress is mirror-symmetric w.r.t. the free surface
            stressxz(:, 2 - j) = -stressxz(:, 1 + j)

        end do
        !$omp end parallel do

        ! Vx
        !$omp parallel do private(i, j, pdxxx, pdzxz, pdxxz, pdzzz) collapse(2) schedule(auto)
        do j = 1, nz + pml
            do i = -pml, nx + pml - 1

                pdxxx = idx*pdxxx_stencil
                pdzxz = idz*pdzxz_stencil*dz_scaling_i(j)

                if (i + 1 <= 1 .or. i + 1 >= nx + 1 .or. j >= nz + 1) then
                    memory_pdxxx(i + 1, j) = axhi(i + 1, j)*pdxxx + bxhi(i + 1, j)*memory_pdxxx(i + 1, j)
                    pdxxx = (pdxxx + memory_pdxxx(i + 1, j))/kxhi(i + 1, j)
                    memory_pdzxz(i + 1, j) = azhi(i + 1, j)*pdzxz + bzhi(i + 1, j)*memory_pdzxz(i + 1, j)
                    pdzxz = (pdzxz + memory_pdzxz(i + 1, j))/kzhi(i + 1, j)
                end if

                vx(i + 1, j) = vx(i + 1, j) + dt/rho_eff_x*(pdxxx + pdzxz)

            end do
        end do
        !$omp end parallel do

        ! Vz
        !$omp parallel do private(i, j, pdxxx, pdzxz, pdxxz, pdzzz) collapse(2) schedule(auto)
        do j = 1, nz + pml - 1
            do i = -pml + 1, nx + pml

                pdxxz = idx*pdxxz_stencil
                pdzzz = idz*pdzzz_stencil*dz_scaling_h(j + 1)

                if (i <= 0 .or. i >= nx + 1 .or. j + 1 >= nz + 1) then
                    memory_pdxxz(i, j + 1) = axih(i, j + 1)*pdxxz + bxih(i, j + 1)*memory_pdxxz(i, j + 1)
                    pdxxz = (pdxxz + memory_pdxxz(i, j + 1))/kxih(i, j + 1)
                    memory_pdzzz(i, j + 1) = azih(i, j + 1)*pdzzz + bzih(i, j + 1)*memory_pdzzz(i, j + 1)
                    pdzzz = (pdzzz + memory_pdzzz(i, j + 1))/kzih(i, j + 1)
                end if

                vz(i, j + 1) = vz(i, j + 1) + dt/rho_eff_z*(pdxxz + pdzzz)

            end do
        end do
        !$omp end parallel do

    end subroutine

    !
    !> Add source
    !
    !> With a free surface, every source comes with its mirror image (see add_source_value_2d). The
    !> image is 2P - S, where S is the source and P the part of it that a traction-free surface
    !> admits, so that on the surface row, where the source and its image coincide, they add up to
    !> 2P on that half cell. A force is fully admitted (image = force). A moment tensor loses its
    !> traction components Mzz and Mxz: sigma_zz = 0 at the surface gives
    !> eps_zz = (Mzz - c13*eps_xx)/c33, which turns Mzz into the horizontal moment -(c13/c33)*Mzz,
    !> so P = (Mxx - c13/c33*Mzz, 0, 0) and the images are Mxx - 2*c13/c33*Mzz, -Mzz and -Mxz.
    !
    subroutine add_source(t)

        integer, intent(in) :: t

        integer :: k
        integer :: nbeg, nend
        integer :: sgx, sgz, shx, shz
        real :: amp, polar, mxx, mzz, mxz, r13

        do k = 1, sgmtr%ns

            nbeg = nint(sgmtr%srcr(k)%t0/dt) + 1
            nend = nbeg + sgmtr%srcr(k)%nt - 1

            sgx = sgmtr%srcr(k)%gx
            sgz = sgmtr%srcr(k)%gz
            shx = sgmtr%srcr(k)%hx
            shz = sgmtr%srcr(k)%hz

            if (t >= nbeg .and. t <= nend) then

                amp = sgmtr%srcr(k)%stf(t - nbeg + 1)*sgmtr%srcr(k)%amp*dt

                associate (s => sgmtr%srcr(k))

                select case (s%mechanism)

                    case ('force')
                        ! Force vector, added as a force density on a cell of the nominal size dx*dz,
                        ! scaled on the refined mesh to the height of each cell (dz_scaling), so that
                        ! the total force is stf*dx*dz at any depth
                        polar = s%polar
                        call add_source_value_2d(vx, shx, sgz, s%interp_hx, s%interp_iz, &
                            sin(polar)*amp/rho(shx, sgz), sin(polar)*amp/rho(shx, sgz), .false., yn_free_surface, &
                            dz_scaling_i)
                        call add_source_value_2d(vz, sgx, shz, s%interp_ix, s%interp_hz, &
                            cos(polar)*amp/rho(sgx, shz), cos(polar)*amp/rho(sgx, shz), .true., yn_free_surface, &
                            dz_scaling_h)

                    case ('explosion', 'mt')
                        ! Moment tensor (explosion: identity) added as stress drop
                        if (s%mechanism == 'explosion') then
                            mxx = 1.0
                            mzz = 1.0
                            mxz = 0.0
                        else
                            mxx = s%moment_tensor(1, 1)
                            mzz = s%moment_tensor(3, 3)
                            mxz = s%moment_tensor(1, 3)
                        end if
                        ! Moment density on a cell of the nominal size dx*dz, scaled on the refined mesh
                        ! to the height of each cell (dz_scaling)
                        amp = amp/(dx*dz)
                        ! Free-surface conversion ratio c13/c33 at the surface above the source
                        r13 = c13(sgx, 1)/c33(sgx, 1)
                        call add_source_value_2d(stressxx, sgx, sgz, s%interp_ix, s%interp_iz, &
                            -amp*mxx, -amp*(mxx - 2.0*r13*mzz), .false., yn_free_surface, dz_scaling_i)
                        call add_source_value_2d(stresszz, sgx, sgz, s%interp_ix, s%interp_iz, &
                            -amp*mzz, amp*mzz, .false., yn_free_surface, dz_scaling_i)
                        call add_source_value_2d(stressxz, shx, shz, s%interp_hx, s%interp_hz, &
                            -amp*mxz, amp*mxz, .true., yn_free_surface, dz_scaling_h)

                end select

                end associate

            end if
        end do

    end subroutine

end module
