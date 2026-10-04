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

module elastic_vhtiort_3d_wavefield

    use libflit
    use elastic_vhtiort_3d_vars
    use elastic_vhtiort_3d_cfspml

    implicit none

    ! FD stencil
#include 'macro_fd_stencil.f90'

    ! Average rho on half_x, half_y or half_z node
#define rho_eff_x (0.5*sum(rho(i:i + 1, j, k)))
#define rho_eff_y (0.5*sum(rho(i, j:j + 1, k)))
#define rho_eff_z (0.5*sum(rho(i, j, k:k + 1)))

    ! Average mu on half-half-integer, integer-half-half, or half-integer-half nodes
#define c44_eff_yz (4.0/sum(1.0/c44(i, j:j + 1, k:k + 1)))
#define c55_eff_xz (4.0/sum(1.0/c55(i:i + 1, j, k:k + 1)))
#define c66_eff_xy (4.0/sum(1.0/c66(i:i + 1, j:j + 1, k)))

    ! Macros for MPML wavefield damping
#define damp_coef(dt, a_, b_, a, b) \
    a = -2.0*abs(dt)*a_/(2.0 + abs(dt)*b_); \
    b = (2.0 - abs(dt)*b_)/(2.0 + abs(dt)*b_);
#define damp_wavefield(w, v, da, db, dk, i, j, k) \
    w(i, j, k) = da*v + db*w(i, j, k); \
    v = (v + w(i, j, k))/dk;
contains

    !
    !> Update wavefields
    !
    subroutine update_wavefield(dt, &
            stressxx, stressyy, stresszz, &
            stressyz, stressxz, stressxy, &
            vx, vy, vz, &
            memory_pdxxx, memory_pdyxy, memory_pdzxz, &
            memory_pdxxy, memory_pdyyy, memory_pdzyz, &
            memory_pdxxz, memory_pdyyz, memory_pdzzz, &
            memory_pdxvx, memory_pdyvx, memory_pdzvx, &
            memory_pdxvy, memory_pdyvy, memory_pdzvy, &
            memory_pdxvz, memory_pdyvz, memory_pdzvz)

        real, intent(in) :: dt
        real, allocatable, dimension(:, :, :), intent(inout) :: stressxx, stressyy, stresszz
        real, allocatable, dimension(:, :, :), intent(inout) :: stressyz, stressxz, stressxy
        real, allocatable, dimension(:, :, :), intent(inout) :: vx, vy, vz
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxxx, memory_pdyxy, memory_pdzxz
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxxy, memory_pdyyy, memory_pdzyz
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxxz, memory_pdyyz, memory_pdzzz
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxvx, memory_pdyvx, memory_pdzvx
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxvy, memory_pdyvy, memory_pdzvy
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxvz, memory_pdyvz, memory_pdzvz

        integer :: i, j, k
        real :: pdxvx, pdyvx, pdzvx
        real :: pdxvy, pdyvy, pdzvy
        real :: pdxvz, pdyvz, pdzvz
        real :: pdxxx, pdyxy, pdzxz
        real :: pdxxy, pdyyy, pdzyz
        real :: pdxxz, pdyyz, pdzzz
        real :: ratiox, ratioy, ratioz
        real :: dax, day, daz
        real :: a, b
        real :: ax, ay, az, bx, by, bz, kx, ky, kz

        call commute_array_group(vx, fdhalf)
        call commute_array_group(vy, fdhalf)
        call commute_array_group(vz, fdhalf)

        !$omp parallel do private(i, j, k, pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1, nx2

                    pdxvx = idx*pdxvx_stencil
                    pdyvy = idy*pdyvy_stencil
                    pdzvz = idz*pdzvz_stencil

                    ! MPML
                    if (i <= 0 .or. i >= nx + 1 &
                            .or. j <= 0 .or. j >= ny + 1 &
                            .or. k <= 0 .or. k >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i <= 0) then
                            ratiox = dampratio_left(j, k)
                        else if (i >= nx + 1) then
                            ratiox = dampratio_right(j, k)
                        end if
                        if (j <= 0) then
                            ratioy = dampratio_front(i, k)
                        else if (j >= ny + 1) then
                            ratioy = dampratio_back(i, k)
                        end if
                        if (k <= 0) then
                            ratioz = dampratio_top(i, j)
                        else if (k >= nz + 1) then
                            ratioz = dampratio_bottom(i, j)
                        end if

                        dax = idax(i, j, k) + ratioy*iday(i, j, k) + ratioz*idaz(i, j, k)
                        a = dax/kappaxi(i, j, k)
                        b = dax/kappaxi(i, j, k) + alphaxi(i, j, k)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxi(i, j, k)
                        damp_wavefield(memory_pdxvx, pdxvx, ax, bx, kx, i, j, k)

                        day = iday(i, j, k) + ratiox*idax(i, j, k) + ratioz*idaz(i, j, k)
                        a = day/kappayi(i, j, k)
                        b = day/kappayi(i, j, k) + alphayi(i, j, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i, j, k)
                        damp_wavefield(memory_pdyvy, pdyvy, ay, by, ky, i, j, k)

                        daz = idaz(i, j, k) + ratiox*idax(i, j, k) + ratioy*iday(i, j, k)
                        a = daz/kappazi(i, j, k)
                        b = daz/kappazi(i, j, k) + alphazi(i, j, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i, j, k)
                        damp_wavefield(memory_pdzvz, pdzvz, az, bz, kz, i, j, k)

                    end if

                    stressxx(i, j, k) = stressxx(i, j, k) + dt*( &
                        c11(i, j, k)*pdxvx &
                        + c12(i, j, k)*pdyvy &
                        + c13(i, j, k)*pdzvz)

                    stressyy(i, j, k) = stressyy(i, j, k) + dt*( &
                        c12(i, j, k)*pdxvx &
                        + c22(i, j, k)*pdyvy &
                        + c23(i, j, k)*pdzvz)

                    stresszz(i, j, k) = stresszz(i, j, k) + dt*( &
                        c13(i, j, k)*pdxvx &
                        + c23(i, j, k)*pdyvy &
                        + c33(i, j, k)*pdzvz)

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1 - 1, nz2 - 1
            do j = ny1 - 1, ny2 - 1
                do i = nx1, nx2

                    pdzvy = idz*pdzvy_stencil
                    pdyvz = idy*pdyvz_stencil

                    ! MPML
                    if (i <= 0 .or. i >= nx + 1 &
                            .or. j + 1 <= 1 .or. j + 1 >= ny + 1 &
                            .or. k + 1 <= 1 .or. k + 1 >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i <= 0) then
                            ratiox = dampratio_left(j + 1, k + 1)
                        else if (i >= nx + 1) then
                            ratiox = dampratio_right(j + 1, k + 1)
                        end if
                        if (j + 1 <= 1) then
                            ratioy = dampratio_front(i, k + 1)
                        else if (j + 1 >= ny + 1) then
                            ratioy = dampratio_back(i, k + 1)
                        end if
                        if (k + 1 <= 1) then
                            ratioz = dampratio_top(i, j + 1)
                        else if (k + 1 >= nz + 1) then
                            ratioz = dampratio_bottom(i, j + 1)
                        end if

                        day = hday(i, j + 1, k + 1) + ratiox*idax(i, j + 1, k + 1) + ratioz*hdaz(i, j + 1, k + 1)
                        a = day/kappayh(i, j + 1, k + 1)
                        b = day/kappayh(i, j + 1, k + 1) + alphayh(i, j + 1, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i, j + 1, k + 1)
                        damp_wavefield(memory_pdyvz, pdyvz, ay, by, ky, i, j + 1, k + 1)

                        daz = hdaz(i, j + 1, k + 1) + ratiox*idax(i, j + 1, k + 1) + ratioy*hday(i, j + 1, k + 1)
                        a = daz/kappazh(i, j + 1, k + 1)
                        b = daz/kappazh(i, j + 1, k + 1) + alphazh(i, j + 1, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i, j + 1, k + 1)
                        damp_wavefield(memory_pdzvy, pdzvy, az, bz, kz, i, j + 1, k + 1)

                    end if

                    stressyz(i, j + 1, k + 1) = stressyz(i, j + 1, k + 1) + dt*c44_eff_yz*(pdzvy + pdyvz)

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1 - 1, nz2 - 1
            do j = ny1, ny2
                do i = nx1 - 1, nx2 - 1

                    pdzvx = idz*pdzvx_stencil
                    pdxvz = idx*pdxvz_stencil

                    ! MPML
                    if (i + 1 <= 1 .or. i + 1 >= nx + 1 &
                            .or. j <= 0 .or. j >= ny + 1 &
                            .or. k + 1 <= 1 .or. k + 1 >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i + 1 <= 1) then
                            ratiox = dampratio_left(j, k + 1)
                        else if (i + 1 >= nx + 1) then
                            ratiox = dampratio_right(j, k + 1)
                        end if
                        if (j <= 0) then
                            ratioy = dampratio_front(i + 1, k + 1)
                        else if (j >= ny + 1) then
                            ratioy = dampratio_back(i + 1, k + 1)
                        end if
                        if (k + 1 <= 1) then
                            ratioz = dampratio_top(i + 1, j)
                        else if (k + 1 >= nz + 1) then
                            ratioz = dampratio_bottom(i + 1, j)
                        end if

                        dax = hdax(i + 1, j, k + 1) + ratioy*iday(i + 1, j, k + 1) + ratioz*hdaz(i + 1, j, k + 1)
                        a = dax/kappaxh(i + 1, j, k + 1)
                        b = dax/kappaxh(i + 1, j, k + 1) + alphaxh(i + 1, j, k + 1)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxh(i + 1, j, k + 1)
                        damp_wavefield(memory_pdxvz, pdxvz, ax, bx, kx, i + 1, j, k + 1)

                        daz = hdaz(i + 1, j, k + 1) + ratiox*hdax(i + 1, j, k + 1) + ratioy*iday(i + 1, j, k + 1)
                        a = daz/kappazh(i + 1, j, k + 1)
                        b = daz/kappazh(i + 1, j, k + 1) + alphazh(i + 1, j, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i + 1, j, k + 1)
                        damp_wavefield(memory_pdzvx, pdzvx, az, bz, kz, i + 1, j, k + 1)

                    end if

                    stressxz(i + 1, j, k + 1) = stressxz(i + 1, j, k + 1) + dt*c55_eff_xz*(pdzvx + pdxvz)

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1 - 1, ny2 - 1
                do i = nx1 - 1, nx2 - 1

                    pdyvx = idy*pdyvx_stencil
                    pdxvy = idx*pdxvy_stencil

                    ! MPML
                    if (i + 1 <= 1 .or. i + 1 >= nx + 1 &
                            .or. j + 1 <= 1 .or. j + 1 >= ny + 1 &
                            .or. k <= 0 .or. k >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i + 1 <= 1) then
                            ratiox = dampratio_left(j + 1, k)
                        else if (i + 1 >= nx + 1) then
                            ratiox = dampratio_right(j + 1, k)
                        end if
                        if (j + 1 <= 1) then
                            ratioy = dampratio_front(i + 1, k)
                        else if (j + 1 >= ny + 1) then
                            ratioy = dampratio_back(i + 1, k)
                        end if
                        if (k <= 0) then
                            ratioz = dampratio_top(i + 1, j + 1)
                        else if (k >= nz + 1) then
                            ratioz = dampratio_bottom(i + 1, j + 1)
                        end if

                        dax = hdax(i + 1, j + 1, k) + ratioy*hday(i + 1, j + 1, k) + ratioz*idaz(i + 1, j + 1, k)
                        a = dax/kappaxh(i + 1, j + 1, k)
                        b = dax/kappaxh(i + 1, j + 1, k) + alphaxh(i + 1, j + 1, k)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxh(i + 1, j + 1, k)
                        damp_wavefield(memory_pdxvy, pdxvy, ax, bx, kx, i + 1, j + 1, k)

                        day = hday(i + 1, j + 1, k) + ratiox*hdax(i + 1, j + 1, k) + ratioz*idaz(i + 1, j + 1, k)
                        a = day/kappayh(i + 1, j + 1, k)
                        b = day/kappayh(i + 1, j + 1, k) + alphayh(i + 1, j + 1, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i + 1, j + 1, k)
                        damp_wavefield(memory_pdyvx, pdyvx, ay, by, ky, i + 1, j + 1, k)

                    end if

                    stressxy(i + 1, j + 1, k) = stressxy(i + 1, j + 1, k) + dt*c66_eff_xy*(pdyvx + pdxvy)

                end do
            end do
        end do
        !$omp end parallel do

        call commute_array_group(stressxx, fdhalf)
        call commute_array_group(stressyy, fdhalf)
        call commute_array_group(stresszz, fdhalf)
        call commute_array_group(stressyz, fdhalf)
        call commute_array_group(stressxz, fdhalf)
        call commute_array_group(stressxy, fdhalf)

        ! vx component
        !$omp parallel do private(i, j, k, pdxxx, pdyxy, pdzxz, pdxxy, pdyyy, pdzyz, pdxxz, pdyyz, pdzzz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1 - 1, nx2 - 1

                    pdxxx = idx*pdxxx_stencil
                    pdyxy = idy*pdyxy_stencil
                    pdzxz = idz*pdzxz_stencil

                    ! MPML
                    if (i + 1 <= 1 .or. i + 1 >= nx + 1 &
                            .or. j <= 0 .or. j >= ny + 1 &
                            .or. k <= 0 .or. k >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i + 1 <= 1) then
                            ratiox = dampratio_left(j, k)
                        else if (i + 1 >= nx + 1) then
                            ratiox = dampratio_right(j, k)
                        end if
                        if (j <= 0) then
                            ratioy = dampratio_front(i + 1, k)
                        else if (j >= ny + 1) then
                            ratioy = dampratio_back(i + 1, k)
                        end if
                        if (k <= 0) then
                            ratioz = dampratio_top(i + 1, j)
                        else if (k >= nz + 1) then
                            ratioz = dampratio_bottom(i + 1, j)
                        end if

                        dax = hdax(i + 1, j, k) + ratioy*iday(i + 1, j, k) + ratioz*idaz(i + 1, j, k)
                        a = dax/kappaxh(i + 1, j, k)
                        b = dax/kappaxh(i + 1, j, k) + alphaxh(i + 1, j, k)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxh(i + 1, j, k)
                        damp_wavefield(memory_pdxxx, pdxxx, ax, bx, kx, i + 1, j, k)

                        day = iday(i + 1, j, k) + ratiox*hdax(i + 1, j, k) + ratioz*idaz(i + 1, j, k)
                        a = day/kappayi(i + 1, j, k)
                        b = day/kappayi(i + 1, j, k) + alphayi(i + 1, j, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i + 1, j, k)
                        damp_wavefield(memory_pdyxy, pdyxy, ay, by, ky, i + 1, j, k)

                        daz = idaz(i + 1, j, k) + ratiox*hdax(i + 1, j, k) + ratioy*iday(i + 1, j, k)
                        a = daz/kappazi(i + 1, j, k)
                        b = daz/kappazi(i + 1, j, k) + alphazi(i + 1, j, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i + 1, j, k)
                        damp_wavefield(memory_pdzxz, pdzxz, az, bz, kz, i + 1, j, k)

                    end if

                    vx(i + 1, j, k) = vx(i + 1, j, k) + dt/rho_eff_x*(pdxxx + pdyxy + pdzxz)

                end do
            end do
        end do
        !$omp end parallel do

        ! vy component
        !$omp parallel do private(i, j, k, pdxxx, pdyxy, pdzxz, pdxxy, pdyyy, pdzyz, pdxxz, pdyyz, pdzzz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1 - 1, ny2 - 1
                do i = nx1, nx2

                    pdxxy = idx*pdxxy_stencil
                    pdyyy = idy*pdyyy_stencil
                    pdzyz = idz*pdzyz_stencil

                    ! MPML
                    if (i <= 0 .or. i >= nx + 1 &
                            .or. j + 1 <= 1 .or. j + 1 >= ny + 1 &
                            .or. k <= 0 .or. k >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i <= 0) then
                            ratiox = dampratio_left(j + 1, k)
                        else if (i >= nx + 1) then
                            ratiox = dampratio_right(j + 1, k)
                        end if
                        if (j + 1 <= 1) then
                            ratioy = dampratio_front(i, k)
                        else if (j + 1 >= ny + 1) then
                            ratioy = dampratio_back(i, k)
                        end if
                        if (k <= 0) then
                            ratioz = dampratio_top(i, j + 1)
                        else if (k >= nz + 1) then
                            ratioz = dampratio_bottom(i, j + 1)
                        end if

                        dax = idax(i, j + 1, k) + ratioy*hday(i, j + 1, k) + ratioz*idaz(i, j + 1, k)
                        a = dax/kappaxi(i, j + 1, k)
                        b = dax/kappaxi(i, j + 1, k) + alphaxi(i, j + 1, k)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxi(i, j + 1, k)
                        damp_wavefield(memory_pdxxy, pdxxy, ax, bx, kx, i, j + 1, k)

                        day = hday(i, j + 1, k) + ratiox*idax(i, j + 1, k) + ratioz*idaz(i, j + 1, k)
                        a = day/kappayh(i, j + 1, k)
                        b = day/kappayh(i, j + 1, k) + alphayh(i, j + 1, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i, j + 1, k)
                        damp_wavefield(memory_pdyyy, pdyyy, ay, by, ky, i, j + 1, k)

                        daz = idaz(i, j + 1, k) + ratiox*idax(i, j + 1, k) + ratioy*hday(i, j + 1, k)
                        a = daz/kappazi(i, j + 1, k)
                        b = daz/kappazi(i, j + 1, k) + alphazi(i, j + 1, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i, j + 1, k)
                        damp_wavefield(memory_pdzyz, pdzyz, az, bz, kz, i, j + 1, k)

                    end if

                    vy(i, j + 1, k) = vy(i, j + 1, k) + dt/rho_eff_y*(pdxxy + pdyyy + pdzyz)

                end do
            end do
        end do
        !$omp end parallel do

        ! vz component
        !$omp parallel do private(i, j, k, pdxxx, pdyxy, pdzxz, pdxxy, pdyyy, pdzyz, pdxxz, pdyyz, pdzzz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1 - 1, nz2 - 1
            do j = ny1, ny2
                do i = nx1, nx2

                    pdxxz = idx*pdxxz_stencil
                    pdyyz = idy*pdyyz_stencil
                    pdzzz = idz*pdzzz_stencil

                    ! MPML
                    if (i <= 0 .or. i >= nx + 1 &
                            .or. j <= 0 .or. j >= ny + 1 &
                            .or. k + 1 <= 1 .or. k + 1 >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i <= 0) then
                            ratiox = dampratio_left(j, k + 1)
                        else if (i >= nx + 1) then
                            ratiox = dampratio_right(j, k + 1)
                        end if
                        if (j <= 0) then
                            ratioy = dampratio_front(i, k + 1)
                        else if (j >= ny + 1) then
                            ratioy = dampratio_back(i, k + 1)
                        end if
                        if (k + 1 <= 1) then
                            ratioz = dampratio_top(i, j)
                        else if (k + 1 >= nz + 1) then
                            ratioz = dampratio_bottom(i, j)
                        end if

                        dax = idax(i, j, k + 1) + ratioy*iday(i, j, k + 1) + ratioz*hdaz(i, j, k + 1)
                        a = dax/kappaxi(i, j, k + 1)
                        b = dax/kappaxi(i, j, k + 1) + alphaxi(i, j, k + 1)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxi(i, j, k + 1)
                        damp_wavefield(memory_pdxxz, pdxxz, ax, bx, kx, i, j, k + 1)

                        day = iday(i, j, k + 1) + ratiox*idax(i, j, k + 1) + ratioz*hdaz(i, j, k + 1)
                        a = day/kappayi(i, j, k + 1)
                        b = day/kappayi(i, j, k + 1) + alphayi(i, j, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i, j, k + 1)
                        damp_wavefield(memory_pdyyz, pdyyz, ay, by, ky, i, j, k + 1)

                        daz = hdaz(i, j, k + 1) + ratiox*idax(i, j, k + 1) + ratioy*iday(i, j, k + 1)
                        a = daz/kappazh(i, j, k + 1)
                        b = daz/kappazh(i, j, k + 1) + alphazh(i, j, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i, j, k + 1)
                        damp_wavefield(memory_pdzzz, pdzzz, az, bz, kz, i, j, k + 1)

                    end if

                    vz(i, j, k + 1) = vz(i, j, k + 1) + dt/rho_eff_z*(pdxxz + pdyyz + pdzzz)

                end do
            end do
        end do
        !$omp end parallel do

    end subroutine

    !
    !> Update wavefields in a 3D elastic medium with free-surface condition
    !
    subroutine update_wavefield_free_surface(dt, &
            stressxx, stressyy, stresszz, &
            stressyz, stressxz, stressxy, &
            vx, vy, vz, &
            memory_pdxxx, memory_pdyxy, memory_pdzxz, &
            memory_pdxxy, memory_pdyyy, memory_pdzyz, &
            memory_pdxxz, memory_pdyyz, memory_pdzzz, &
            memory_pdxvx, memory_pdyvx, memory_pdzvx, &
            memory_pdxvy, memory_pdyvy, memory_pdzvy, &
            memory_pdxvz, memory_pdyvz, memory_pdzvz)

        real, intent(in) :: dt
        real, allocatable, dimension(:, :, :), intent(inout) :: stressxx, stressyy, stresszz
        real, allocatable, dimension(:, :, :), intent(inout) :: stressyz, stressxz, stressxy
        real, allocatable, dimension(:, :, :), intent(inout) :: vx, vy, vz
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxxx, memory_pdyxy, memory_pdzxz
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxxy, memory_pdyyy, memory_pdzyz
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxxz, memory_pdyyz, memory_pdzzz
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxvx, memory_pdyvx, memory_pdzvx
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxvy, memory_pdyvy, memory_pdzvy
        real, allocatable, dimension(:, :, :), intent(inout) :: memory_pdxvz, memory_pdyvz, memory_pdzvz

        integer :: i, j, k
        real :: pdxvx, pdyvx, pdzvx
        real :: pdxvy, pdyvy, pdzvy
        real :: pdxvz, pdyvz, pdzvz
        real :: pdxxx, pdyxy, pdzxz
        real :: pdxxy, pdyyy, pdzyz
        real :: pdxxz, pdyyz, pdzzz
        real :: ratiox, ratioy, ratioz
        real :: dax, day, daz
        real :: a, b, ax, ay, az, bx, by, bz, kx, ky, kz

        ! Particle velocities above free surface are zero
        !$omp parallel do private(k) schedule(auto)
        do k = 1, fdhalf
            if (1 - k >= nz1 .and. 1 - k <= nz2) then
                vx(:, :, 1 - k) = 0.0
                vy(:, :, 1 - k) = 0.0
            end if
            if (2 - k >= nz1 .and. 2 - k <= nz2) then
                vz(:, :, 2 - k) = 0.0
            end if
        end do
        !$omp end parallel do

        call commute_array_group(vx, fdhalf)
        call commute_array_group(vy, fdhalf)
        call commute_array_group(vz, fdhalf)

        !$omp parallel do private(i, j, k, pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1), nz2
            do j = ny1, ny2
                do i = nx1, nx2

                    pdxvx = idx*pdxvx_stencil
                    pdyvy = idy*pdyvy_stencil

                    ! MPML
                    if (k == 1) then

                        if (i <= 0 .or. i >= nx + 1 .or. j <= 0 .or. j >= ny + 1) then

                            ratiox = 0.0
                            ratioy = 0.0
                            ratioz = 0.0
                            if (i <= 0) then
                                ratiox = dampratio_left(j, k)
                            else if (i >= nx + 1) then
                                ratiox = dampratio_right(j, k)
                            end if
                            if (j <= 0) then
                                ratioy = dampratio_front(i, k)
                            else if (j >= ny + 1) then
                                ratioy = dampratio_back(i, k)
                            end if
                            ratioz = dampratio_bottom(i, j)

                            dax = idax(i, j, k) + ratioy*iday(i, j, k) + ratioz*idaz(i, j, k)
                            a = dax/kappaxi(i, j, k)
                            b = dax/kappaxi(i, j, k) + alphaxi(i, j, k)
                            damp_coef(dt, a, b, ax, bx)
                            kx = kappaxi(i, j, k)
                            damp_wavefield(memory_pdxvx, pdxvx, ax, bx, kx, i, j, k)

                            day = iday(i, j, k) + ratiox*idax(i, j, k) + ratioz*idaz(i, j, k)
                            a = day/kappayi(i, j, k)
                            b = day/kappayi(i, j, k) + alphayi(i, j, k)
                            damp_coef(dt, a, b, ay, by)
                            ky = kappayi(i, j, k)
                            damp_wavefield(memory_pdyvy, pdyvy, ay, by, ky, i, j, k)

                        end if

                        pdzvz = -c13(i, j, k)/c33(i, j, k)*pdxvx - c23(i, j, k)/c33(i, j, k)*pdyvy

                    else

                        pdzvz = idz*pdzvz_stencil*dz_scaling_i(k)

                        if (i <= 0 .or. i >= nx + 1 .or. j <= 0 .or. j >= ny + 1 .or. k >= nz + 1) then

                            ratiox = 0.0
                            ratioy = 0.0
                            ratioz = 0.0
                            if (i <= 0) then
                                ratiox = dampratio_left(j, k)
                            else if (i >= nx + 1) then
                                ratiox = dampratio_right(j, k)
                            end if
                            if (j <= 0) then
                                ratioy = dampratio_front(i, k)
                            else if (j >= ny + 1) then
                                ratioy = dampratio_back(i, k)
                            end if
                            ratioz = dampratio_bottom(i, j)

                            dax = idax(i, j, k) + ratioy*iday(i, j, k) + ratioz*idaz(i, j, k)
                            a = dax/kappaxi(i, j, k)
                            b = dax/kappaxi(i, j, k) + alphaxi(i, j, k)
                            damp_coef(dt, a, b, ax, bx)
                            kx = kappaxi(i, j, k)
                            damp_wavefield(memory_pdxvx, pdxvx, ax, bx, kx, i, j, k)

                            day = iday(i, j, k) + ratiox*idax(i, j, k) + ratioz*idaz(i, j, k)
                            a = day/kappayi(i, j, k)
                            b = day/kappayi(i, j, k) + alphayi(i, j, k)
                            damp_coef(dt, a, b, ay, by)
                            ky = kappayi(i, j, k)
                            damp_wavefield(memory_pdyvy, pdyvy, ay, by, ky, i, j, k)

                            daz = idaz(i, j, k) + ratiox*idax(i, j, k) + ratioy*iday(i, j, k)
                            a = daz/kappazi(i, j, k)
                            b = daz/kappazi(i, j, k) + alphazi(i, j, k)
                            damp_coef(dt, a, b, az, bz)
                            kz = kappazi(i, j, k)
                            damp_wavefield(memory_pdzvz, pdzvz, az, bz, kz, i, j, k)

                        end if

                    end if

                    stressxx(i, j, k) = stressxx(i, j, k) + dt*( &
                        c11(i, j, k)*pdxvx &
                        + c12(i, j, k)*pdyvy &
                        + c13(i, j, k)*pdzvz)

                    stressyy(i, j, k) = stressyy(i, j, k) + dt*( &
                        c12(i, j, k)*pdxvx &
                        + c22(i, j, k)*pdyvy &
                        + c23(i, j, k)*pdzvz)

                    stresszz(i, j, k) = stresszz(i, j, k) + dt*( &
                        c13(i, j, k)*pdxvx &
                        + c23(i, j, k)*pdyvy &
                        + c33(i, j, k)*pdzvz)

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1 - 1), nz2 - 1
            do j = ny1 - 1, ny2 - 1
                do i = nx1, nx2

                    pdzvy = idz*pdzvy_stencil*dz_scaling_h(k + 1)
                    pdyvz = idy*pdyvz_stencil

                    ! MPML
                    if (i <= 0 .or. i >= nx + 1 &
                            .or. j + 1 <= 1 .or. j + 1 >= ny + 1 .or. k + 1 >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i <= 0) then
                            ratiox = dampratio_left(j + 1, k + 1)
                        else if (i >= nx + 1) then
                            ratiox = dampratio_right(j + 1, k + 1)
                        end if
                        if (j + 1 <= 1) then
                            ratioy = dampratio_front(i, k + 1)
                        else if (j + 1 >= ny + 1) then
                            ratioy = dampratio_back(i, k + 1)
                        end if
                        ratioz = dampratio_bottom(i, j + 1)

                        day = hday(i, j + 1, k + 1) + ratiox*idax(i, j + 1, k + 1) + ratioz*hdaz(i, j + 1, k + 1)
                        a = day/kappayh(i, j + 1, k + 1)
                        b = day/kappayh(i, j + 1, k + 1) + alphayh(i, j + 1, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i, j + 1, k + 1)
                        damp_wavefield(memory_pdyvz, pdyvz, ay, by, ky, i, j + 1, k + 1)

                        daz = hdaz(i, j + 1, k + 1) + ratiox*idax(i, j + 1, k + 1) + ratioy*hday(i, j + 1, k + 1)
                        a = daz/kappazh(i, j + 1, k + 1)
                        b = daz/kappazh(i, j + 1, k + 1) + alphazh(i, j + 1, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i, j + 1, k + 1)
                        damp_wavefield(memory_pdzvy, pdzvy, az, bz, kz, i, j + 1, k + 1)

                    end if

                    stressyz(i, j + 1, k + 1) = stressyz(i, j + 1, k + 1) + dt*c44_eff_yz*(pdzvy + pdyvz)

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1 - 1), nz2 - 1
            do j = ny1, ny2
                do i = nx1 - 1, nx2 - 1

                    pdzvx = idz*pdzvx_stencil*dz_scaling_h(k + 1)
                    pdxvz = idx*pdxvz_stencil

                    ! MPML
                    if (i + 1 <= 1 .or. i + 1 >= nx + 1 &
                            .or. j <= 0 .or. j >= ny + 1 .or. k + 1 >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i + 1 <= 1) then
                            ratiox = dampratio_left(j, k + 1)
                        else if (i + 1 >= nx + 1) then
                            ratiox = dampratio_right(j, k + 1)
                        end if
                        if (j <= 0) then
                            ratioy = dampratio_front(i + 1, k + 1)
                        else if (j >= ny + 1) then
                            ratioy = dampratio_back(i + 1, k + 1)
                        end if
                        ratioz = dampratio_bottom(i + 1, j)

                        dax = hdax(i + 1, j, k + 1) + ratioy*iday(i + 1, j, k + 1) + ratioz*hdaz(i + 1, j, k + 1)
                        a = dax/kappaxh(i + 1, j, k + 1)
                        b = dax/kappaxh(i + 1, j, k + 1) + alphaxh(i + 1, j, k + 1)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxh(i + 1, j, k + 1)
                        damp_wavefield(memory_pdxvz, pdxvz, ax, bx, kx, i + 1, j, k + 1)

                        daz = hdaz(i + 1, j, k + 1) + ratiox*hdax(i + 1, j, k + 1) + ratioy*iday(i + 1, j, k + 1)
                        a = daz/kappazh(i + 1, j, k + 1)
                        b = daz/kappazh(i + 1, j, k + 1) + alphazh(i + 1, j, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i + 1, j, k + 1)
                        damp_wavefield(memory_pdzvx, pdzvx, az, bz, kz, i + 1, j, k + 1)

                    end if

                    stressxz(i + 1, j, k + 1) = stressxz(i + 1, j, k + 1) + dt*c55_eff_xz*(pdzvx + pdxvz)

                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1), nz2
            do j = ny1 - 1, ny2 - 1
                do i = nx1 - 1, nx2 - 1

                    pdyvx = idy*pdyvx_stencil
                    pdxvy = idx*pdxvy_stencil

                    ! MPML
                    if (i + 1 <= 1 .or. i + 1 >= nx + 1 &
                            .or. j + 1 <= 1 .or. j + 1 >= ny + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i + 1 <= 1) then
                            ratiox = dampratio_left(j + 1, k)
                        else if (i + 1 >= nx + 1) then
                            ratiox = dampratio_right(j + 1, k)
                        end if
                        if (j + 1 <= 1) then
                            ratioy = dampratio_front(i + 1, k)
                        else if (j + 1 >= ny + 1) then
                            ratioy = dampratio_back(i + 1, k)
                        end if
                        ratioz = dampratio_bottom(i + 1, j + 1)

                        dax = hdax(i + 1, j + 1, k) + ratioy*hday(i + 1, j + 1, k) + ratioz*idaz(i + 1, j + 1, k)
                        a = dax/kappaxh(i + 1, j + 1, k)
                        b = dax/kappaxh(i + 1, j + 1, k) + alphaxh(i + 1, j + 1, k)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxh(i + 1, j + 1, k)
                        damp_wavefield(memory_pdxvy, pdxvy, ax, bx, kx, i + 1, j + 1, k)

                        day = hday(i + 1, j + 1, k) + ratiox*hdax(i + 1, j + 1, k) + ratioz*idaz(i + 1, j + 1, k)
                        a = day/kappayh(i + 1, j + 1, k)
                        b = day/kappayh(i + 1, j + 1, k) + alphayh(i + 1, j + 1, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i + 1, j + 1, k)
                        damp_wavefield(memory_pdyvx, pdyvx, ay, by, ky, i + 1, j + 1, k)

                    end if

                    stressxy(i + 1, j + 1, k) = stressxy(i + 1, j + 1, k) + dt*c66_eff_xy*(pdyvx + pdxvy)

                end do
            end do
        end do
        !$omp end parallel do

        ! Free surface conditions
        !$omp parallel do private(k) schedule(auto)
        do k = 1, fdhalf

            if (1 - k >= nz1 .and. 1 - k <= nz2 .and. 1 + k >= nz1 - fdhalf + 1 .and. 1 + k <= nz2 + fdhalf) then

                ! The stress are mirror-symmetric w.r.t. the free surface, i.e., k = 1
                stresszz(:, :, 1 - k) = -stresszz(:, :, 1 + k)

                ! At the free surface, normal stress is strictly zero
                stresszz(:, :, 1) = 0.0

                ! The following are redundant but necessary for implementing near-surface explosion/moment tensor source
                ! In the case of flat free surface, there are no z derivatives of these field variables
                stressxx(:, :, 1 - k) = -stressxx(:, :, 1 + k)
                stressyy(:, :, 1 - k) = -stressyy(:, :, 1 + k)
                stressxy(:, :, 1 - k) = -stressxy(:, :, 1 + k)

            end if

            if (2 - k >= nz1 .and. 2 - k <= nz2 .and. 1 + k >= nz1 - fdhalf + 1 .and. 1 + k <= nz2 + fdhalf) then

                ! The shear stresses sigma_?z are mirror-symmetric
                stressyz(:, :, 2 - k) = -stressyz(:, :, 1 + k)
                stressxz(:, :, 2 - k) = -stressxz(:, :, 1 + k)

            end if

        end do
        !$omp end parallel do

        call commute_array_group(stressxx, fdhalf)
        call commute_array_group(stressyy, fdhalf)
        call commute_array_group(stresszz, fdhalf)
        call commute_array_group(stressyz, fdhalf)
        call commute_array_group(stressxz, fdhalf)
        call commute_array_group(stressxy, fdhalf)

        ! vx component
        !$omp parallel do private(i, j, k, pdxxx, pdyxy, pdzxz, pdxxy, pdyyy, pdzyz, pdxxz, pdyyz, pdzzz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1), nz2
            do j = ny1, ny2
                do i = nx1 - 1, nx2 - 1

                    pdxxx = idx*pdxxx_stencil
                    pdyxy = idy*pdyxy_stencil
                    pdzxz = idz*pdzxz_stencil*dz_scaling_i(k)

                    ! MPML
                    if (i + 1 <= 1 .or. i + 1 >= nx + 1 &
                            .or. j <= 0 .or. j >= ny + 1 .or. k >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i + 1 <= 1) then
                            ratiox = dampratio_left(j, k)
                        else if (i + 1 >= nx + 1) then
                            ratiox = dampratio_right(j, k)
                        end if
                        if (j <= 0) then
                            ratioy = dampratio_front(i + 1, k)
                        else if (j >= ny + 1) then
                            ratioy = dampratio_back(i + 1, k)
                        end if
                        ratioz = dampratio_bottom(i + 1, j)

                        dax = hdax(i + 1, j, k) + ratioy*iday(i + 1, j, k) + ratioz*idaz(i + 1, j, k)
                        a = dax/kappaxh(i + 1, j, k)
                        b = dax/kappaxh(i + 1, j, k) + alphaxh(i + 1, j, k)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxh(i + 1, j, k)
                        damp_wavefield(memory_pdxxx, pdxxx, ax, bx, kx, i + 1, j, k)

                        day = iday(i + 1, j, k) + ratiox*hdax(i + 1, j, k) + ratioz*idaz(i + 1, j, k)
                        a = day/kappayi(i + 1, j, k)
                        b = day/kappayi(i + 1, j, k) + alphayi(i + 1, j, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i + 1, j, k)
                        damp_wavefield(memory_pdyxy, pdyxy, ay, by, ky, i + 1, j, k)

                        daz = idaz(i + 1, j, k) + ratiox*hdax(i + 1, j, k) + ratioy*iday(i + 1, j, k)
                        a = daz/kappazi(i + 1, j, k)
                        b = daz/kappazi(i + 1, j, k) + alphazi(i + 1, j, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i + 1, j, k)
                        damp_wavefield(memory_pdzxz, pdzxz, az, bz, kz, i + 1, j, k)

                    end if

                    vx(i + 1, j, k) = vx(i + 1, j, k) + dt/rho_eff_x*(pdxxx + pdyxy + pdzxz)

                end do
            end do
        end do
        !$omp end parallel do

        ! vy component
        !$omp parallel do private(i, j, k, pdxxx, pdyxy, pdzxz, pdxxy, pdyyy, pdzyz, pdxxz, pdyyz, pdzzz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1), nz2
            do j = ny1 - 1, ny2 - 1
                do i = nx1, nx2

                    pdxxy = idx*pdxxy_stencil
                    pdyyy = idy*pdyyy_stencil
                    pdzyz = idz*pdzyz_stencil*dz_scaling_i(k)

                    ! MPML
                    if (i <= 0 .or. i >= nx + 1 &
                            .or. j + 1 <= 1 .or. j + 1 >= ny + 1 .or. k >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i <= 0) then
                            ratiox = dampratio_left(j + 1, k)
                        else if (i >= nx + 1) then
                            ratiox = dampratio_right(j + 1, k)
                        end if
                        if (j + 1 <= 1) then
                            ratioy = dampratio_front(i, k)
                        else if (j + 1 >= ny + 1) then
                            ratioy = dampratio_back(i, k)
                        end if
                        ratioz = dampratio_bottom(i, j + 1)

                        dax = idax(i, j + 1, k) + ratioy*hday(i, j + 1, k) + ratioz*idaz(i, j + 1, k)
                        a = dax/kappaxi(i, j + 1, k)
                        b = dax/kappaxi(i, j + 1, k) + alphaxi(i, j + 1, k)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxi(i, j + 1, k)
                        damp_wavefield(memory_pdxxy, pdxxy, ax, bx, kx, i, j + 1, k)

                        day = hday(i, j + 1, k) + ratiox*idax(i, j + 1, k) + ratioz*idaz(i, j + 1, k)
                        a = day/kappayh(i, j + 1, k)
                        b = day/kappayh(i, j + 1, k) + alphayh(i, j + 1, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i, j + 1, k)
                        damp_wavefield(memory_pdyyy, pdyyy, ay, by, ky, i, j + 1, k)

                        daz = idaz(i, j + 1, k) + ratiox*idax(i, j + 1, k) + ratioy*hday(i, j + 1, k)
                        a = daz/kappazi(i, j + 1, k)
                        b = daz/kappazi(i, j + 1, k) + alphazi(i, j + 1, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i, j + 1, k)
                        damp_wavefield(memory_pdzyz, pdzyz, az, bz, kz, i, j + 1, k)

                    end if

                    vy(i, j + 1, k) = vy(i, j + 1, k) + dt/rho_eff_y*(pdxxy + pdyyy + pdzyz)

                end do
            end do
        end do
        !$omp end parallel do

        ! vz component
        !$omp parallel do private(i, j, k, pdxxx, pdyxy, pdzxz, pdxxy, pdyyy, pdzyz, pdxxz, pdyyz, pdzzz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1 - 1), nz2 - 1
            do j = ny1, ny2
                do i = nx1, nx2

                    pdxxz = idx*pdxxz_stencil
                    pdyyz = idy*pdyyz_stencil
                    pdzzz = idz*pdzzz_stencil*dz_scaling_h(k + 1)

                    ! MPML
                    if (i <= 0 .or. i >= nx + 1 &
                            .or. j <= 0 .or. j >= ny + 1 .or. k + 1 >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0
                        if (i <= 0) then
                            ratiox = dampratio_left(j, k + 1)
                        else if (i >= nx + 1) then
                            ratiox = dampratio_right(j, k + 1)
                        end if
                        if (j <= 0) then
                            ratioy = dampratio_front(i, k + 1)
                        else if (j >= ny + 1) then
                            ratioy = dampratio_back(i, k + 1)
                        end if
                        ratioz = dampratio_bottom(i, j)

                        dax = idax(i, j, k + 1) + ratioy*iday(i, j, k + 1) + ratioz*hdaz(i, j, k + 1)
                        a = dax/kappaxi(i, j, k + 1)
                        b = dax/kappaxi(i, j, k + 1) + alphaxi(i, j, k + 1)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxi(i, j, k + 1)
                        damp_wavefield(memory_pdxxz, pdxxz, ax, bx, kx, i, j, k + 1)

                        day = iday(i, j, k + 1) + ratiox*idax(i, j, k + 1) + ratioz*hdaz(i, j, k + 1)
                        a = day/kappayi(i, j, k + 1)
                        b = day/kappayi(i, j, k + 1) + alphayi(i, j, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i, j, k + 1)
                        damp_wavefield(memory_pdyyz, pdyyz, ay, by, ky, i, j, k + 1)

                        daz = hdaz(i, j, k + 1) + ratiox*idax(i, j, k + 1) + ratioy*iday(i, j, k + 1)
                        a = daz/kappazh(i, j, k + 1)
                        b = daz/kappazh(i, j, k + 1) + alphazh(i, j, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i, j, k + 1)
                        damp_wavefield(memory_pdzzz, pdzzz, az, bz, kz, i, j, k + 1)

                    end if

                    vz(i, j, k + 1) = vz(i, j, k + 1) + dt/rho_eff_z*(pdxxz + pdyyz + pdzzz)

                end do
            end do
        end do
        !$omp end parallel do

    end subroutine

    !
    !> Source term
    !
    subroutine add_source(t)

        integer, intent(in) :: t

        integer :: sgx, sgy, sgz
        real :: polar, azimuth, amp
        integer :: k, nbeg, nend
        real :: m11, m12, m13, m22, m23, m33
        integer :: irx, iry, irz
        real :: dz_s
        real :: rho_s(1:1)

        do k = 1, sgmtr%ns

            if (yn_free_surface) then
                dz_s = dz_i(sgmtr%srcr(k)%gz)
            else
                dz_s = dz
            end if

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
                        amp = amp/(dx*dy*dz_s)

                        sgx = sgmtr%srcr(k)%gx
                        sgy = sgmtr%srcr(k)%gy
                        sgz = sgmtr%srcr(k)%gz
                        !$omp parallel do private(irx, iry, irz) collapse(3) schedule(auto)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxx(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxx(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressyy(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyy(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stresszz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stresszz(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                    case ('mt')
                        ! Moment tensor implemented as stress drop on stress components
                        m11 = sgmtr%srcr(k)%moment_tensor(1, 1)
                        m12 = sgmtr%srcr(k)%moment_tensor(1, 2)
                        m13 = sgmtr%srcr(k)%moment_tensor(1, 3)
                        m22 = sgmtr%srcr(k)%moment_tensor(2, 2)
                        m23 = sgmtr%srcr(k)%moment_tensor(2, 3)
                        m33 = sgmtr%srcr(k)%moment_tensor(3, 3)

                        amp = amp/(dx*dy*dz_s)

                        sgx = sgmtr%srcr(k)%gx
                        sgy = sgmtr%srcr(k)%gy
                        sgz = sgmtr%srcr(k)%gz
                        !$omp parallel do private(irx, iry, irz) collapse(3) schedule(auto)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxx(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxx(sgx + irx, sgy + iry, sgz + irz) - amp*m11 &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressyy(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyy(sgx + irx, sgy + iry, sgz + irz) - amp*m22 &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stresszz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stresszz(sgx + irx, sgy + iry, sgz + irz) - amp*m33 &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        sgx = sgmtr%srcr(k)%hx
                        sgy = sgmtr%srcr(k)%hy
                        sgz = sgmtr%srcr(k)%gz
                        !$omp parallel do private(irx, iry, irz) collapse(3) schedule(auto)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxy(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxy(sgx + irx, sgy + iry, sgz + irz) - amp*m12 &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        sgx = sgmtr%srcr(k)%hx
                        sgy = sgmtr%srcr(k)%gy
                        sgz = sgmtr%srcr(k)%hz
                        !$omp parallel do private(irx, iry, irz) collapse(3) schedule(auto)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxz(sgx + irx, sgy + iry, sgz + irz) - amp*m13 &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        sgx = sgmtr%srcr(k)%gx
                        sgy = sgmtr%srcr(k)%hy
                        sgz = sgmtr%srcr(k)%hz
                        !$omp parallel do private(irx, iry, irz) collapse(3) schedule(auto)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressyz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyz(sgx + irx, sgy + iry, sgz + irz) - amp*m23 &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                end select
            end if
        end do

    end subroutine

end module
