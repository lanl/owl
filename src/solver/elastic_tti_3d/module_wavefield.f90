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


module elastic_tti_3d_wavefield

    use libflit
    use elastic_tti_3d_vars
    use elastic_tti_3d_cfspml

    implicit none

    ! FD stencil
#include 'macro_fd_stencil.f90'

    ! Macros for MPML wavefield damping
#define damp_coef(dt, a_, b_, a, b) \
    a = -2.0*abs(dt)*a_/(2.0 + abs(dt)*b_); \
    b = (2.0 - abs(dt)*b_)/(2.0 + abs(dt)*b_);

#define damp_wavefield(w, v, da, db, dk, i, j, k) \
    w(i, j, k) = da*v + db*w(i, j, k); \
    v = (v + w(i, j, k))/dk;

contains

    !
    !> Update wavefields in 3D general anisotropic linear elastic media
    !
    subroutine update_wavefield(dt, &
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

        ! Arguments
        real, intent(in) :: dt
        real, allocatable, dimension(:, :, :), intent(inout) :: &
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
            memory_pdzzz_hxhyhz

        integer :: i, j, k
        real :: pdxvx, pdyvx, pdzvx
        real :: pdxvy, pdyvy, pdzvy
        real :: pdxvz, pdyvz, pdzvz
        real :: pdxxx, pdyxy, pdzxz
        real :: pdxxy, pdyyy, pdzyz
        real :: pdxxz, pdyyz, pdzzz
        real :: pd1, pd2, pd3, pd4, pd5, pd6
        real :: c11_ixiyiz, c12_ixiyiz, c13_ixiyiz, c14_ixiyiz, c15_ixiyiz, c16_ixiyiz
        real :: c22_ixiyiz, c23_ixiyiz, c24_ixiyiz, c25_ixiyiz, c26_ixiyiz
        real :: c33_ixiyiz, c34_ixiyiz, c35_ixiyiz, c36_ixiyiz
        real :: c44_ixiyiz, c45_ixiyiz, c46_ixiyiz
        real :: c55_ixiyiz, c56_ixiyiz
        real :: c66_ixiyiz
        real :: c11_hxhyiz, c12_hxhyiz, c13_hxhyiz, c14_hxhyiz, c15_hxhyiz, c16_hxhyiz
        real :: c22_hxhyiz, c23_hxhyiz, c24_hxhyiz, c25_hxhyiz, c26_hxhyiz
        real :: c33_hxhyiz, c34_hxhyiz, c35_hxhyiz, c36_hxhyiz
        real :: c44_hxhyiz, c45_hxhyiz, c46_hxhyiz
        real :: c55_hxhyiz, c56_hxhyiz
        real :: c66_hxhyiz
        real :: c11_hxiyhz, c12_hxiyhz, c13_hxiyhz, c14_hxiyhz, c15_hxiyhz, c16_hxiyhz
        real :: c22_hxiyhz, c23_hxiyhz, c24_hxiyhz, c25_hxiyhz, c26_hxiyhz
        real :: c33_hxiyhz, c34_hxiyhz, c35_hxiyhz, c36_hxiyhz
        real :: c44_hxiyhz, c45_hxiyhz, c46_hxiyhz
        real :: c55_hxiyhz, c56_hxiyhz
        real :: c66_hxiyhz
        real :: c11_ixhyhz, c12_ixhyhz, c13_ixhyhz, c14_ixhyhz, c15_ixhyhz, c16_ixhyhz
        real :: c22_ixhyhz, c23_ixhyhz, c24_ixhyhz, c25_ixhyhz, c26_ixhyhz
        real :: c33_ixhyhz, c34_ixhyhz, c35_ixhyhz, c36_ixhyhz
        real :: c44_ixhyhz, c45_ixhyhz, c46_ixhyhz
        real :: c55_ixhyhz, c56_ixhyhz
        real :: c66_ixhyhz
        real :: rho_hxiyiz, rho_ixhyiz, rho_ixiyhz, rho_hxhyhz
        real :: ratiox, ratioy, ratioz
        real :: dax, day, daz
        real :: a, b, ax, ay, az, bx, by, bz, kx, ky, kz

        ! Exchange boundary wavefields
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

        ! Stress set a (i, j, k)
        !$omp parallel do private(i, j, k, &
            !$omp pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, &
            !$omp c11_ixiyiz, c12_ixiyiz, c13_ixiyiz, c14_ixiyiz, c15_ixiyiz, c16_ixiyiz, &
            !$omp c22_ixiyiz, c23_ixiyiz, c24_ixiyiz, c25_ixiyiz, c26_ixiyiz, &
            !$omp c33_ixiyiz, c34_ixiyiz, c35_ixiyiz, c36_ixiyiz, &
            !$omp c44_ixiyiz, c45_ixiyiz, c46_ixiyiz, &
            !$omp c55_ixiyiz, c56_ixiyiz, &
            !$omp c66_ixiyiz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1, nx2

                    c11_ixiyiz = c11(i, j, k)
                    c12_ixiyiz = c12(i, j, k)
                    c13_ixiyiz = c13(i, j, k)
                    c14_ixiyiz = c14(i, j, k)
                    c15_ixiyiz = c15(i, j, k)
                    c16_ixiyiz = c16(i, j, k)
                    c22_ixiyiz = c22(i, j, k)
                    c23_ixiyiz = c23(i, j, k)
                    c24_ixiyiz = c24(i, j, k)
                    c25_ixiyiz = c25(i, j, k)
                    c26_ixiyiz = c26(i, j, k)
                    c33_ixiyiz = c33(i, j, k)
                    c34_ixiyiz = c34(i, j, k)
                    c35_ixiyiz = c35(i, j, k)
                    c36_ixiyiz = c36(i, j, k)
                    c44_ixiyiz = c44(i, j, k)
                    c45_ixiyiz = c45(i, j, k)
                    c46_ixiyiz = c46(i, j, k)
                    c55_ixiyiz = c55(i, j, k)
                    c56_ixiyiz = c56(i, j, k)
                    c66_ixiyiz = c66(i, j, k)

                    pdxvx = pdxw_stencil(vx_hxiyiz, i, j, k)/dx
                    pdyvx = pdyw_stencil(vx_ixhyiz, i, j, k)/dy
                    pdzvx = pdzw_stencil(vx_ixiyhz, i, j, k)/dz
                    pdxvy = pdxw_stencil(vy_hxiyiz, i, j, k)/dx
                    pdyvy = pdyw_stencil(vy_ixhyiz, i, j, k)/dy
                    pdzvy = pdzw_stencil(vy_ixiyhz, i, j, k)/dz
                    pdxvz = pdxw_stencil(vz_hxiyiz, i, j, k)/dx
                    pdyvz = pdyw_stencil(vz_ixhyiz, i, j, k)/dy
                    pdzvz = pdzw_stencil(vz_ixiyhz, i, j, k)/dz

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
                        damp_wavefield(memory_pdxvx_ixiyiz, pdxvx, ax, bx, kx, i, j, k)
                        damp_wavefield(memory_pdxvy_ixiyiz, pdxvy, ax, bx, kx, i, j, k)
                        damp_wavefield(memory_pdxvz_ixiyiz, pdxvz, ax, bx, kx, i, j, k)

                        day = iday(i, j, k) + ratiox*idax(i, j, k) + ratioz*idaz(i, j, k)
                        a = day/kappayi(i, j, k)
                        b = day/kappayi(i, j, k) + alphayi(i, j, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i, j, k)
                        damp_wavefield(memory_pdyvx_ixiyiz, pdyvx, ay, by, ky, i, j, k)
                        damp_wavefield(memory_pdyvy_ixiyiz, pdyvy, ay, by, ky, i, j, k)
                        damp_wavefield(memory_pdyvz_ixiyiz, pdyvz, ay, by, ky, i, j, k)

                        daz = idaz(i, j, k) + ratiox*idax(i, j, k) + ratioy*iday(i, j, k)
                        a = daz/kappazi(i, j, k)
                        b = daz/kappazi(i, j, k) + alphazi(i, j, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i, j, k)
                        damp_wavefield(memory_pdzvx_ixiyiz, pdzvx, az, bz, kz, i, j, k)
                        damp_wavefield(memory_pdzvy_ixiyiz, pdzvy, az, bz, kz, i, j, k)
                        damp_wavefield(memory_pdzvz_ixiyiz, pdzvz, az, bz, kz, i, j, k)

                    end if

                    pd1 = pdxvx
                    pd2 = pdyvy
                    pd3 = pdzvz
                    pd4 = pdzvy + pdyvz
                    pd5 = pdzvx + pdxvz
                    pd6 = pdyvx + pdxvy

                    stressxx_ixiyiz(i, j, k) = stressxx_ixiyiz(i, j, k) + dt*( &
                        c11_ixiyiz*pd1 &
                        + c12_ixiyiz*pd2 &
                        + c13_ixiyiz*pd3 &
                        + c14_ixiyiz*pd4 &
                        + c15_ixiyiz*pd5 &
                        + c16_ixiyiz*pd6)

                    stressyy_ixiyiz(i, j, k) = stressyy_ixiyiz(i, j, k) + dt*( &
                        c12_ixiyiz*pd1 &
                        + c22_ixiyiz*pd2 &
                        + c23_ixiyiz*pd3 &
                        + c24_ixiyiz*pd4 &
                        + c25_ixiyiz*pd5 &
                        + c26_ixiyiz*pd6)

                    stresszz_ixiyiz(i, j, k) = stresszz_ixiyiz(i, j, k) + dt*( &
                        c13_ixiyiz*pd1 &
                        + c23_ixiyiz*pd2 &
                        + c33_ixiyiz*pd3 &
                        + c34_ixiyiz*pd4 &
                        + c35_ixiyiz*pd5 &
                        + c36_ixiyiz*pd6)

                    stressyz_ixiyiz(i, j, k) = stressyz_ixiyiz(i, j, k) + dt*( &
                        c14_ixiyiz*pd1 &
                        + c24_ixiyiz*pd2 &
                        + c34_ixiyiz*pd3 &
                        + c44_ixiyiz*pd4 &
                        + c45_ixiyiz*pd5 &
                        + c46_ixiyiz*pd6)

                    stressxz_ixiyiz(i, j, k) = stressxz_ixiyiz(i, j, k) + dt*( &
                        c15_ixiyiz*pd1 &
                        + c25_ixiyiz*pd2 &
                        + c35_ixiyiz*pd3 &
                        + c45_ixiyiz*pd4 &
                        + c55_ixiyiz*pd5 &
                        + c56_ixiyiz*pd6)

                    stressxy_ixiyiz(i, j, k) = stressxy_ixiyiz(i, j, k) + dt*( &
                        c16_ixiyiz*pd1 &
                        + c26_ixiyiz*pd2 &
                        + c36_ixiyiz*pd3 &
                        + c46_ixiyiz*pd4 &
                        + c56_ixiyiz*pd5 &
                        + c66_ixiyiz*pd6)

                end do
            end do
        end do
        !$omp end parallel do

        ! Stress set b (i+-1/2, j+-1/2, k)
        !$omp parallel do private(i, j, k, &
            !$omp pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, &
            !$omp c11_hxhyiz, c12_hxhyiz, c13_hxhyiz, c14_hxhyiz, c15_hxhyiz, c16_hxhyiz, &
            !$omp c22_hxhyiz, c23_hxhyiz, c24_hxhyiz, c25_hxhyiz, c26_hxhyiz, &
            !$omp c33_hxhyiz, c34_hxhyiz, c35_hxhyiz, c36_hxhyiz, &
            !$omp c44_hxhyiz, c45_hxhyiz, c46_hxhyiz, &
            !$omp c55_hxhyiz, c56_hxhyiz, &
            !$omp c66_hxhyiz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1 - 1, ny2 - 1
                do i = nx1 - 1, nx2 - 1

                    c11_hxhyiz = 0.25*sum(c11(i:i + 1, j:j + 1, k))
                    c12_hxhyiz = 0.25*sum(c12(i:i + 1, j:j + 1, k))
                    c13_hxhyiz = 0.25*sum(c13(i:i + 1, j:j + 1, k))
                    c14_hxhyiz = 0.25*sum(c14(i:i + 1, j:j + 1, k))
                    c15_hxhyiz = 0.25*sum(c15(i:i + 1, j:j + 1, k))
                    c16_hxhyiz = 0.25*sum(c16(i:i + 1, j:j + 1, k))
                    c22_hxhyiz = 0.25*sum(c22(i:i + 1, j:j + 1, k))
                    c23_hxhyiz = 0.25*sum(c23(i:i + 1, j:j + 1, k))
                    c24_hxhyiz = 0.25*sum(c24(i:i + 1, j:j + 1, k))
                    c25_hxhyiz = 0.25*sum(c25(i:i + 1, j:j + 1, k))
                    c26_hxhyiz = 0.25*sum(c26(i:i + 1, j:j + 1, k))
                    c33_hxhyiz = 0.25*sum(c33(i:i + 1, j:j + 1, k))
                    c34_hxhyiz = 0.25*sum(c34(i:i + 1, j:j + 1, k))
                    c35_hxhyiz = 0.25*sum(c35(i:i + 1, j:j + 1, k))
                    c36_hxhyiz = 0.25*sum(c36(i:i + 1, j:j + 1, k))
                    c44_hxhyiz = 0.25*sum(c44(i:i + 1, j:j + 1, k))
                    c45_hxhyiz = 0.25*sum(c45(i:i + 1, j:j + 1, k))
                    c46_hxhyiz = 0.25*sum(c46(i:i + 1, j:j + 1, k))
                    c55_hxhyiz = 0.25*sum(c55(i:i + 1, j:j + 1, k))
                    c56_hxhyiz = 0.25*sum(c56(i:i + 1, j:j + 1, k))
                    c66_hxhyiz = 0.25*sum(c66(i:i + 1, j:j + 1, k))

                    pdxvx = pdxw_stencil(vx_ixhyiz, i, j + 1, k)/dx
                    pdyvx = pdyw_stencil(vx_hxiyiz, i + 1, j, k)/dy
                    pdzvx = pdzw_stencil(vx_hxhyhz, i + 1, j + 1, k)/dz
                    pdxvy = pdxw_stencil(vy_ixhyiz, i, j + 1, k)/dx
                    pdyvy = pdyw_stencil(vy_hxiyiz, i + 1, j, k)/dy
                    pdzvy = pdzw_stencil(vy_hxhyhz, i + 1, j + 1, k)/dz
                    pdxvz = pdxw_stencil(vz_ixhyiz, i, j + 1, k)/dx
                    pdyvz = pdyw_stencil(vz_hxiyiz, i + 1, j, k)/dy
                    pdzvz = pdzw_stencil(vz_hxhyhz, i + 1, j + 1, k)/dz

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
                        damp_wavefield(memory_pdxvx_hxhyiz, pdxvx, ax, bx, kx, i + 1, j + 1, k)
                        damp_wavefield(memory_pdxvy_hxhyiz, pdxvy, ax, bx, kx, i + 1, j + 1, k)
                        damp_wavefield(memory_pdxvz_hxhyiz, pdxvz, ax, bx, kx, i + 1, j + 1, k)

                        day = hday(i + 1, j + 1, k) + ratiox*hdax(i + 1, j + 1, k) + ratioz*idaz(i + 1, j + 1, k)
                        a = day/kappayh(i + 1, j + 1, k)
                        b = day/kappayh(i + 1, j + 1, k) + alphayh(i + 1, j + 1, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i + 1, j + 1, k)
                        damp_wavefield(memory_pdyvx_hxhyiz, pdyvx, ay, by, ky, i + 1, j + 1, k)
                        damp_wavefield(memory_pdyvy_hxhyiz, pdyvy, ay, by, ky, i + 1, j + 1, k)
                        damp_wavefield(memory_pdyvz_hxhyiz, pdyvz, ay, by, ky, i + 1, j + 1, k)

                        daz = idaz(i + 1, j + 1, k) + ratiox*hdax(i + 1, j + 1, k) + ratioy*hday(i + 1, j + 1, k)
                        a = daz/kappazi(i + 1, j + 1, k)
                        b = daz/kappazi(i + 1, j + 1, k) + alphazi(i + 1, j + 1, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i + 1, j + 1, k)
                        damp_wavefield(memory_pdzvx_hxhyiz, pdzvx, az, bz, kz, i + 1, j + 1, k)
                        damp_wavefield(memory_pdzvy_hxhyiz, pdzvy, az, bz, kz, i + 1, j + 1, k)
                        damp_wavefield(memory_pdzvz_hxhyiz, pdzvz, az, bz, kz, i + 1, j + 1, k)

                    end if

                    pd1 = pdxvx
                    pd2 = pdyvy
                    pd3 = pdzvz
                    pd4 = pdzvy + pdyvz
                    pd5 = pdzvx + pdxvz
                    pd6 = pdyvx + pdxvy

                    stressxx_hxhyiz(i + 1, j + 1, k) = stressxx_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c11_hxhyiz*pd1 &
                        + c12_hxhyiz*pd2 &
                        + c13_hxhyiz*pd3 &
                        + c14_hxhyiz*pd4 &
                        + c15_hxhyiz*pd5 &
                        + c16_hxhyiz*pd6)

                    stressyy_hxhyiz(i + 1, j + 1, k) = stressyy_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c12_hxhyiz*pd1 &
                        + c22_hxhyiz*pd2 &
                        + c23_hxhyiz*pd3 &
                        + c24_hxhyiz*pd4 &
                        + c25_hxhyiz*pd5 &
                        + c26_hxhyiz*pd6)

                    stresszz_hxhyiz(i + 1, j + 1, k) = stresszz_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c13_hxhyiz*pd1 &
                        + c23_hxhyiz*pd2 &
                        + c33_hxhyiz*pd3 &
                        + c34_hxhyiz*pd4 &
                        + c35_hxhyiz*pd5 &
                        + c36_hxhyiz*pd6)

                    stressyz_hxhyiz(i + 1, j + 1, k) = stressyz_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c14_hxhyiz*pd1 &
                        + c24_hxhyiz*pd2 &
                        + c34_hxhyiz*pd3 &
                        + c44_hxhyiz*pd4 &
                        + c45_hxhyiz*pd5 &
                        + c46_hxhyiz*pd6)

                    stressxz_hxhyiz(i + 1, j + 1, k) = stressxz_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c15_hxhyiz*pd1 &
                        + c25_hxhyiz*pd2 &
                        + c35_hxhyiz*pd3 &
                        + c45_hxhyiz*pd4 &
                        + c55_hxhyiz*pd5 &
                        + c56_hxhyiz*pd6)

                    stressxy_hxhyiz(i + 1, j + 1, k) = stressxy_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c16_hxhyiz*pd1 &
                        + c26_hxhyiz*pd2 &
                        + c36_hxhyiz*pd3 &
                        + c46_hxhyiz*pd4 &
                        + c56_hxhyiz*pd5 &
                        + c66_hxhyiz*pd6)

                end do
            end do
        end do
        !$omp end parallel do

        ! Stress set c (i+-1/2, j, k+-1/2)
        !$omp parallel do private(i, j, k, &
            !$omp pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, &
            !$omp c11_hxiyhz, c12_hxiyhz, c13_hxiyhz, c14_hxiyhz, c15_hxiyhz, c16_hxiyhz, &
            !$omp c22_hxiyhz, c23_hxiyhz, c24_hxiyhz, c25_hxiyhz, c26_hxiyhz, &
            !$omp c33_hxiyhz, c34_hxiyhz, c35_hxiyhz, c36_hxiyhz, &
            !$omp c44_hxiyhz, c45_hxiyhz, c46_hxiyhz, &
            !$omp c55_hxiyhz, c56_hxiyhz, &
            !$omp c66_hxiyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1 - 1, nz2 - 1
            do j = ny1, ny2
                do i = nx1 - 1, nx2 - 1

                    c11_hxiyhz = 0.25*sum(c11(i:i + 1, j, k:k + 1))
                    c12_hxiyhz = 0.25*sum(c12(i:i + 1, j, k:k + 1))
                    c13_hxiyhz = 0.25*sum(c13(i:i + 1, j, k:k + 1))
                    c14_hxiyhz = 0.25*sum(c14(i:i + 1, j, k:k + 1))
                    c15_hxiyhz = 0.25*sum(c15(i:i + 1, j, k:k + 1))
                    c16_hxiyhz = 0.25*sum(c16(i:i + 1, j, k:k + 1))
                    c22_hxiyhz = 0.25*sum(c22(i:i + 1, j, k:k + 1))
                    c23_hxiyhz = 0.25*sum(c23(i:i + 1, j, k:k + 1))
                    c24_hxiyhz = 0.25*sum(c24(i:i + 1, j, k:k + 1))
                    c25_hxiyhz = 0.25*sum(c25(i:i + 1, j, k:k + 1))
                    c26_hxiyhz = 0.25*sum(c26(i:i + 1, j, k:k + 1))
                    c33_hxiyhz = 0.25*sum(c33(i:i + 1, j, k:k + 1))
                    c34_hxiyhz = 0.25*sum(c34(i:i + 1, j, k:k + 1))
                    c35_hxiyhz = 0.25*sum(c35(i:i + 1, j, k:k + 1))
                    c36_hxiyhz = 0.25*sum(c36(i:i + 1, j, k:k + 1))
                    c44_hxiyhz = 0.25*sum(c44(i:i + 1, j, k:k + 1))
                    c45_hxiyhz = 0.25*sum(c45(i:i + 1, j, k:k + 1))
                    c46_hxiyhz = 0.25*sum(c46(i:i + 1, j, k:k + 1))
                    c55_hxiyhz = 0.25*sum(c55(i:i + 1, j, k:k + 1))
                    c56_hxiyhz = 0.25*sum(c56(i:i + 1, j, k:k + 1))
                    c66_hxiyhz = 0.25*sum(c66(i:i + 1, j, k:k + 1))

                    ! Spatial derivatives of particle velocity components
                    pdxvx = pdxw_stencil(vx_ixiyhz, i, j, k + 1)/dx
                    pdyvx = pdyw_stencil(vx_hxhyhz, i + 1, j, k + 1)/dy
                    pdzvx = pdzw_stencil(vx_hxiyiz, i + 1, j, k)/dz
                    pdxvy = pdxw_stencil(vy_ixiyhz, i, j, k + 1)/dx
                    pdyvy = pdyw_stencil(vy_hxhyhz, i + 1, j, k + 1)/dy
                    pdzvy = pdzw_stencil(vy_hxiyiz, i + 1, j, k)/dz
                    pdxvz = pdxw_stencil(vz_ixiyhz, i, j, k + 1)/dx
                    pdyvz = pdyw_stencil(vz_hxhyhz, i + 1, j, k + 1)/dy
                    pdzvz = pdzw_stencil(vz_hxiyiz, i + 1, j, k)/dz

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
                        damp_wavefield(memory_pdxvx_hxiyhz, pdxvx, ax, bx, kx, i + 1, j, k + 1)
                        damp_wavefield(memory_pdxvy_hxiyhz, pdxvy, ax, bx, kx, i + 1, j, k + 1)
                        damp_wavefield(memory_pdxvz_hxiyhz, pdxvz, ax, bx, kx, i + 1, j, k + 1)

                        day = iday(i + 1, j, k + 1) + ratiox*hdax(i + 1, j, k + 1) + ratioz*hdaz(i + 1, j, k + 1)
                        a = day/kappayi(i + 1, j, k + 1)
                        b = day/kappayi(i + 1, j, k + 1) + alphayi(i + 1, j, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i + 1, j, k + 1)
                        damp_wavefield(memory_pdyvx_hxiyhz, pdyvx, ay, by, ky, i + 1, j, k + 1)
                        damp_wavefield(memory_pdyvy_hxiyhz, pdyvy, ay, by, ky, i + 1, j, k + 1)
                        damp_wavefield(memory_pdyvz_hxiyhz, pdyvz, ay, by, ky, i + 1, j, k + 1)

                        daz = hdaz(i + 1, j, k + 1) + ratiox*hdax(i + 1, j, k + 1) + ratioy*iday(i + 1, j, k + 1)
                        a = daz/kappazh(i + 1, j, k + 1)
                        b = daz/kappazh(i + 1, j, k + 1) + alphazh(i + 1, j, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i + 1, j, k + 1)
                        damp_wavefield(memory_pdzvx_hxiyhz, pdzvx, az, bz, kz, i + 1, j, k + 1)
                        damp_wavefield(memory_pdzvy_hxiyhz, pdzvy, az, bz, kz, i + 1, j, k + 1)
                        damp_wavefield(memory_pdzvz_hxiyhz, pdzvz, az, bz, kz, i + 1, j, k + 1)

                    end if

                    pd1 = pdxvx
                    pd2 = pdyvy
                    pd3 = pdzvz
                    pd4 = pdzvy + pdyvz
                    pd5 = pdzvx + pdxvz
                    pd6 = pdyvx + pdxvy

                    stressxx_hxiyhz(i + 1, j, k + 1) = stressxx_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c11_hxiyhz*pd1 &
                        + c12_hxiyhz*pd2 &
                        + c13_hxiyhz*pd3 &
                        + c14_hxiyhz*pd4 &
                        + c15_hxiyhz*pd5 &
                        + c16_hxiyhz*pd6)

                    stressyy_hxiyhz(i + 1, j, k + 1) = stressyy_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c12_hxiyhz*pd1 &
                        + c22_hxiyhz*pd2 &
                        + c23_hxiyhz*pd3 &
                        + c24_hxiyhz*pd4 &
                        + c25_hxiyhz*pd5 &
                        + c26_hxiyhz*pd6)

                    stresszz_hxiyhz(i + 1, j, k + 1) = stresszz_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c13_hxiyhz*pd1 &
                        + c23_hxiyhz*pd2 &
                        + c33_hxiyhz*pd3 &
                        + c34_hxiyhz*pd4 &
                        + c35_hxiyhz*pd5 &
                        + c36_hxiyhz*pd6)

                    stressyz_hxiyhz(i + 1, j, k + 1) = stressyz_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c14_hxiyhz*pd1 &
                        + c24_hxiyhz*pd2 &
                        + c34_hxiyhz*pd3 &
                        + c44_hxiyhz*pd4 &
                        + c45_hxiyhz*pd5 &
                        + c46_hxiyhz*pd6)

                    stressxz_hxiyhz(i + 1, j, k + 1) = stressxz_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c15_hxiyhz*pd1 &
                        + c25_hxiyhz*pd2 &
                        + c35_hxiyhz*pd3 &
                        + c45_hxiyhz*pd4 &
                        + c55_hxiyhz*pd5 &
                        + c56_hxiyhz*pd6)

                    stressxy_hxiyhz(i + 1, j, k + 1) = stressxy_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c16_hxiyhz*pd1 &
                        + c26_hxiyhz*pd2 &
                        + c36_hxiyhz*pd3 &
                        + c46_hxiyhz*pd4 &
                        + c56_hxiyhz*pd5 &
                        + c66_hxiyhz*pd6)

                end do
            end do
        end do
        !$omp end parallel do

        ! Stress set d (i, j+-1/2, k+-1/2)
        !$omp parallel do private(i, j, k, &
            !$omp pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, &
            !$omp c11_ixhyhz, c12_ixhyhz, c13_ixhyhz, c14_ixhyhz, c15_ixhyhz, c16_ixhyhz, &
            !$omp c22_ixhyhz, c23_ixhyhz, c24_ixhyhz, c25_ixhyhz, c26_ixhyhz, &
            !$omp c33_ixhyhz, c34_ixhyhz, c35_ixhyhz, c36_ixhyhz, &
            !$omp c44_ixhyhz, c45_ixhyhz, c46_ixhyhz, &
            !$omp c55_ixhyhz, c56_ixhyhz, &
            !$omp c66_ixhyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1 - 1, nz2 - 1
            do j = ny1 - 1, ny2 - 1
                do i = nx1, nx2

                    c11_ixhyhz = 0.25*sum(c11(i, j:j + 1, k:k + 1))
                    c12_ixhyhz = 0.25*sum(c12(i, j:j + 1, k:k + 1))
                    c13_ixhyhz = 0.25*sum(c13(i, j:j + 1, k:k + 1))
                    c14_ixhyhz = 0.25*sum(c14(i, j:j + 1, k:k + 1))
                    c15_ixhyhz = 0.25*sum(c15(i, j:j + 1, k:k + 1))
                    c16_ixhyhz = 0.25*sum(c16(i, j:j + 1, k:k + 1))
                    c22_ixhyhz = 0.25*sum(c22(i, j:j + 1, k:k + 1))
                    c23_ixhyhz = 0.25*sum(c23(i, j:j + 1, k:k + 1))
                    c24_ixhyhz = 0.25*sum(c24(i, j:j + 1, k:k + 1))
                    c25_ixhyhz = 0.25*sum(c25(i, j:j + 1, k:k + 1))
                    c26_ixhyhz = 0.25*sum(c26(i, j:j + 1, k:k + 1))
                    c33_ixhyhz = 0.25*sum(c33(i, j:j + 1, k:k + 1))
                    c34_ixhyhz = 0.25*sum(c34(i, j:j + 1, k:k + 1))
                    c35_ixhyhz = 0.25*sum(c35(i, j:j + 1, k:k + 1))
                    c36_ixhyhz = 0.25*sum(c36(i, j:j + 1, k:k + 1))
                    c44_ixhyhz = 0.25*sum(c44(i, j:j + 1, k:k + 1))
                    c45_ixhyhz = 0.25*sum(c45(i, j:j + 1, k:k + 1))
                    c46_ixhyhz = 0.25*sum(c46(i, j:j + 1, k:k + 1))
                    c55_ixhyhz = 0.25*sum(c55(i, j:j + 1, k:k + 1))
                    c56_ixhyhz = 0.25*sum(c56(i, j:j + 1, k:k + 1))
                    c66_ixhyhz = 0.25*sum(c66(i, j:j + 1, k:k + 1))

                    pdxvx = pdxw_stencil(vx_hxhyhz, i, j + 1, k + 1)/dx
                    pdyvx = pdyw_stencil(vx_ixiyhz, i, j, k + 1)/dy
                    pdzvx = pdzw_stencil(vx_ixhyiz, i, j + 1, k)/dz
                    pdxvy = pdxw_stencil(vy_hxhyhz, i, j + 1, k + 1)/dx
                    pdyvy = pdyw_stencil(vy_ixiyhz, i, j, k + 1)/dy
                    pdzvy = pdzw_stencil(vy_ixhyiz, i, j + 1, k)/dz
                    pdxvz = pdxw_stencil(vz_hxhyhz, i, j + 1, k + 1)/dx
                    pdyvz = pdyw_stencil(vz_ixiyhz, i, j, k + 1)/dy
                    pdzvz = pdzw_stencil(vz_ixhyiz, i, j + 1, k)/dz

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

                        dax = idax(i, j + 1, k + 1) + ratioy*hday(i, j + 1, k + 1) + ratioz*hdaz(i, j + 1, k + 1)
                        a = dax/kappaxi(i, j + 1, k + 1)
                        b = dax/kappaxi(i, j + 1, k + 1) + alphaxi(i, j + 1, k + 1)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxi(i, j + 1, k + 1)
                        damp_wavefield(memory_pdxvx_ixhyhz, pdxvx, ax, bx, kx, i, j + 1, k + 1)
                        damp_wavefield(memory_pdxvy_ixhyhz, pdxvy, ax, bx, kx, i, j + 1, k + 1)
                        damp_wavefield(memory_pdxvz_ixhyhz, pdxvz, ax, bx, kx, i, j + 1, k + 1)

                        day = hday(i, j + 1, k + 1) + ratiox*idax(i, j + 1, k + 1) + ratioz*hdaz(i, j + 1, k + 1)
                        a = day/kappayh(i, j + 1, k + 1)
                        b = day/kappayh(i, j + 1, k + 1) + alphayh(i, j + 1, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i, j + 1, k + 1)
                        damp_wavefield(memory_pdyvx_ixhyhz, pdyvx, ay, by, ky, i, j + 1, k + 1)
                        damp_wavefield(memory_pdyvy_ixhyhz, pdyvy, ay, by, ky, i, j + 1, k + 1)
                        damp_wavefield(memory_pdyvz_ixhyhz, pdyvz, ay, by, ky, i, j + 1, k + 1)

                        daz = hdaz(i, j + 1, k + 1) + ratiox*idax(i, j + 1, k + 1) + ratioy*hday(i, j + 1, k + 1)
                        a = daz/kappazh(i, j + 1, k + 1)
                        b = daz/kappazh(i, j + 1, k + 1) + alphazh(i, j + 1, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i, j + 1, k + 1)
                        damp_wavefield(memory_pdzvx_ixhyhz, pdzvx, az, bz, kz, i, j + 1, k + 1)
                        damp_wavefield(memory_pdzvy_ixhyhz, pdzvy, az, bz, kz, i, j + 1, k + 1)
                        damp_wavefield(memory_pdzvz_ixhyhz, pdzvz, az, bz, kz, i, j + 1, k + 1)

                    end if

                    pd1 = pdxvx
                    pd2 = pdyvy
                    pd3 = pdzvz
                    pd4 = pdzvy + pdyvz
                    pd5 = pdzvx + pdxvz
                    pd6 = pdyvx + pdxvy

                    stressxx_ixhyhz(i, j + 1, k + 1) = stressxx_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c11_ixhyhz*pd1 &
                        + c12_ixhyhz*pd2 &
                        + c13_ixhyhz*pd3 &
                        + c14_ixhyhz*pd4 &
                        + c15_ixhyhz*pd5 &
                        + c16_ixhyhz*pd6)

                    stressyy_ixhyhz(i, j + 1, k + 1) = stressyy_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c12_ixhyhz*pd1 &
                        + c22_ixhyhz*pd2 &
                        + c23_ixhyhz*pd3 &
                        + c24_ixhyhz*pd4 &
                        + c25_ixhyhz*pd5 &
                        + c26_ixhyhz*pd6)

                    stresszz_ixhyhz(i, j + 1, k + 1) = stresszz_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c13_ixhyhz*pd1 &
                        + c23_ixhyhz*pd2 &
                        + c33_ixhyhz*pd3 &
                        + c34_ixhyhz*pd4 &
                        + c35_ixhyhz*pd5 &
                        + c36_ixhyhz*pd6)

                    stressyz_ixhyhz(i, j + 1, k + 1) = stressyz_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c14_ixhyhz*pd1 &
                        + c24_ixhyhz*pd2 &
                        + c34_ixhyhz*pd3 &
                        + c44_ixhyhz*pd4 &
                        + c45_ixhyhz*pd5 &
                        + c46_ixhyhz*pd6)

                    stressxz_ixhyhz(i, j + 1, k + 1) = stressxz_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c15_ixhyhz*pd1 &
                        + c25_ixhyhz*pd2 &
                        + c35_ixhyhz*pd3 &
                        + c45_ixhyhz*pd4 &
                        + c55_ixhyhz*pd5 &
                        + c56_ixhyhz*pd6)

                    stressxy_ixhyhz(i, j + 1, k + 1) = stressxy_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c16_ixhyhz*pd1 &
                        + c26_ixhyhz*pd2 &
                        + c36_ixhyhz*pd3 &
                        + c46_ixhyhz*pd4 &
                        + c56_ixhyhz*pd5 &
                        + c66_ixhyhz*pd6)

                end do
            end do
        end do
        !$omp end parallel do

        ! Exchange boundary wavefields
        call commute_array_group(stressxx_ixiyiz, fdhalf)
        call commute_array_group(stressyy_ixiyiz, fdhalf)
        call commute_array_group(stresszz_ixiyiz, fdhalf)
        call commute_array_group(stressxy_ixiyiz, fdhalf)
        call commute_array_group(stressxz_ixiyiz, fdhalf)
        call commute_array_group(stressyz_ixiyiz, fdhalf)
        call commute_array_group(stressxx_hxhyiz, fdhalf)
        call commute_array_group(stressyy_hxhyiz, fdhalf)
        call commute_array_group(stresszz_hxhyiz, fdhalf)
        call commute_array_group(stressxy_hxhyiz, fdhalf)
        call commute_array_group(stressxz_hxhyiz, fdhalf)
        call commute_array_group(stressyz_hxhyiz, fdhalf)
        call commute_array_group(stressxx_hxiyhz, fdhalf)
        call commute_array_group(stressyy_hxiyhz, fdhalf)
        call commute_array_group(stresszz_hxiyhz, fdhalf)
        call commute_array_group(stressxy_hxiyhz, fdhalf)
        call commute_array_group(stressxz_hxiyhz, fdhalf)
        call commute_array_group(stressyz_hxiyhz, fdhalf)
        call commute_array_group(stressxx_ixhyhz, fdhalf)
        call commute_array_group(stressyy_ixhyhz, fdhalf)
        call commute_array_group(stresszz_ixhyhz, fdhalf)
        call commute_array_group(stressxy_ixhyhz, fdhalf)
        call commute_array_group(stressxz_ixhyhz, fdhalf)
        call commute_array_group(stressyz_ixhyhz, fdhalf)

        ! Particle velocity set a (i+-1/2, j, k)
        !$omp parallel do private(i, j, k, &
            !$omp pdxxx, pdyxy, pdzxz, &
            !$omp pdxxy, pdyyy, pdzyz, &
            !$omp pdxxz, pdyyz, pdzzz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, &
            !$omp rho_hxiyiz, rho_ixhyiz, rho_ixiyhz, rho_hxhyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1 - 1, nx2 - 1

                    rho_hxiyiz = 0.5*sum(rho(i:i + 1, j, k))

                    pdxxx = pdxw_stencil(stressxx_ixiyiz, i, j, k)/dx
                    pdxxy = pdxw_stencil(stressxy_ixiyiz, i, j, k)/dx
                    pdxxz = pdxw_stencil(stressxz_ixiyiz, i, j, k)/dx
                    pdyxy = pdyw_stencil(stressxy_hxhyiz, i + 1, j, k)/dy
                    pdyyy = pdyw_stencil(stressyy_hxhyiz, i + 1, j, k)/dy
                    pdyyz = pdyw_stencil(stressyz_hxhyiz, i + 1, j, k)/dy
                    pdzxz = pdzw_stencil(stressxz_hxiyhz, i + 1, j, k)/dz
                    pdzyz = pdzw_stencil(stressyz_hxiyhz, i + 1, j, k)/dz
                    pdzzz = pdzw_stencil(stresszz_hxiyhz, i + 1, j, k)/dz

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
                        damp_wavefield(memory_pdxxx_hxiyiz, pdxxx, ax, bx, kx, i + 1, j, k)
                        damp_wavefield(memory_pdxxy_hxiyiz, pdxxy, ax, bx, kx, i + 1, j, k)
                        damp_wavefield(memory_pdxxz_hxiyiz, pdxxz, ax, bx, kx, i + 1, j, k)

                        day = iday(i + 1, j, k) + ratiox*hdax(i + 1, j, k) + ratioz*idaz(i + 1, j, k)
                        a = day/kappayi(i + 1, j, k)
                        b = day/kappayi(i + 1, j, k) + alphayi(i + 1, j, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i + 1, j, k)
                        damp_wavefield(memory_pdyxy_hxiyiz, pdyxy, ay, by, ky, i + 1, j, k)
                        damp_wavefield(memory_pdyyy_hxiyiz, pdyyy, ay, by, ky, i + 1, j, k)
                        damp_wavefield(memory_pdyyz_hxiyiz, pdyyz, ay, by, ky, i + 1, j, k)

                        daz = idaz(i + 1, j, k) + ratiox*hdax(i + 1, j, k) + ratioy*iday(i + 1, j, k)
                        a = daz/kappazi(i + 1, j, k)
                        b = daz/kappazi(i + 1, j, k) + alphazi(i + 1, j, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i + 1, j, k)
                        damp_wavefield(memory_pdzxz_hxiyiz, pdzxz, az, bz, kz, i + 1, j, k)
                        damp_wavefield(memory_pdzyz_hxiyiz, pdzyz, az, bz, kz, i + 1, j, k)
                        damp_wavefield(memory_pdzzz_hxiyiz, pdzzz, az, bz, kz, i + 1, j, k)

                    end if

                    pd1 = pdxxx &
                        + pdyxy &
                        + pdzxz
                    pd2 = pdxxy &
                        + pdyyy &
                        + pdzyz
                    pd3 = pdxxz &
                        + pdyyz &
                        + pdzzz

                    vx_hxiyiz(i + 1, j, k) = vx_hxiyiz(i + 1, j, k) + dt/rho_hxiyiz*pd1
                    vy_hxiyiz(i + 1, j, k) = vy_hxiyiz(i + 1, j, k) + dt/rho_hxiyiz*pd2
                    vz_hxiyiz(i + 1, j, k) = vz_hxiyiz(i + 1, j, k) + dt/rho_hxiyiz*pd3

                end do
            end do
        end do
        !$omp end parallel do

        ! Particle velocity set b (i, j+-1/2, k)
        !$omp parallel do private(i, j, k, &
            !$omp pdxxx, pdyxy, pdzxz, &
            !$omp pdxxy, pdyyy, pdzyz, &
            !$omp pdxxz, pdyyz, pdzzz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, &
            !$omp rho_hxiyiz, rho_ixhyiz, rho_ixiyhz, rho_hxhyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1 - 1, ny2 - 1
                do i = nx1, nx2

                    rho_ixhyiz = 0.5*sum(rho(i, j:j + 1, k))

                    pdxxx = pdxw_stencil(stressxx_hxhyiz, i, j + 1, k)/dx
                    pdxxy = pdxw_stencil(stressxy_hxhyiz, i, j + 1, k)/dx
                    pdxxz = pdxw_stencil(stressxz_hxhyiz, i, j + 1, k)/dx
                    pdyxy = pdyw_stencil(stressxy_ixiyiz, i, j, k)/dy
                    pdyyy = pdyw_stencil(stressyy_ixiyiz, i, j, k)/dy
                    pdyyz = pdyw_stencil(stressyz_ixiyiz, i, j, k)/dy
                    pdzxz = pdzw_stencil(stressxz_ixhyhz, i, j + 1, k)/dz
                    pdzyz = pdzw_stencil(stressyz_ixhyhz, i, j + 1, k)/dz
                    pdzzz = pdzw_stencil(stresszz_ixhyhz, i, j + 1, k)/dz

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
                        damp_wavefield(memory_pdxxx_ixhyiz, pdxxx, ax, bx, kx, i, j + 1, k)
                        damp_wavefield(memory_pdxxy_ixhyiz, pdxxy, ax, bx, kx, i, j + 1, k)
                        damp_wavefield(memory_pdxxz_ixhyiz, pdxxz, ax, bx, kx, i, j + 1, k)

                        day = hday(i, j + 1, k) + ratiox*idax(i, j + 1, k) + ratioz*idaz(i, j + 1, k)
                        a = day/kappayh(i, j + 1, k)
                        b = day/kappayh(i, j + 1, k) + alphayh(i, j + 1, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i, j + 1, k)
                        damp_wavefield(memory_pdyxy_ixhyiz, pdyxy, ay, by, ky, i, j + 1, k)
                        damp_wavefield(memory_pdyyy_ixhyiz, pdyyy, ay, by, ky, i, j + 1, k)
                        damp_wavefield(memory_pdyyz_ixhyiz, pdyyz, ay, by, ky, i, j + 1, k)

                        daz = idaz(i, j + 1, k) + ratiox*idax(i, j + 1, k) + ratioy*hday(i, j + 1, k)
                        a = daz/kappazi(i, j + 1, k)
                        b = daz/kappazi(i, j + 1, k) + alphazi(i, j + 1, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i, j + 1, k)
                        damp_wavefield(memory_pdzxz_ixhyiz, pdzxz, az, bz, kz, i, j + 1, k)
                        damp_wavefield(memory_pdzyz_ixhyiz, pdzyz, az, bz, kz, i, j + 1, k)
                        damp_wavefield(memory_pdzzz_ixhyiz, pdzzz, az, bz, kz, i, j + 1, k)

                    end if

                    pd1 = pdxxx &
                        + pdyxy &
                        + pdzxz
                    pd2 = pdxxy &
                        + pdyyy &
                        + pdzyz
                    pd3 = pdxxz &
                        + pdyyz &
                        + pdzzz

                    vx_ixhyiz(i, j + 1, k) = vx_ixhyiz(i, j + 1, k) + dt/rho_ixhyiz*pd1
                    vy_ixhyiz(i, j + 1, k) = vy_ixhyiz(i, j + 1, k) + dt/rho_ixhyiz*pd2
                    vz_ixhyiz(i, j + 1, k) = vz_ixhyiz(i, j + 1, k) + dt/rho_ixhyiz*pd3

                end do
            end do
        end do
        !$omp end parallel do

        ! Particle velocity set c (i, j, k+-1/2)
        !$omp parallel do private(i, j, k, &
            !$omp pdxxx, pdyxy, pdzxz, &
            !$omp pdxxy, pdyyy, pdzyz, &
            !$omp pdxxz, pdyyz, pdzzz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, &
            !$omp rho_hxiyiz, rho_ixhyiz, rho_ixiyhz, rho_hxhyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1 - 1, nz2 - 1
            do j = ny1, ny2
                do i = nx1, nx2

                    rho_ixiyhz = 0.5*sum(rho(i, j, k:k + 1))

                    pdxxx = pdxw_stencil(stressxx_hxiyhz, i, j, k + 1)/dx
                    pdxxy = pdxw_stencil(stressxy_hxiyhz, i, j, k + 1)/dx
                    pdxxz = pdxw_stencil(stressxz_hxiyhz, i, j, k + 1)/dx
                    pdyxy = pdyw_stencil(stressxy_ixhyhz, i, j, k + 1)/dy
                    pdyyy = pdyw_stencil(stressyy_ixhyhz, i, j, k + 1)/dy
                    pdyyz = pdyw_stencil(stressyz_ixhyhz, i, j, k + 1)/dy
                    pdzxz = pdzw_stencil(stressxz_ixiyiz, i, j, k)/dz
                    pdzyz = pdzw_stencil(stressyz_ixiyiz, i, j, k)/dz
                    pdzzz = pdzw_stencil(stresszz_ixiyiz, i, j, k)/dz

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
                        damp_wavefield(memory_pdxxx_ixiyhz, pdxxx, ax, bx, kx, i, j, k + 1)
                        damp_wavefield(memory_pdxxy_ixiyhz, pdxxy, ax, bx, kx, i, j, k + 1)
                        damp_wavefield(memory_pdxxz_ixiyhz, pdxxz, ax, bx, kx, i, j, k + 1)

                        day = iday(i, j, k + 1) + ratiox*idax(i, j, k + 1) + ratioz*hdaz(i, j, k + 1)
                        a = day/kappayi(i, j, k + 1)
                        b = day/kappayi(i, j, k + 1) + alphayi(i, j, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i, j, k + 1)
                        damp_wavefield(memory_pdyxy_ixiyhz, pdyxy, ay, by, ky, i, j, k + 1)
                        damp_wavefield(memory_pdyyy_ixiyhz, pdyyy, ay, by, ky, i, j, k + 1)
                        damp_wavefield(memory_pdyyz_ixiyhz, pdyyz, ay, by, ky, i, j, k + 1)

                        daz = hdaz(i, j, k + 1) + ratiox*idax(i, j, k + 1) + ratioy*iday(i, j, k + 1)
                        a = daz/kappazh(i, j, k + 1)
                        b = daz/kappazh(i, j, k + 1) + alphazh(i, j, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i, j, k + 1)
                        damp_wavefield(memory_pdzxz_ixiyhz, pdzxz, az, bz, kz, i, j, k + 1)
                        damp_wavefield(memory_pdzyz_ixiyhz, pdzyz, az, bz, kz, i, j, k + 1)
                        damp_wavefield(memory_pdzzz_ixiyhz, pdzzz, az, bz, kz, i, j, k + 1)

                    end if

                    pd1 = pdxxx &
                        + pdyxy &
                        + pdzxz
                    pd2 = pdxxy &
                        + pdyyy &
                        + pdzyz
                    pd3 = pdxxz &
                        + pdyyz &
                        + pdzzz

                    vx_ixiyhz(i, j, k + 1) = vx_ixiyhz(i, j, k + 1) + dt/rho_ixiyhz*pd1
                    vy_ixiyhz(i, j, k + 1) = vy_ixiyhz(i, j, k + 1) + dt/rho_ixiyhz*pd2
                    vz_ixiyhz(i, j, k + 1) = vz_ixiyhz(i, j, k + 1) + dt/rho_ixiyhz*pd3

                end do
            end do
        end do
        !$omp end parallel do

        ! Particle velocity set d (i+-1/2, j+-1/2, k+-1/2)
        !$omp parallel do private(i, j, k, &
            !$omp pdxxx, pdyxy, pdzxz, &
            !$omp pdxxy, pdyyy, pdzyz, &
            !$omp pdxxz, pdyyz, pdzzz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, &
            !$omp rho_hxiyiz, rho_ixhyiz, rho_ixiyhz, rho_hxhyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = nz1 - 1, nz2 - 1
            do j = ny1 - 1, ny2 - 1
                do i = nx1 - 1, nx2 - 1

                    rho_hxhyhz = 0.125*sum(rho(i:i + 1, j:j + 1, k:k + 1))

                    pdxxx = pdxw_stencil(stressxx_ixhyhz, i, j + 1, k + 1)/dx
                    pdxxy = pdxw_stencil(stressxy_ixhyhz, i, j + 1, k + 1)/dx
                    pdxxz = pdxw_stencil(stressxz_ixhyhz, i, j + 1, k + 1)/dx
                    pdyxy = pdyw_stencil(stressxy_hxiyhz, i + 1, j, k + 1)/dy
                    pdyyy = pdyw_stencil(stressyy_hxiyhz, i + 1, j, k + 1)/dy
                    pdyyz = pdyw_stencil(stressyz_hxiyhz, i + 1, j, k + 1)/dy
                    pdzxz = pdzw_stencil(stressxz_hxhyiz, i + 1, j + 1, k)/dz
                    pdzyz = pdzw_stencil(stressyz_hxhyiz, i + 1, j + 1, k)/dz
                    pdzzz = pdzw_stencil(stresszz_hxhyiz, i + 1, j + 1, k)/dz

                    ! MPML
                    if (i + 1 <= 1 .or. i + 1 >= nx + 1 &
                            .or. j + 1 <= 1 .or. j + 1 >= ny + 1 &
                            .or. k + 1 <= 1 .or. k + 1 >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0

                        if (i + 1 <= 1) then
                            ratiox = dampratio_left(j + 1, k + 1)
                        else if (i + 1 >= nx + 1) then
                            ratiox = dampratio_right(j + 1, k + 1)
                        end if
                        if (j + 1 <= 1) then
                            ratioy = dampratio_front(i + 1, k + 1)
                        else if (j + 1 >= ny + 1) then
                            ratioy = dampratio_back(i + 1, k + 1)
                        end if
                        if (k + 1 <= 1) then
                            ratioz = dampratio_top(i + 1, j + 1)
                        else if (k + 1 >= nz + 1) then
                            ratioz = dampratio_bottom(i + 1, j + 1)
                        end if

                        dax = hdax(i + 1, j + 1, k + 1) + ratioy*hday(i + 1, j + 1, k + 1) + ratioz*hdaz(i + 1, j + 1, k + 1)
                        a = dax/kappaxh(i + 1, j + 1, k + 1)
                        b = dax/kappaxh(i + 1, j + 1, k + 1) + alphaxh(i + 1, j + 1, k + 1)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxh(i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdxxx_hxhyhz, pdxxx, ax, bx, kx, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdxxy_hxhyhz, pdxxy, ax, bx, kx, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdxxz_hxhyhz, pdxxz, ax, bx, kx, i + 1, j + 1, k + 1)

                        day = hday(i + 1, j + 1, k + 1) + ratiox*hdax(i + 1, j + 1, k + 1) + ratioz*hdaz(i + 1, j + 1, k + 1)
                        a = day/kappayh(i + 1, j + 1, k + 1)
                        b = day/kappayh(i + 1, j + 1, k + 1) + alphayh(i + 1, j + 1, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdyxy_hxhyhz, pdyxy, ay, by, ky, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdyyy_hxhyhz, pdyyy, ay, by, ky, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdyyz_hxhyhz, pdyyz, ay, by, ky, i + 1, j + 1, k + 1)

                        daz = hdaz(i + 1, j + 1, k + 1) + ratiox*hdax(i + 1, j + 1, k + 1) + ratioy*hday(i + 1, j + 1, k + 1)
                        a = daz/kappazh(i + 1, j + 1, k + 1)
                        b = daz/kappazh(i + 1, j + 1, k + 1) + alphazh(i + 1, j + 1, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdzxz_hxhyhz, pdzxz, az, bz, kz, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdzyz_hxhyhz, pdzyz, az, bz, kz, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdzzz_hxhyhz, pdzzz, az, bz, kz, i + 1, j + 1, k + 1)

                    end if

                    pd1 = pdxxx &
                        + pdyxy &
                        + pdzxz
                    pd2 = pdxxy &
                        + pdyyy &
                        + pdzyz
                    pd3 = pdxxz &
                        + pdyyz &
                        + pdzzz

                    vx_hxhyhz(i + 1, j + 1, k + 1) = vx_hxhyhz(i + 1, j + 1, k + 1) + dt/rho_hxhyhz*pd1
                    vy_hxhyhz(i + 1, j + 1, k + 1) = vy_hxhyhz(i + 1, j + 1, k + 1) + dt/rho_hxhyhz*pd2
                    vz_hxhyhz(i + 1, j + 1, k + 1) = vz_hxhyhz(i + 1, j + 1, k + 1) + dt/rho_hxhyhz*pd3

                end do
            end do
        end do
        !$omp end parallel do

    end subroutine

    !
    !> Update wavefields in 3D general anisotropic linear elastic media
    !> with a topographic free surface
    !
    subroutine update_wavefield_free_surface(dt, &
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

        ! Arguments
        real, intent(in) :: dt
        real, allocatable, dimension(:, :, :), intent(inout) :: &
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
            memory_pdzzz_hxhyhz

        integer :: i, j, k
        real :: pdxvx, pdyvx, pdzvx
        real :: pdxvy, pdyvy, pdzvy
        real :: pdxvz, pdyvz, pdzvz
        real :: pdxxx, pdyxy, pdzxx, pdzxy, pdzxz
        real :: pdxxy, pdyyy, pdzyy, pdzyz
        real :: pdxxz, pdyyz, pdzzz
        real :: pd1, pd2, pd3, pd4, pd5, pd6
        real :: tpaz, tpbz, tpcz, eta
        real :: c11_ixiyiz, c12_ixiyiz, c13_ixiyiz, c14_ixiyiz, c15_ixiyiz, c16_ixiyiz
        real :: c22_ixiyiz, c23_ixiyiz, c24_ixiyiz, c25_ixiyiz, c26_ixiyiz
        real :: c33_ixiyiz, c34_ixiyiz, c35_ixiyiz, c36_ixiyiz
        real :: c44_ixiyiz, c45_ixiyiz, c46_ixiyiz
        real :: c55_ixiyiz, c56_ixiyiz
        real :: c66_ixiyiz
        real :: c11_hxhyiz, c12_hxhyiz, c13_hxhyiz, c14_hxhyiz, c15_hxhyiz, c16_hxhyiz
        real :: c22_hxhyiz, c23_hxhyiz, c24_hxhyiz, c25_hxhyiz, c26_hxhyiz
        real :: c33_hxhyiz, c34_hxhyiz, c35_hxhyiz, c36_hxhyiz
        real :: c44_hxhyiz, c45_hxhyiz, c46_hxhyiz
        real :: c55_hxhyiz, c56_hxhyiz
        real :: c66_hxhyiz
        real :: c11_hxiyhz, c12_hxiyhz, c13_hxiyhz, c14_hxiyhz, c15_hxiyhz, c16_hxiyhz
        real :: c22_hxiyhz, c23_hxiyhz, c24_hxiyhz, c25_hxiyhz, c26_hxiyhz
        real :: c33_hxiyhz, c34_hxiyhz, c35_hxiyhz, c36_hxiyhz
        real :: c44_hxiyhz, c45_hxiyhz, c46_hxiyhz
        real :: c55_hxiyhz, c56_hxiyhz
        real :: c66_hxiyhz
        real :: c11_ixhyhz, c12_ixhyhz, c13_ixhyhz, c14_ixhyhz, c15_ixhyhz, c16_ixhyhz
        real :: c22_ixhyhz, c23_ixhyhz, c24_ixhyhz, c25_ixhyhz, c26_ixhyhz
        real :: c33_ixhyhz, c34_ixhyhz, c35_ixhyhz, c36_ixhyhz
        real :: c44_ixhyhz, c45_ixhyhz, c46_ixhyhz
        real :: c55_ixhyhz, c56_ixhyhz
        real :: c66_ixhyhz
        real :: rho_hxiyiz, rho_ixhyiz, rho_ixiyhz, rho_hxhyhz
        real :: ratiox, ratioy, ratioz
        real :: dax, day, daz
        real :: a, b, ax, ay, az, bx, by, bz, kx, ky, kz

        ! Set particle velocities above the free surface to zero
        !$omp parallel do private(k) schedule(auto)
        do k = 1, fdhalf
            if (1 - k >= nz1 .and. 1 - k <= nz2) then
                vx_hxiyiz(:, :, 1 - k) = 0.0
                vy_hxiyiz(:, :, 1 - k) = 0.0
                vz_hxiyiz(:, :, 1 - k) = 0.0
                vx_ixhyiz(:, :, 1 - k) = 0.0
                vy_ixhyiz(:, :, 1 - k) = 0.0
                vz_ixhyiz(:, :, 1 - k) = 0.0
            end if
            if (2 - k >= nz1 .and. 2 - k <= nz2) then
                vx_ixiyhz(:, :, 2 - k) = 0.0
                vy_ixiyhz(:, :, 2 - k) = 0.0
                vz_ixiyhz(:, :, 2 - k) = 0.0
                vx_hxhyhz(:, :, 2 - k) = 0.0
                vy_hxhyhz(:, :, 2 - k) = 0.0
                vz_hxhyhz(:, :, 2 - k) = 0.0
            end if
        end do
        !$omp end parallel do

        ! Exchange boundary wavefields
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

        ! Stress set a (i, j, k)
        !$omp parallel do private(i, j, k, &
            !$omp pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, tpaz, tpbz, tpcz, eta, &
            !$omp c11_ixiyiz, c12_ixiyiz, c13_ixiyiz, c14_ixiyiz, c15_ixiyiz, c16_ixiyiz, &
            !$omp c22_ixiyiz, c23_ixiyiz, c24_ixiyiz, c25_ixiyiz, c26_ixiyiz, &
            !$omp c33_ixiyiz, c34_ixiyiz, c35_ixiyiz, c36_ixiyiz, &
            !$omp c44_ixiyiz, c45_ixiyiz, c46_ixiyiz, &
            !$omp c55_ixiyiz, c56_ixiyiz, &
            !$omp c66_ixiyiz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1), nz2
            do j = ny1, ny2
                do i = nx1, nx2

                    eta = eta_zz_i(k)
                    tpaz = max(0.0, eta_max - eta)/(topo_ixiy(i, j) + depth_max)*slopex_ixiy(i, j)
                    tpbz = max(0.0, eta_max - eta)/(topo_ixiy(i, j) + depth_max)*slopey_ixiy(i, j)
                    tpcz = eta_max/(topo_ixiy(i, j) + depth_max)

                    c11_ixiyiz = c11(i, j, k)
                    c12_ixiyiz = c12(i, j, k)
                    c13_ixiyiz = c13(i, j, k)
                    c14_ixiyiz = c14(i, j, k)
                    c15_ixiyiz = c15(i, j, k)
                    c16_ixiyiz = c16(i, j, k)
                    c22_ixiyiz = c22(i, j, k)
                    c23_ixiyiz = c23(i, j, k)
                    c24_ixiyiz = c24(i, j, k)
                    c25_ixiyiz = c25(i, j, k)
                    c26_ixiyiz = c26(i, j, k)
                    c33_ixiyiz = c33(i, j, k)
                    c34_ixiyiz = c34(i, j, k)
                    c35_ixiyiz = c35(i, j, k)
                    c36_ixiyiz = c36(i, j, k)
                    c44_ixiyiz = c44(i, j, k)
                    c45_ixiyiz = c45(i, j, k)
                    c46_ixiyiz = c46(i, j, k)
                    c55_ixiyiz = c55(i, j, k)
                    c56_ixiyiz = c56(i, j, k)
                    c66_ixiyiz = c66(i, j, k)

                    pdxvx = pdxw_stencil(vx_hxiyiz, i, j, k)/dx
                    pdyvx = pdyw_stencil(vx_ixhyiz, i, j, k)/dy
                    pdxvy = pdxw_stencil(vy_hxiyiz, i, j, k)/dx
                    pdyvy = pdyw_stencil(vy_ixhyiz, i, j, k)/dy
                    pdxvz = pdxw_stencil(vz_hxiyiz, i, j, k)/dx
                    pdyvz = pdyw_stencil(vz_ixhyiz, i, j, k)/dy

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

                            dax = idax(i, j, k) + ratioy*iday(i, j, k) + ratioz*idaz(i, j, k)
                            a = dax/kappaxi(i, j, k)
                            b = dax/kappaxi(i, j, k) + alphaxi(i, j, k)
                            damp_coef(dt, a, b, ax, bx)
                            kx = kappaxi(i, j, k)
                            damp_wavefield(memory_pdxvx_ixiyiz, pdxvx, ax, bx, kx, i, j, k)
                            damp_wavefield(memory_pdxvy_ixiyiz, pdxvy, ax, bx, kx, i, j, k)
                            damp_wavefield(memory_pdxvz_ixiyiz, pdxvz, ax, bx, kx, i, j, k)

                            day = iday(i, j, k) + ratiox*idax(i, j, k) + ratioz*idaz(i, j, k)
                            a = day/kappayi(i, j, k)
                            b = day/kappayi(i, j, k) + alphayi(i, j, k)
                            damp_coef(dt, a, b, ay, by)
                            ky = kappayi(i, j, k)
                            damp_wavefield(memory_pdyvx_ixiyiz, pdyvx, ay, by, ky, i, j, k)
                            damp_wavefield(memory_pdyvy_ixiyiz, pdyvy, ay, by, ky, i, j, k)
                            damp_wavefield(memory_pdyvz_ixiyiz, pdyvz, ay, by, ky, i, j, k)

                        end if

                        pdzvx = coefiix(i, j)%array(1, 1)*pdxvx &
                            + coefiix(i, j)%array(1, 2)*pdxvy &
                            + coefiix(i, j)%array(1, 3)*pdxvz &
                            + coefiiy(i, j)%array(1, 1)*pdyvx &
                            + coefiiy(i, j)%array(1, 2)*pdyvy &
                            + coefiiy(i, j)%array(1, 3)*pdyvz
                        pdzvy = coefiix(i, j)%array(2, 1)*pdxvx &
                            + coefiix(i, j)%array(2, 2)*pdxvy &
                            + coefiix(i, j)%array(2, 3)*pdxvz &
                            + coefiiy(i, j)%array(2, 1)*pdyvx &
                            + coefiiy(i, j)%array(2, 2)*pdyvy &
                            + coefiiy(i, j)%array(2, 3)*pdyvz
                        pdzvz = coefiix(i, j)%array(3, 1)*pdxvx &
                            + coefiix(i, j)%array(3, 2)*pdxvy &
                            + coefiix(i, j)%array(3, 3)*pdxvz &
                            + coefiiy(i, j)%array(3, 1)*pdyvx &
                            + coefiiy(i, j)%array(3, 2)*pdyvy &
                            + coefiiy(i, j)%array(3, 3)*pdyvz

                    else

                        pdzvx = pdzw_stencil(vx_ixiyhz, i, j, k)/dz*eta_dz_scaling_i(k)
                        pdzvy = pdzw_stencil(vy_ixiyhz, i, j, k)/dz*eta_dz_scaling_i(k)
                        pdzvz = pdzw_stencil(vz_ixiyhz, i, j, k)/dz*eta_dz_scaling_i(k)

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
                            if (k >= nz + 1) then
                                ratioz = dampratio_bottom(i, j)
                            end if

                            dax = idax(i, j, k) + ratioy*iday(i, j, k) + ratioz*idaz(i, j, k)
                            a = dax/kappaxi(i, j, k)
                            b = dax/kappaxi(i, j, k) + alphaxi(i, j, k)
                            damp_coef(dt, a, b, ax, bx)
                            kx = kappaxi(i, j, k)
                            damp_wavefield(memory_pdxvx_ixiyiz, pdxvx, ax, bx, kx, i, j, k)
                            damp_wavefield(memory_pdxvy_ixiyiz, pdxvy, ax, bx, kx, i, j, k)
                            damp_wavefield(memory_pdxvz_ixiyiz, pdxvz, ax, bx, kx, i, j, k)

                            day = iday(i, j, k) + ratiox*idax(i, j, k) + ratioz*idaz(i, j, k)
                            a = day/kappayi(i, j, k)
                            b = day/kappayi(i, j, k) + alphayi(i, j, k)
                            damp_coef(dt, a, b, ay, by)
                            ky = kappayi(i, j, k)
                            damp_wavefield(memory_pdyvx_ixiyiz, pdyvx, ay, by, ky, i, j, k)
                            damp_wavefield(memory_pdyvy_ixiyiz, pdyvy, ay, by, ky, i, j, k)
                            damp_wavefield(memory_pdyvz_ixiyiz, pdyvz, ay, by, ky, i, j, k)

                            daz = idaz(i, j, k) + ratiox*idax(i, j, k) + ratioy*iday(i, j, k)
                            a = daz/kappazi(i, j, k)
                            b = daz/kappazi(i, j, k) + alphazi(i, j, k)
                            damp_coef(dt, a, b, az, bz)
                            kz = kappazi(i, j, k)
                            damp_wavefield(memory_pdzvx_ixiyiz, pdzvx, az, bz, kz, i, j, k)
                            damp_wavefield(memory_pdzvy_ixiyiz, pdzvy, az, bz, kz, i, j, k)
                            damp_wavefield(memory_pdzvz_ixiyiz, pdzvz, az, bz, kz, i, j, k)

                        end if

                    end if

                    pd1 = pdxvx + tpaz*pdzvx
                    pd2 = pdyvy + tpbz*pdzvy
                    pd3 = tpcz*pdzvz
                    pd4 = tpcz*pdzvy + pdyvz + tpbz*pdzvz
                    pd5 = tpcz*pdzvx + pdxvz + tpaz*pdzvz
                    pd6 = pdyvx + tpbz*pdzvx + pdxvy + tpaz*pdzvy

                    stressxx_ixiyiz(i, j, k) = stressxx_ixiyiz(i, j, k) + dt*( &
                        c11_ixiyiz*pd1 &
                        + c12_ixiyiz*pd2 &
                        + c13_ixiyiz*pd3 &
                        + c14_ixiyiz*pd4 &
                        + c15_ixiyiz*pd5 &
                        + c16_ixiyiz*pd6)

                    stressyy_ixiyiz(i, j, k) = stressyy_ixiyiz(i, j, k) + dt*( &
                        c12_ixiyiz*pd1 &
                        + c22_ixiyiz*pd2 &
                        + c23_ixiyiz*pd3 &
                        + c24_ixiyiz*pd4 &
                        + c25_ixiyiz*pd5 &
                        + c26_ixiyiz*pd6)

                    stresszz_ixiyiz(i, j, k) = stresszz_ixiyiz(i, j, k) + dt*( &
                        c13_ixiyiz*pd1 &
                        + c23_ixiyiz*pd2 &
                        + c33_ixiyiz*pd3 &
                        + c34_ixiyiz*pd4 &
                        + c35_ixiyiz*pd5 &
                        + c36_ixiyiz*pd6)

                    stressyz_ixiyiz(i, j, k) = stressyz_ixiyiz(i, j, k) + dt*( &
                        c14_ixiyiz*pd1 &
                        + c24_ixiyiz*pd2 &
                        + c34_ixiyiz*pd3 &
                        + c44_ixiyiz*pd4 &
                        + c45_ixiyiz*pd5 &
                        + c46_ixiyiz*pd6)

                    stressxz_ixiyiz(i, j, k) = stressxz_ixiyiz(i, j, k) + dt*( &
                        c15_ixiyiz*pd1 &
                        + c25_ixiyiz*pd2 &
                        + c35_ixiyiz*pd3 &
                        + c45_ixiyiz*pd4 &
                        + c55_ixiyiz*pd5 &
                        + c56_ixiyiz*pd6)

                    stressxy_ixiyiz(i, j, k) = stressxy_ixiyiz(i, j, k) + dt*( &
                        c16_ixiyiz*pd1 &
                        + c26_ixiyiz*pd2 &
                        + c36_ixiyiz*pd3 &
                        + c46_ixiyiz*pd4 &
                        + c56_ixiyiz*pd5 &
                        + c66_ixiyiz*pd6)

                end do
            end do
        end do
        !$omp end parallel do

        ! Stress set b (i+-1/2, j+-1/2, k)
        !$omp parallel do private(i, j, k, &
            !$omp pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, tpaz, tpbz, tpcz, eta, &
            !$omp c11_hxhyiz, c12_hxhyiz, c13_hxhyiz, c14_hxhyiz, c15_hxhyiz, c16_hxhyiz, &
            !$omp c22_hxhyiz, c23_hxhyiz, c24_hxhyiz, c25_hxhyiz, c26_hxhyiz, &
            !$omp c33_hxhyiz, c34_hxhyiz, c35_hxhyiz, c36_hxhyiz, &
            !$omp c44_hxhyiz, c45_hxhyiz, c46_hxhyiz, &
            !$omp c55_hxhyiz, c56_hxhyiz, &
            !$omp c66_hxhyiz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1), nz2
            do j = ny1 - 1, ny2 - 1
                do i = nx1 - 1, nx2 - 1

                    eta = eta_zz_i(k)
                    tpaz = max(0.0, eta_max - eta)/(topo_hxhy(i + 1, j + 1) + depth_max)*slopex_hxhy(i + 1, j + 1)
                    tpbz = max(0.0, eta_max - eta)/(topo_hxhy(i + 1, j + 1) + depth_max)*slopey_hxhy(i + 1, j + 1)
                    tpcz = eta_max/(topo_hxhy(i + 1, j + 1) + depth_max)

                    c11_hxhyiz = 0.25*sum(c11(i:i + 1, j:j + 1, k))
                    c12_hxhyiz = 0.25*sum(c12(i:i + 1, j:j + 1, k))
                    c13_hxhyiz = 0.25*sum(c13(i:i + 1, j:j + 1, k))
                    c14_hxhyiz = 0.25*sum(c14(i:i + 1, j:j + 1, k))
                    c15_hxhyiz = 0.25*sum(c15(i:i + 1, j:j + 1, k))
                    c16_hxhyiz = 0.25*sum(c16(i:i + 1, j:j + 1, k))
                    c22_hxhyiz = 0.25*sum(c22(i:i + 1, j:j + 1, k))
                    c23_hxhyiz = 0.25*sum(c23(i:i + 1, j:j + 1, k))
                    c24_hxhyiz = 0.25*sum(c24(i:i + 1, j:j + 1, k))
                    c25_hxhyiz = 0.25*sum(c25(i:i + 1, j:j + 1, k))
                    c26_hxhyiz = 0.25*sum(c26(i:i + 1, j:j + 1, k))
                    c33_hxhyiz = 0.25*sum(c33(i:i + 1, j:j + 1, k))
                    c34_hxhyiz = 0.25*sum(c34(i:i + 1, j:j + 1, k))
                    c35_hxhyiz = 0.25*sum(c35(i:i + 1, j:j + 1, k))
                    c36_hxhyiz = 0.25*sum(c36(i:i + 1, j:j + 1, k))
                    c44_hxhyiz = 0.25*sum(c44(i:i + 1, j:j + 1, k))
                    c45_hxhyiz = 0.25*sum(c45(i:i + 1, j:j + 1, k))
                    c46_hxhyiz = 0.25*sum(c46(i:i + 1, j:j + 1, k))
                    c55_hxhyiz = 0.25*sum(c55(i:i + 1, j:j + 1, k))
                    c56_hxhyiz = 0.25*sum(c56(i:i + 1, j:j + 1, k))
                    c66_hxhyiz = 0.25*sum(c66(i:i + 1, j:j + 1, k))

                    pdxvx = pdxw_stencil(vx_ixhyiz, i, j + 1, k)/dx
                    pdyvx = pdyw_stencil(vx_hxiyiz, i + 1, j, k)/dy
                    pdxvy = pdxw_stencil(vy_ixhyiz, i, j + 1, k)/dx
                    pdyvy = pdyw_stencil(vy_hxiyiz, i + 1, j, k)/dy
                    pdxvz = pdxw_stencil(vz_ixhyiz, i, j + 1, k)/dx
                    pdyvz = pdyw_stencil(vz_hxiyiz, i + 1, j, k)/dy

                    if (k == 1) then

                        ! MPML
                        if (i + 1 <= 1 .or. i + 1 >= nx + 1 .or. j + 1 <= 1 .or. j + 1 >= ny + 1) then

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
                            if (k >= nz + 1) then
                                ratioz = dampratio_bottom(i + 1, j + 1)
                            end if

                            dax = hdax(i + 1, j + 1, k) + ratioy*hday(i + 1, j + 1, k) + ratioz*idaz(i + 1, j + 1, k)
                            a = dax/kappaxh(i + 1, j + 1, k)
                            b = dax/kappaxh(i + 1, j + 1, k) + alphaxh(i + 1, j + 1, k)
                            damp_coef(dt, a, b, ax, bx)
                            kx = kappaxh(i + 1, j + 1, k)
                            damp_wavefield(memory_pdxvx_hxhyiz, pdxvx, ax, bx, kx, i + 1, j + 1, k)
                            damp_wavefield(memory_pdxvy_hxhyiz, pdxvy, ax, bx, kx, i + 1, j + 1, k)
                            damp_wavefield(memory_pdxvz_hxhyiz, pdxvz, ax, bx, kx, i + 1, j + 1, k)

                            day = hday(i + 1, j + 1, k) + ratiox*hdax(i + 1, j + 1, k) + ratioz*idaz(i + 1, j + 1, k)
                            a = day/kappayh(i + 1, j + 1, k)
                            b = day/kappayh(i + 1, j + 1, k) + alphayh(i + 1, j + 1, k)
                            damp_coef(dt, a, b, ay, by)
                            ky = kappayh(i + 1, j + 1, k)
                            damp_wavefield(memory_pdyvx_hxhyiz, pdyvx, ay, by, ky, i + 1, j + 1, k)
                            damp_wavefield(memory_pdyvy_hxhyiz, pdyvy, ay, by, ky, i + 1, j + 1, k)
                            damp_wavefield(memory_pdyvz_hxhyiz, pdyvz, ay, by, ky, i + 1, j + 1, k)

                        end if

                        pdzvx = coefhhx(i + 1, j + 1)%array(1, 1)*pdxvx &
                            + coefhhx(i + 1, j + 1)%array(1, 2)*pdxvy &
                            + coefhhx(i + 1, j + 1)%array(1, 3)*pdxvz &
                            + coefhhy(i + 1, j + 1)%array(1, 1)*pdyvx &
                            + coefhhy(i + 1, j + 1)%array(1, 2)*pdyvy &
                            + coefhhy(i + 1, j + 1)%array(1, 3)*pdyvz
                        pdzvy = coefhhx(i + 1, j + 1)%array(2, 1)*pdxvx &
                            + coefhhx(i + 1, j + 1)%array(2, 2)*pdxvy &
                            + coefhhx(i + 1, j + 1)%array(2, 3)*pdxvz &
                            + coefhhy(i + 1, j + 1)%array(2, 1)*pdyvx &
                            + coefhhy(i + 1, j + 1)%array(2, 2)*pdyvy &
                            + coefhhy(i + 1, j + 1)%array(2, 3)*pdyvz
                        pdzvz = coefhhx(i + 1, j + 1)%array(3, 1)*pdxvx &
                            + coefhhx(i + 1, j + 1)%array(3, 2)*pdxvy &
                            + coefhhx(i + 1, j + 1)%array(3, 3)*pdxvz &
                            + coefhhy(i + 1, j + 1)%array(3, 1)*pdyvx &
                            + coefhhy(i + 1, j + 1)%array(3, 2)*pdyvy &
                            + coefhhy(i + 1, j + 1)%array(3, 3)*pdyvz

                    else

                        pdzvx = pdzw_stencil(vx_hxhyhz, i + 1, j + 1, k)/dz*eta_dz_scaling_i(k)
                        pdzvy = pdzw_stencil(vy_hxhyhz, i + 1, j + 1, k)/dz*eta_dz_scaling_i(k)
                        pdzvz = pdzw_stencil(vz_hxhyhz, i + 1, j + 1, k)/dz*eta_dz_scaling_i(k)

                        ! MPML
                        if (i + 1 <= 1 .or. i + 1 >= nx + 1 .or. j + 1 <= 1 .or. j + 1 >= ny + 1 .or. k >= nz + 1) then

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
                            if (k >= nz + 1) then
                                ratioz = dampratio_bottom(i + 1, j + 1)
                            end if

                            dax = hdax(i + 1, j + 1, k) + ratioy*hday(i + 1, j + 1, k) + ratioz*idaz(i + 1, j + 1, k)
                            a = dax/kappaxh(i + 1, j + 1, k)
                            b = dax/kappaxh(i + 1, j + 1, k) + alphaxh(i + 1, j + 1, k)
                            damp_coef(dt, a, b, ax, bx)
                            kx = kappaxh(i + 1, j + 1, k)
                            damp_wavefield(memory_pdxvx_hxhyiz, pdxvx, ax, bx, kx, i + 1, j + 1, k)
                            damp_wavefield(memory_pdxvy_hxhyiz, pdxvy, ax, bx, kx, i + 1, j + 1, k)
                            damp_wavefield(memory_pdxvz_hxhyiz, pdxvz, ax, bx, kx, i + 1, j + 1, k)

                            day = hday(i + 1, j + 1, k) + ratiox*hdax(i + 1, j + 1, k) + ratioz*idaz(i + 1, j + 1, k)
                            a = day/kappayh(i + 1, j + 1, k)
                            b = day/kappayh(i + 1, j + 1, k) + alphayh(i + 1, j + 1, k)
                            damp_coef(dt, a, b, ay, by)
                            ky = kappayh(i + 1, j + 1, k)
                            damp_wavefield(memory_pdyvx_hxhyiz, pdyvx, ay, by, ky, i + 1, j + 1, k)
                            damp_wavefield(memory_pdyvy_hxhyiz, pdyvy, ay, by, ky, i + 1, j + 1, k)
                            damp_wavefield(memory_pdyvz_hxhyiz, pdyvz, ay, by, ky, i + 1, j + 1, k)

                            daz = idaz(i + 1, j + 1, k) + ratiox*hdax(i + 1, j + 1, k) + ratioy*hday(i + 1, j + 1, k)
                            a = daz/kappazi(i + 1, j + 1, k)
                            b = daz/kappazi(i + 1, j + 1, k) + alphazi(i + 1, j + 1, k)
                            damp_coef(dt, a, b, az, bz)
                            kz = kappazi(i + 1, j + 1, k)
                            damp_wavefield(memory_pdzvx_hxhyiz, pdzvx, az, bz, kz, i + 1, j + 1, k)
                            damp_wavefield(memory_pdzvy_hxhyiz, pdzvy, az, bz, kz, i + 1, j + 1, k)
                            damp_wavefield(memory_pdzvz_hxhyiz, pdzvz, az, bz, kz, i + 1, j + 1, k)

                        end if

                    end if

                    pd1 = pdxvx + tpaz*pdzvx
                    pd2 = pdyvy + tpbz*pdzvy
                    pd3 = tpcz*pdzvz
                    pd4 = tpcz*pdzvy + pdyvz + tpbz*pdzvz
                    pd5 = tpcz*pdzvx + pdxvz + tpaz*pdzvz
                    pd6 = pdyvx + tpbz*pdzvx + pdxvy + tpaz*pdzvy

                    stressxx_hxhyiz(i + 1, j + 1, k) = stressxx_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c11_hxhyiz*pd1 &
                        + c12_hxhyiz*pd2 &
                        + c13_hxhyiz*pd3 &
                        + c14_hxhyiz*pd4 &
                        + c15_hxhyiz*pd5 &
                        + c16_hxhyiz*pd6)

                    stressyy_hxhyiz(i + 1, j + 1, k) = stressyy_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c12_hxhyiz*pd1 &
                        + c22_hxhyiz*pd2 &
                        + c23_hxhyiz*pd3 &
                        + c24_hxhyiz*pd4 &
                        + c25_hxhyiz*pd5 &
                        + c26_hxhyiz*pd6)

                    stresszz_hxhyiz(i + 1, j + 1, k) = stresszz_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c13_hxhyiz*pd1 &
                        + c23_hxhyiz*pd2 &
                        + c33_hxhyiz*pd3 &
                        + c34_hxhyiz*pd4 &
                        + c35_hxhyiz*pd5 &
                        + c36_hxhyiz*pd6)

                    stressyz_hxhyiz(i + 1, j + 1, k) = stressyz_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c14_hxhyiz*pd1 &
                        + c24_hxhyiz*pd2 &
                        + c34_hxhyiz*pd3 &
                        + c44_hxhyiz*pd4 &
                        + c45_hxhyiz*pd5 &
                        + c46_hxhyiz*pd6)

                    stressxz_hxhyiz(i + 1, j + 1, k) = stressxz_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c15_hxhyiz*pd1 &
                        + c25_hxhyiz*pd2 &
                        + c35_hxhyiz*pd3 &
                        + c45_hxhyiz*pd4 &
                        + c55_hxhyiz*pd5 &
                        + c56_hxhyiz*pd6)

                    stressxy_hxhyiz(i + 1, j + 1, k) = stressxy_hxhyiz(i + 1, j + 1, k) + dt*( &
                        c16_hxhyiz*pd1 &
                        + c26_hxhyiz*pd2 &
                        + c36_hxhyiz*pd3 &
                        + c46_hxhyiz*pd4 &
                        + c56_hxhyiz*pd5 &
                        + c66_hxhyiz*pd6)

                end do
            end do
        end do
        !$omp end parallel do

        ! Stress set c (i+-1/2, j, k+-1/2)
        !$omp parallel do private(i, j, k, &
            !$omp pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, tpaz, tpbz, tpcz, eta, &
            !$omp c11_hxiyhz, c12_hxiyhz, c13_hxiyhz, c14_hxiyhz, c15_hxiyhz, c16_hxiyhz, &
            !$omp c22_hxiyhz, c23_hxiyhz, c24_hxiyhz, c25_hxiyhz, c26_hxiyhz, &
            !$omp c33_hxiyhz, c34_hxiyhz, c35_hxiyhz, c36_hxiyhz, &
            !$omp c44_hxiyhz, c45_hxiyhz, c46_hxiyhz, &
            !$omp c55_hxiyhz, c56_hxiyhz, &
            !$omp c66_hxiyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1 - 1), nz2 - 1
            do j = ny1, ny2
                do i = nx1 - 1, nx2 - 1

                    eta = eta_zz_h(k + 1)
                    tpaz = max(0.0, eta_max - eta)/(topo_hxiy(i + 1, j) + depth_max)*slopex_hxiy(i + 1, j)
                    tpbz = max(0.0, eta_max - eta)/(topo_hxiy(i + 1, j) + depth_max)*slopey_hxiy(i + 1, j)
                    tpcz = eta_max/(topo_hxiy(i + 1, j) + depth_max)

                    c11_hxiyhz = 0.25*sum(c11(i:i + 1, j, k:k + 1))
                    c12_hxiyhz = 0.25*sum(c12(i:i + 1, j, k:k + 1))
                    c13_hxiyhz = 0.25*sum(c13(i:i + 1, j, k:k + 1))
                    c14_hxiyhz = 0.25*sum(c14(i:i + 1, j, k:k + 1))
                    c15_hxiyhz = 0.25*sum(c15(i:i + 1, j, k:k + 1))
                    c16_hxiyhz = 0.25*sum(c16(i:i + 1, j, k:k + 1))
                    c22_hxiyhz = 0.25*sum(c22(i:i + 1, j, k:k + 1))
                    c23_hxiyhz = 0.25*sum(c23(i:i + 1, j, k:k + 1))
                    c24_hxiyhz = 0.25*sum(c24(i:i + 1, j, k:k + 1))
                    c25_hxiyhz = 0.25*sum(c25(i:i + 1, j, k:k + 1))
                    c26_hxiyhz = 0.25*sum(c26(i:i + 1, j, k:k + 1))
                    c33_hxiyhz = 0.25*sum(c33(i:i + 1, j, k:k + 1))
                    c34_hxiyhz = 0.25*sum(c34(i:i + 1, j, k:k + 1))
                    c35_hxiyhz = 0.25*sum(c35(i:i + 1, j, k:k + 1))
                    c36_hxiyhz = 0.25*sum(c36(i:i + 1, j, k:k + 1))
                    c44_hxiyhz = 0.25*sum(c44(i:i + 1, j, k:k + 1))
                    c45_hxiyhz = 0.25*sum(c45(i:i + 1, j, k:k + 1))
                    c46_hxiyhz = 0.25*sum(c46(i:i + 1, j, k:k + 1))
                    c55_hxiyhz = 0.25*sum(c55(i:i + 1, j, k:k + 1))
                    c56_hxiyhz = 0.25*sum(c56(i:i + 1, j, k:k + 1))
                    c66_hxiyhz = 0.25*sum(c66(i:i + 1, j, k:k + 1))

                    ! Spatial derivatives of particle velocity components
                    pdxvx = pdxw_stencil(vx_ixiyhz, i, j, k + 1)/dx
                    pdyvx = pdyw_stencil(vx_hxhyhz, i + 1, j, k + 1)/dy
                    pdzvx = pdzw_stencil(vx_hxiyiz, i + 1, j, k)/dz*eta_dz_scaling_h(k + 1)
                    pdxvy = pdxw_stencil(vy_ixiyhz, i, j, k + 1)/dx
                    pdyvy = pdyw_stencil(vy_hxhyhz, i + 1, j, k + 1)/dy
                    pdzvy = pdzw_stencil(vy_hxiyiz, i + 1, j, k)/dz*eta_dz_scaling_h(k + 1)
                    pdxvz = pdxw_stencil(vz_ixiyhz, i, j, k + 1)/dx
                    pdyvz = pdyw_stencil(vz_hxhyhz, i + 1, j, k + 1)/dy
                    pdzvz = pdzw_stencil(vz_hxiyiz, i + 1, j, k)/dz*eta_dz_scaling_h(k + 1)

                    ! MPML
                    if (i + 1 <= 1 .or. i + 1 >= nx + 1 &
                            .or. j <= 0 .or. j >= ny + 1 &
                            .or. k + 1 >= nz + 1) then

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
                        if (k + 1 >= nz + 1) then
                            ratioz = dampratio_bottom(i + 1, j)
                        end if

                        dax = hdax(i + 1, j, k + 1) + ratioy*iday(i + 1, j, k + 1) + ratioz*hdaz(i + 1, j, k + 1)
                        a = dax/kappaxh(i + 1, j, k + 1)
                        b = dax/kappaxh(i + 1, j, k + 1) + alphaxh(i + 1, j, k + 1)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxh(i + 1, j, k + 1)
                        damp_wavefield(memory_pdxvx_hxiyhz, pdxvx, ax, bx, kx, i + 1, j, k + 1)
                        damp_wavefield(memory_pdxvy_hxiyhz, pdxvy, ax, bx, kx, i + 1, j, k + 1)
                        damp_wavefield(memory_pdxvz_hxiyhz, pdxvz, ax, bx, kx, i + 1, j, k + 1)

                        day = iday(i + 1, j, k + 1) + ratiox*hdax(i + 1, j, k + 1) + ratioz*hdaz(i + 1, j, k + 1)
                        a = day/kappayi(i + 1, j, k + 1)
                        b = day/kappayi(i + 1, j, k + 1) + alphayi(i + 1, j, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i + 1, j, k + 1)
                        damp_wavefield(memory_pdyvx_hxiyhz, pdyvx, ay, by, ky, i + 1, j, k + 1)
                        damp_wavefield(memory_pdyvy_hxiyhz, pdyvy, ay, by, ky, i + 1, j, k + 1)
                        damp_wavefield(memory_pdyvz_hxiyhz, pdyvz, ay, by, ky, i + 1, j, k + 1)

                        daz = hdaz(i + 1, j, k + 1) + ratiox*hdax(i + 1, j, k + 1) + ratioy*iday(i + 1, j, k + 1)
                        a = daz/kappazh(i + 1, j, k + 1)
                        b = daz/kappazh(i + 1, j, k + 1) + alphazh(i + 1, j, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i + 1, j, k + 1)
                        damp_wavefield(memory_pdzvx_hxiyhz, pdzvx, az, bz, kz, i + 1, j, k + 1)
                        damp_wavefield(memory_pdzvy_hxiyhz, pdzvy, az, bz, kz, i + 1, j, k + 1)
                        damp_wavefield(memory_pdzvz_hxiyhz, pdzvz, az, bz, kz, i + 1, j, k + 1)

                    end if

                    pd1 = pdxvx + tpaz*pdzvx
                    pd2 = pdyvy + tpbz*pdzvy
                    pd3 = tpcz*pdzvz
                    pd4 = tpcz*pdzvy + pdyvz + tpbz*pdzvz
                    pd5 = tpcz*pdzvx + pdxvz + tpaz*pdzvz
                    pd6 = pdyvx + tpbz*pdzvx + pdxvy + tpaz*pdzvy

                    stressxx_hxiyhz(i + 1, j, k + 1) = stressxx_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c11_hxiyhz*pd1 &
                        + c12_hxiyhz*pd2 &
                        + c13_hxiyhz*pd3 &
                        + c14_hxiyhz*pd4 &
                        + c15_hxiyhz*pd5 &
                        + c16_hxiyhz*pd6)

                    stressyy_hxiyhz(i + 1, j, k + 1) = stressyy_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c12_hxiyhz*pd1 &
                        + c22_hxiyhz*pd2 &
                        + c23_hxiyhz*pd3 &
                        + c24_hxiyhz*pd4 &
                        + c25_hxiyhz*pd5 &
                        + c26_hxiyhz*pd6)

                    stresszz_hxiyhz(i + 1, j, k + 1) = stresszz_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c13_hxiyhz*pd1 &
                        + c23_hxiyhz*pd2 &
                        + c33_hxiyhz*pd3 &
                        + c34_hxiyhz*pd4 &
                        + c35_hxiyhz*pd5 &
                        + c36_hxiyhz*pd6)

                    stressyz_hxiyhz(i + 1, j, k + 1) = stressyz_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c14_hxiyhz*pd1 &
                        + c24_hxiyhz*pd2 &
                        + c34_hxiyhz*pd3 &
                        + c44_hxiyhz*pd4 &
                        + c45_hxiyhz*pd5 &
                        + c46_hxiyhz*pd6)

                    stressxz_hxiyhz(i + 1, j, k + 1) = stressxz_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c15_hxiyhz*pd1 &
                        + c25_hxiyhz*pd2 &
                        + c35_hxiyhz*pd3 &
                        + c45_hxiyhz*pd4 &
                        + c55_hxiyhz*pd5 &
                        + c56_hxiyhz*pd6)

                    stressxy_hxiyhz(i + 1, j, k + 1) = stressxy_hxiyhz(i + 1, j, k + 1) + dt*( &
                        c16_hxiyhz*pd1 &
                        + c26_hxiyhz*pd2 &
                        + c36_hxiyhz*pd3 &
                        + c46_hxiyhz*pd4 &
                        + c56_hxiyhz*pd5 &
                        + c66_hxiyhz*pd6)

                end do
            end do
        end do
        !$omp end parallel do

        ! Stress set d (i, j+-1/2, k+-1/2)
        !$omp parallel do private(i, j, k, &
            !$omp pdxvx, pdyvx, pdzvx, pdxvy, pdyvy, pdzvy, pdxvz, pdyvz, pdzvz, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, tpaz, tpbz, tpcz, eta, &
            !$omp c11_ixhyhz, c12_ixhyhz, c13_ixhyhz, c14_ixhyhz, c15_ixhyhz, c16_ixhyhz, &
            !$omp c22_ixhyhz, c23_ixhyhz, c24_ixhyhz, c25_ixhyhz, c26_ixhyhz, &
            !$omp c33_ixhyhz, c34_ixhyhz, c35_ixhyhz, c36_ixhyhz, &
            !$omp c44_ixhyhz, c45_ixhyhz, c46_ixhyhz, &
            !$omp c55_ixhyhz, c56_ixhyhz, &
            !$omp c66_ixhyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1 - 1), nz2 - 1
            do j = ny1 - 1, ny2 - 1
                do i = nx1, nx2

                    eta = eta_zz_h(k + 1)
                    tpaz = max(0.0, eta_max - eta)/(topo_ixhy(i, j + 1) + depth_max)*slopex_ixhy(i, j + 1)
                    tpbz = max(0.0, eta_max - eta)/(topo_ixhy(i, j + 1) + depth_max)*slopey_ixhy(i, j + 1)
                    tpcz = eta_max/(topo_ixhy(i, j + 1) + depth_max)

                    c11_ixhyhz = 0.25*sum(c11(i, j:j + 1, k:k + 1))
                    c12_ixhyhz = 0.25*sum(c12(i, j:j + 1, k:k + 1))
                    c13_ixhyhz = 0.25*sum(c13(i, j:j + 1, k:k + 1))
                    c14_ixhyhz = 0.25*sum(c14(i, j:j + 1, k:k + 1))
                    c15_ixhyhz = 0.25*sum(c15(i, j:j + 1, k:k + 1))
                    c16_ixhyhz = 0.25*sum(c16(i, j:j + 1, k:k + 1))
                    c22_ixhyhz = 0.25*sum(c22(i, j:j + 1, k:k + 1))
                    c23_ixhyhz = 0.25*sum(c23(i, j:j + 1, k:k + 1))
                    c24_ixhyhz = 0.25*sum(c24(i, j:j + 1, k:k + 1))
                    c25_ixhyhz = 0.25*sum(c25(i, j:j + 1, k:k + 1))
                    c26_ixhyhz = 0.25*sum(c26(i, j:j + 1, k:k + 1))
                    c33_ixhyhz = 0.25*sum(c33(i, j:j + 1, k:k + 1))
                    c34_ixhyhz = 0.25*sum(c34(i, j:j + 1, k:k + 1))
                    c35_ixhyhz = 0.25*sum(c35(i, j:j + 1, k:k + 1))
                    c36_ixhyhz = 0.25*sum(c36(i, j:j + 1, k:k + 1))
                    c44_ixhyhz = 0.25*sum(c44(i, j:j + 1, k:k + 1))
                    c45_ixhyhz = 0.25*sum(c45(i, j:j + 1, k:k + 1))
                    c46_ixhyhz = 0.25*sum(c46(i, j:j + 1, k:k + 1))
                    c55_ixhyhz = 0.25*sum(c55(i, j:j + 1, k:k + 1))
                    c56_ixhyhz = 0.25*sum(c56(i, j:j + 1, k:k + 1))
                    c66_ixhyhz = 0.25*sum(c66(i, j:j + 1, k:k + 1))

                    pdxvx = pdxw_stencil(vx_hxhyhz, i, j + 1, k + 1)/dx
                    pdyvx = pdyw_stencil(vx_ixiyhz, i, j, k + 1)/dy
                    pdzvx = pdzw_stencil(vx_ixhyiz, i, j + 1, k)/dz*eta_dz_scaling_h(k + 1)
                    pdxvy = pdxw_stencil(vy_hxhyhz, i, j + 1, k + 1)/dx
                    pdyvy = pdyw_stencil(vy_ixiyhz, i, j, k + 1)/dy
                    pdzvy = pdzw_stencil(vy_ixhyiz, i, j + 1, k)/dz*eta_dz_scaling_h(k + 1)
                    pdxvz = pdxw_stencil(vz_hxhyhz, i, j + 1, k + 1)/dx
                    pdyvz = pdyw_stencil(vz_ixiyhz, i, j, k + 1)/dy
                    pdzvz = pdzw_stencil(vz_ixhyiz, i, j + 1, k)/dz*eta_dz_scaling_h(k + 1)

                    ! MPML
                    if (i <= 0 .or. i >= nx + 1 &
                            .or. j + 1 <= 1 .or. j + 1 >= ny + 1 &
                            .or. k + 1 >= nz + 1) then

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
                        if (k + 1 >= nz + 1) then
                            ratioz = dampratio_bottom(i, j + 1)
                        end if

                        dax = idax(i, j + 1, k + 1) + ratioy*hday(i, j + 1, k + 1) + ratioz*hdaz(i, j + 1, k + 1)
                        a = dax/kappaxi(i, j + 1, k + 1)
                        b = dax/kappaxi(i, j + 1, k + 1) + alphaxi(i, j + 1, k + 1)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxi(i, j + 1, k + 1)
                        damp_wavefield(memory_pdxvx_ixhyhz, pdxvx, ax, bx, kx, i, j + 1, k + 1)
                        damp_wavefield(memory_pdxvy_ixhyhz, pdxvy, ax, bx, kx, i, j + 1, k + 1)
                        damp_wavefield(memory_pdxvz_ixhyhz, pdxvz, ax, bx, kx, i, j + 1, k + 1)

                        day = hday(i, j + 1, k + 1) + ratiox*idax(i, j + 1, k + 1) + ratioz*hdaz(i, j + 1, k + 1)
                        a = day/kappayh(i, j + 1, k + 1)
                        b = day/kappayh(i, j + 1, k + 1) + alphayh(i, j + 1, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i, j + 1, k + 1)
                        damp_wavefield(memory_pdyvx_ixhyhz, pdyvx, ay, by, ky, i, j + 1, k + 1)
                        damp_wavefield(memory_pdyvy_ixhyhz, pdyvy, ay, by, ky, i, j + 1, k + 1)
                        damp_wavefield(memory_pdyvz_ixhyhz, pdyvz, ay, by, ky, i, j + 1, k + 1)

                        daz = hdaz(i, j + 1, k + 1) + ratiox*idax(i, j + 1, k + 1) + ratioy*hday(i, j + 1, k + 1)
                        a = daz/kappazh(i, j + 1, k + 1)
                        b = daz/kappazh(i, j + 1, k + 1) + alphazh(i, j + 1, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i, j + 1, k + 1)
                        damp_wavefield(memory_pdzvx_ixhyhz, pdzvx, az, bz, kz, i, j + 1, k + 1)
                        damp_wavefield(memory_pdzvy_ixhyhz, pdzvy, az, bz, kz, i, j + 1, k + 1)
                        damp_wavefield(memory_pdzvz_ixhyhz, pdzvz, az, bz, kz, i, j + 1, k + 1)

                    end if

                    pd1 = pdxvx + tpaz*pdzvx
                    pd2 = pdyvy + tpbz*pdzvy
                    pd3 = tpcz*pdzvz
                    pd4 = tpcz*pdzvy + pdyvz + tpbz*pdzvz
                    pd5 = tpcz*pdzvx + pdxvz + tpaz*pdzvz
                    pd6 = pdyvx + tpbz*pdzvx + pdxvy + tpaz*pdzvy

                    stressxx_ixhyhz(i, j + 1, k + 1) = stressxx_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c11_ixhyhz*pd1 &
                        + c12_ixhyhz*pd2 &
                        + c13_ixhyhz*pd3 &
                        + c14_ixhyhz*pd4 &
                        + c15_ixhyhz*pd5 &
                        + c16_ixhyhz*pd6)

                    stressyy_ixhyhz(i, j + 1, k + 1) = stressyy_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c12_ixhyhz*pd1 &
                        + c22_ixhyhz*pd2 &
                        + c23_ixhyhz*pd3 &
                        + c24_ixhyhz*pd4 &
                        + c25_ixhyhz*pd5 &
                        + c26_ixhyhz*pd6)

                    stresszz_ixhyhz(i, j + 1, k + 1) = stresszz_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c13_ixhyhz*pd1 &
                        + c23_ixhyhz*pd2 &
                        + c33_ixhyhz*pd3 &
                        + c34_ixhyhz*pd4 &
                        + c35_ixhyhz*pd5 &
                        + c36_ixhyhz*pd6)

                    stressyz_ixhyhz(i, j + 1, k + 1) = stressyz_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c14_ixhyhz*pd1 &
                        + c24_ixhyhz*pd2 &
                        + c34_ixhyhz*pd3 &
                        + c44_ixhyhz*pd4 &
                        + c45_ixhyhz*pd5 &
                        + c46_ixhyhz*pd6)

                    stressxz_ixhyhz(i, j + 1, k + 1) = stressxz_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c15_ixhyhz*pd1 &
                        + c25_ixhyhz*pd2 &
                        + c35_ixhyhz*pd3 &
                        + c45_ixhyhz*pd4 &
                        + c55_ixhyhz*pd5 &
                        + c56_ixhyhz*pd6)

                    stressxy_ixhyhz(i, j + 1, k + 1) = stressxy_ixhyhz(i, j + 1, k + 1) + dt*( &
                        c16_ixhyhz*pd1 &
                        + c26_ixhyhz*pd2 &
                        + c36_ixhyhz*pd3 &
                        + c46_ixhyhz*pd4 &
                        + c56_ixhyhz*pd5 &
                        + c66_ixhyhz*pd6)

                end do
            end do
        end do
        !$omp end parallel do

        ! Apply mirror boundary condition to mimic free surface
        ! It is possible that 1 - k or 2 - k are within the block, but their mirror points
        ! are not. In this case, the communication must come before the mirroring.
        !$omp parallel do private(k) schedule(auto)
        do k = 1, fdhalf

            if (1 - k >= nz1 .and. 1 - k <= nz2 .and. 1 + k >= nz1 - fdhalf + 1 .and. 1 + k <= nz2 + fdhalf) then
                stressxx_ixiyiz(:, :, 1 - k) = -stressxx_ixiyiz(:, :, 1 + k)
                stressyy_ixiyiz(:, :, 1 - k) = -stressyy_ixiyiz(:, :, 1 + k)
                stresszz_ixiyiz(:, :, 1 - k) = -stresszz_ixiyiz(:, :, 1 + k)
                stressyz_ixiyiz(:, :, 1 - k) = -stressyz_ixiyiz(:, :, 1 + k)
                stressxz_ixiyiz(:, :, 1 - k) = -stressxz_ixiyiz(:, :, 1 + k)
                stressxy_ixiyiz(:, :, 1 - k) = -stressxy_ixiyiz(:, :, 1 + k)

                stressxx_hxhyiz(:, :, 1 - k) = -stressxx_hxhyiz(:, :, 1 + k)
                stressyy_hxhyiz(:, :, 1 - k) = -stressyy_hxhyiz(:, :, 1 + k)
                stresszz_hxhyiz(:, :, 1 - k) = -stresszz_hxhyiz(:, :, 1 + k)
                stressyz_hxhyiz(:, :, 1 - k) = -stressyz_hxhyiz(:, :, 1 + k)
                stressxz_hxhyiz(:, :, 1 - k) = -stressxz_hxhyiz(:, :, 1 + k)
                stressxy_hxhyiz(:, :, 1 - k) = -stressxy_hxhyiz(:, :, 1 + k)
            end if

            if (2 - k >= nz1 .and. 2 - k <= nz2 .and. 1 + k >= nz1 - fdhalf + 1 .and. 1 + k <= nz2 + fdhalf) then
                stressxx_hxiyhz(:, :, 2 - k) = -stressxx_hxiyhz(:, :, 1 + k)
                stressyy_hxiyhz(:, :, 2 - k) = -stressyy_hxiyhz(:, :, 1 + k)
                stresszz_hxiyhz(:, :, 2 - k) = -stresszz_hxiyhz(:, :, 1 + k)
                stressyz_hxiyhz(:, :, 2 - k) = -stressyz_hxiyhz(:, :, 1 + k)
                stressxz_hxiyhz(:, :, 2 - k) = -stressxz_hxiyhz(:, :, 1 + k)
                stressxy_hxiyhz(:, :, 2 - k) = -stressxy_hxiyhz(:, :, 1 + k)

                stressxx_ixhyhz(:, :, 2 - k) = -stressxx_ixhyhz(:, :, 1 + k)
                stressyy_ixhyhz(:, :, 2 - k) = -stressyy_ixhyhz(:, :, 1 + k)
                stresszz_ixhyhz(:, :, 2 - k) = -stresszz_ixhyhz(:, :, 1 + k)
                stressyz_ixhyhz(:, :, 2 - k) = -stressyz_ixhyhz(:, :, 1 + k)
                stressxz_ixhyhz(:, :, 2 - k) = -stressxz_ixhyhz(:, :, 1 + k)
                stressxy_ixhyhz(:, :, 2 - k) = -stressxy_ixhyhz(:, :, 1 + k)
            end if

        end do
        !$omp end parallel do

        ! Exchange boundary wavefields
        call commute_array_group(stressxx_ixiyiz, fdhalf)
        call commute_array_group(stressyy_ixiyiz, fdhalf)
        call commute_array_group(stresszz_ixiyiz, fdhalf)
        call commute_array_group(stressxy_ixiyiz, fdhalf)
        call commute_array_group(stressxz_ixiyiz, fdhalf)
        call commute_array_group(stressyz_ixiyiz, fdhalf)
        call commute_array_group(stressxx_hxhyiz, fdhalf)
        call commute_array_group(stressyy_hxhyiz, fdhalf)
        call commute_array_group(stresszz_hxhyiz, fdhalf)
        call commute_array_group(stressxy_hxhyiz, fdhalf)
        call commute_array_group(stressxz_hxhyiz, fdhalf)
        call commute_array_group(stressyz_hxhyiz, fdhalf)
        call commute_array_group(stressxx_hxiyhz, fdhalf)
        call commute_array_group(stressyy_hxiyhz, fdhalf)
        call commute_array_group(stresszz_hxiyhz, fdhalf)
        call commute_array_group(stressxy_hxiyhz, fdhalf)
        call commute_array_group(stressxz_hxiyhz, fdhalf)
        call commute_array_group(stressyz_hxiyhz, fdhalf)
        call commute_array_group(stressxx_ixhyhz, fdhalf)
        call commute_array_group(stressyy_ixhyhz, fdhalf)
        call commute_array_group(stresszz_ixhyhz, fdhalf)
        call commute_array_group(stressxy_ixhyhz, fdhalf)
        call commute_array_group(stressxz_ixhyhz, fdhalf)
        call commute_array_group(stressyz_ixhyhz, fdhalf)

        ! Particle velocity set a (i+-1/2, j, k)
        !$omp parallel do private(i, j, k, &
            !$omp pdxxx, pdyxy, pdzxz, &
            !$omp pdxxy, pdyyy, pdzyz, &
            !$omp pdxxz, pdyyz, pdzzz, &
            !$omp pdzxx, pdzxy, pdzyy, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, tpaz, tpbz, tpcz, eta, &
            !$omp rho_hxiyiz, rho_ixhyiz, rho_ixiyhz, rho_hxhyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1), nz2
            do j = ny1, ny2
                do i = nx1 - 1, nx2 - 1

                    eta = eta_zz_i(k)
                    tpaz = max(0.0, eta_max - eta)/(topo_hxiy(i + 1, j) + depth_max)*slopex_hxiy(i + 1, j)
                    tpbz = max(0.0, eta_max - eta)/(topo_hxiy(i + 1, j) + depth_max)*slopey_hxiy(i + 1, j)
                    tpcz = eta_max/(topo_hxiy(i + 1, j) + depth_max)

                    rho_hxiyiz = 0.5*sum(rho(i:i + 1, j, k))

                    pdxxx = pdxw_stencil(stressxx_ixiyiz, i, j, k)/dx
                    pdxxy = pdxw_stencil(stressxy_ixiyiz, i, j, k)/dx
                    pdxxz = pdxw_stencil(stressxz_ixiyiz, i, j, k)/dx
                    pdyxy = pdyw_stencil(stressxy_hxhyiz, i + 1, j, k)/dy
                    pdyyy = pdyw_stencil(stressyy_hxhyiz, i + 1, j, k)/dy
                    pdyyz = pdyw_stencil(stressyz_hxhyiz, i + 1, j, k)/dy

                    pdzxx = pdzw_stencil(stressxx_hxiyhz, i + 1, j, k)/dz*eta_dz_scaling_i(k)
                    pdzxy = pdzw_stencil(stressxy_hxiyhz, i + 1, j, k)/dz*eta_dz_scaling_i(k)
                    pdzxz = pdzw_stencil(stressxz_hxiyhz, i + 1, j, k)/dz*eta_dz_scaling_i(k)
                    pdzyy = pdzw_stencil(stressyy_hxiyhz, i + 1, j, k)/dz*eta_dz_scaling_i(k)
                    pdzyz = pdzw_stencil(stressyz_hxiyhz, i + 1, j, k)/dz*eta_dz_scaling_i(k)
                    pdzzz = pdzw_stencil(stresszz_hxiyhz, i + 1, j, k)/dz*eta_dz_scaling_i(k)

                    ! MPML
                    if (i + 1 <= 1 .or. i + 1 >= nx + 1 &
                            .or. j <= 0 .or. j >= ny + 1 &
                            .or. k >= nz + 1) then

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
                        if (k >= nz + 1) then
                            ratioz = dampratio_bottom(i + 1, j)
                        end if

                        dax = hdax(i + 1, j, k) + ratioy*iday(i + 1, j, k) + ratioz*idaz(i + 1, j, k)
                        a = dax/kappaxh(i + 1, j, k)
                        b = dax/kappaxh(i + 1, j, k) + alphaxh(i + 1, j, k)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxh(i + 1, j, k)
                        damp_wavefield(memory_pdxxx_hxiyiz, pdxxx, ax, bx, kx, i + 1, j, k)
                        damp_wavefield(memory_pdxxy_hxiyiz, pdxxy, ax, bx, kx, i + 1, j, k)
                        damp_wavefield(memory_pdxxz_hxiyiz, pdxxz, ax, bx, kx, i + 1, j, k)

                        day = iday(i + 1, j, k) + ratiox*hdax(i + 1, j, k) + ratioz*idaz(i + 1, j, k)
                        a = day/kappayi(i + 1, j, k)
                        b = day/kappayi(i + 1, j, k) + alphayi(i + 1, j, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i + 1, j, k)
                        damp_wavefield(memory_pdyxy_hxiyiz, pdyxy, ay, by, ky, i + 1, j, k)
                        damp_wavefield(memory_pdyyy_hxiyiz, pdyyy, ay, by, ky, i + 1, j, k)
                        damp_wavefield(memory_pdyyz_hxiyiz, pdyyz, ay, by, ky, i + 1, j, k)

                        daz = idaz(i + 1, j, k) + ratiox*hdax(i + 1, j, k) + ratioy*iday(i + 1, j, k)
                        a = daz/kappazi(i + 1, j, k)
                        b = daz/kappazi(i + 1, j, k) + alphazi(i + 1, j, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i + 1, j, k)
                        damp_wavefield(memory_pdzxx_hxiyiz, pdzxx, az, bz, kz, i + 1, j, k)
                        damp_wavefield(memory_pdzxy_hxiyiz, pdzxy, az, bz, kz, i + 1, j, k)
                        damp_wavefield(memory_pdzxz_hxiyiz, pdzxz, az, bz, kz, i + 1, j, k)
                        damp_wavefield(memory_pdzyy_hxiyiz, pdzyy, az, bz, kz, i + 1, j, k)
                        damp_wavefield(memory_pdzyz_hxiyiz, pdzyz, az, bz, kz, i + 1, j, k)
                        damp_wavefield(memory_pdzzz_hxiyiz, pdzzz, az, bz, kz, i + 1, j, k)

                    end if

                    pd1 = pdxxx + tpaz*pdzxx &
                        + pdyxy + tpbz*pdzxy &
                        + tpcz*pdzxz
                    pd2 = pdxxy + tpaz*pdzxy &
                        + pdyyy + tpbz*pdzyy &
                        + tpcz*pdzyz
                    pd3 = pdxxz + tpaz*pdzxz &
                        + pdyyz + tpbz*pdzyz &
                        + tpcz*pdzzz

                    vx_hxiyiz(i + 1, j, k) = vx_hxiyiz(i + 1, j, k) + dt/rho_hxiyiz*pd1
                    vy_hxiyiz(i + 1, j, k) = vy_hxiyiz(i + 1, j, k) + dt/rho_hxiyiz*pd2
                    vz_hxiyiz(i + 1, j, k) = vz_hxiyiz(i + 1, j, k) + dt/rho_hxiyiz*pd3

                end do
            end do
        end do
        !$omp end parallel do

        ! Particle velocity set b (i, j+-1/2, k)
        !$omp parallel do private(i, j, k, &
            !$omp pdxxx, pdyxy, pdzxz, &
            !$omp pdxxy, pdyyy, pdzyz, &
            !$omp pdxxz, pdyyz, pdzzz, &
            !$omp pdzxx, pdzxy, pdzyy, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, tpaz, tpbz, tpcz, eta, &
            !$omp rho_hxiyiz, rho_ixhyiz, rho_ixiyhz, rho_hxhyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1), nz2
            do j = ny1 - 1, ny2 - 1
                do i = nx1, nx2

                    eta = eta_zz_i(k)
                    tpaz = max(0.0, eta_max - eta)/(topo_ixhy(i, j + 1) + depth_max)*slopex_ixhy(i, j + 1)
                    tpbz = max(0.0, eta_max - eta)/(topo_ixhy(i, j + 1) + depth_max)*slopey_ixhy(i, j + 1)
                    tpcz = eta_max/(topo_ixhy(i, j + 1) + depth_max)

                    rho_ixhyiz = 0.5*sum(rho(i, j:j + 1, k))

                    pdxxx = pdxw_stencil(stressxx_hxhyiz, i, j + 1, k)/dx
                    pdxxy = pdxw_stencil(stressxy_hxhyiz, i, j + 1, k)/dx
                    pdxxz = pdxw_stencil(stressxz_hxhyiz, i, j + 1, k)/dx
                    pdyxy = pdyw_stencil(stressxy_ixiyiz, i, j, k)/dy
                    pdyyy = pdyw_stencil(stressyy_ixiyiz, i, j, k)/dy
                    pdyyz = pdyw_stencil(stressyz_ixiyiz, i, j, k)/dy

                    pdzxx = pdzw_stencil(stressxx_ixhyhz, i, j + 1, k)/dz*eta_dz_scaling_i(k)
                    pdzxy = pdzw_stencil(stressxy_ixhyhz, i, j + 1, k)/dz*eta_dz_scaling_i(k)
                    pdzxz = pdzw_stencil(stressxz_ixhyhz, i, j + 1, k)/dz*eta_dz_scaling_i(k)
                    pdzyy = pdzw_stencil(stressyy_ixhyhz, i, j + 1, k)/dz*eta_dz_scaling_i(k)
                    pdzyz = pdzw_stencil(stressyz_ixhyhz, i, j + 1, k)/dz*eta_dz_scaling_i(k)
                    pdzzz = pdzw_stencil(stresszz_ixhyhz, i, j + 1, k)/dz*eta_dz_scaling_i(k)

                    ! MPML
                    if (i <= 0 .or. i >= nx + 1 &
                            .or. j + 1 <= 1 .or. j + 1 >= ny + 1 &
                            .or. k >= nz + 1) then

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
                        if (k >= nz + 1) then
                            ratioz = dampratio_bottom(i, j + 1)
                        end if

                        dax = idax(i, j + 1, k) + ratioy*hday(i, j + 1, k) + ratioz*idaz(i, j + 1, k)
                        a = dax/kappaxi(i, j + 1, k)
                        b = dax/kappaxi(i, j + 1, k) + alphaxi(i, j + 1, k)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxi(i, j + 1, k)
                        damp_wavefield(memory_pdxxx_ixhyiz, pdxxx, ax, bx, kx, i, j + 1, k)
                        damp_wavefield(memory_pdxxy_ixhyiz, pdxxy, ax, bx, kx, i, j + 1, k)
                        damp_wavefield(memory_pdxxz_ixhyiz, pdxxz, ax, bx, kx, i, j + 1, k)

                        day = hday(i, j + 1, k) + ratiox*idax(i, j + 1, k) + ratioz*idaz(i, j + 1, k)
                        a = day/kappayh(i, j + 1, k)
                        b = day/kappayh(i, j + 1, k) + alphayh(i, j + 1, k)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i, j + 1, k)
                        damp_wavefield(memory_pdyxy_ixhyiz, pdyxy, ay, by, ky, i, j + 1, k)
                        damp_wavefield(memory_pdyyy_ixhyiz, pdyyy, ay, by, ky, i, j + 1, k)
                        damp_wavefield(memory_pdyyz_ixhyiz, pdyyz, ay, by, ky, i, j + 1, k)

                        daz = idaz(i, j + 1, k) + ratiox*idax(i, j + 1, k) + ratioy*hday(i, j + 1, k)
                        a = daz/kappazi(i, j + 1, k)
                        b = daz/kappazi(i, j + 1, k) + alphazi(i, j + 1, k)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazi(i, j + 1, k)
                        damp_wavefield(memory_pdzxx_ixhyiz, pdzxx, az, bz, kz, i, j + 1, k)
                        damp_wavefield(memory_pdzxy_ixhyiz, pdzxy, az, bz, kz, i, j + 1, k)
                        damp_wavefield(memory_pdzxz_ixhyiz, pdzxz, az, bz, kz, i, j + 1, k)
                        damp_wavefield(memory_pdzyy_ixhyiz, pdzyy, az, bz, kz, i, j + 1, k)
                        damp_wavefield(memory_pdzyz_ixhyiz, pdzyz, az, bz, kz, i, j + 1, k)
                        damp_wavefield(memory_pdzzz_ixhyiz, pdzzz, az, bz, kz, i, j + 1, k)

                    end if

                    pd1 = pdxxx + tpaz*pdzxx &
                        + pdyxy + tpbz*pdzxy &
                        + tpcz*pdzxz
                    pd2 = pdxxy + tpaz*pdzxy &
                        + pdyyy + tpbz*pdzyy &
                        + tpcz*pdzyz
                    pd3 = pdxxz + tpaz*pdzxz &
                        + pdyyz + tpbz*pdzyz &
                        + tpcz*pdzzz

                    vx_ixhyiz(i, j + 1, k) = vx_ixhyiz(i, j + 1, k) + dt/rho_ixhyiz*pd1
                    vy_ixhyiz(i, j + 1, k) = vy_ixhyiz(i, j + 1, k) + dt/rho_ixhyiz*pd2
                    vz_ixhyiz(i, j + 1, k) = vz_ixhyiz(i, j + 1, k) + dt/rho_ixhyiz*pd3

                end do
            end do
        end do
        !$omp end parallel do

        ! Particle velocity set c (i, j, k+-1/2)
        !$omp parallel do private(i, j, k, &
            !$omp pdxxx, pdyxy, pdzxz, &
            !$omp pdxxy, pdyyy, pdzyz, &
            !$omp pdxxz, pdyyz, pdzzz, &
            !$omp pdzxx, pdzxy, pdzyy, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, tpaz, tpbz, tpcz, eta, &
            !$omp rho_hxiyiz, rho_ixhyiz, rho_ixiyhz, rho_hxhyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1 - 1), nz2 - 1
            do j = ny1, ny2
                do i = nx1, nx2

                    eta = eta_zz_h(k + 1)
                    tpaz = max(0.0, eta_max - eta)/(topo_ixiy(i, j) + depth_max)*slopex_ixiy(i, j)
                    tpbz = max(0.0, eta_max - eta)/(topo_ixiy(i, j) + depth_max)*slopey_ixiy(i, j)
                    tpcz = eta_max/(topo_ixiy(i, j) + depth_max)

                    rho_ixiyhz = 0.5*sum(rho(i, j, k:k + 1))

                    pdxxx = pdxw_stencil(stressxx_hxiyhz, i, j, k + 1)/dx
                    pdxxy = pdxw_stencil(stressxy_hxiyhz, i, j, k + 1)/dx
                    pdxxz = pdxw_stencil(stressxz_hxiyhz, i, j, k + 1)/dx
                    pdyxy = pdyw_stencil(stressxy_ixhyhz, i, j, k + 1)/dy
                    pdyyy = pdyw_stencil(stressyy_ixhyhz, i, j, k + 1)/dy
                    pdyyz = pdyw_stencil(stressyz_ixhyhz, i, j, k + 1)/dy

                    pdzxx = pdzw_stencil(stressxx_ixiyiz, i, j, k)/dz*eta_dz_scaling_h(k + 1)
                    pdzxy = pdzw_stencil(stressxy_ixiyiz, i, j, k)/dz*eta_dz_scaling_h(k + 1)
                    pdzxz = pdzw_stencil(stressxz_ixiyiz, i, j, k)/dz*eta_dz_scaling_h(k + 1)
                    pdzyy = pdzw_stencil(stressyy_ixiyiz, i, j, k)/dz*eta_dz_scaling_h(k + 1)
                    pdzyz = pdzw_stencil(stressyz_ixiyiz, i, j, k)/dz*eta_dz_scaling_h(k + 1)
                    pdzzz = pdzw_stencil(stresszz_ixiyiz, i, j, k)/dz*eta_dz_scaling_h(k + 1)

                    ! MPML
                    if (i <= 0 .or. i >= nx + 1 &
                            .or. j <= 0 .or. j >= ny + 1 &
                            .or. k + 1 >= nz + 1) then

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
                        if (k + 1 >= nz + 1) then
                            ratioz = dampratio_bottom(i, j)
                        end if

                        dax = idax(i, j, k + 1) + ratioy*iday(i, j, k + 1) + ratioz*hdaz(i, j, k + 1)
                        a = dax/kappaxi(i, j, k + 1)
                        b = dax/kappaxi(i, j, k + 1) + alphaxi(i, j, k + 1)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxi(i, j, k + 1)
                        damp_wavefield(memory_pdxxx_ixiyhz, pdxxx, ax, bx, kx, i, j, k + 1)
                        damp_wavefield(memory_pdxxy_ixiyhz, pdxxy, ax, bx, kx, i, j, k + 1)
                        damp_wavefield(memory_pdxxz_ixiyhz, pdxxz, ax, bx, kx, i, j, k + 1)

                        day = iday(i, j, k + 1) + ratiox*idax(i, j, k + 1) + ratioz*hdaz(i, j, k + 1)
                        a = day/kappayi(i, j, k + 1)
                        b = day/kappayi(i, j, k + 1) + alphayi(i, j, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayi(i, j, k + 1)
                        damp_wavefield(memory_pdyxy_ixiyhz, pdyxy, ay, by, ky, i, j, k + 1)
                        damp_wavefield(memory_pdyyy_ixiyhz, pdyyy, ay, by, ky, i, j, k + 1)
                        damp_wavefield(memory_pdyyz_ixiyhz, pdyyz, ay, by, ky, i, j, k + 1)

                        daz = hdaz(i, j, k + 1) + ratiox*idax(i, j, k + 1) + ratioy*iday(i, j, k + 1)
                        a = daz/kappazh(i, j, k + 1)
                        b = daz/kappazh(i, j, k + 1) + alphazh(i, j, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i, j, k + 1)
                        damp_wavefield(memory_pdzxx_ixiyhz, pdzxx, az, bz, kz, i, j, k + 1)
                        damp_wavefield(memory_pdzxy_ixiyhz, pdzxy, az, bz, kz, i, j, k + 1)
                        damp_wavefield(memory_pdzxz_ixiyhz, pdzxz, az, bz, kz, i, j, k + 1)
                        damp_wavefield(memory_pdzyy_ixiyhz, pdzyy, az, bz, kz, i, j, k + 1)
                        damp_wavefield(memory_pdzyz_ixiyhz, pdzyz, az, bz, kz, i, j, k + 1)
                        damp_wavefield(memory_pdzzz_ixiyhz, pdzzz, az, bz, kz, i, j, k + 1)

                    end if

                    pd1 = pdxxx + tpaz*pdzxx &
                        + pdyxy + tpbz*pdzxy &
                        + tpcz*pdzxz
                    pd2 = pdxxy + tpaz*pdzxy &
                        + pdyyy + tpbz*pdzyy &
                        + tpcz*pdzyz
                    pd3 = pdxxz + tpaz*pdzxz &
                        + pdyyz + tpbz*pdzyz &
                        + tpcz*pdzzz

                    vx_ixiyhz(i, j, k + 1) = vx_ixiyhz(i, j, k + 1) + dt/rho_ixiyhz*pd1
                    vy_ixiyhz(i, j, k + 1) = vy_ixiyhz(i, j, k + 1) + dt/rho_ixiyhz*pd2
                    vz_ixiyhz(i, j, k + 1) = vz_ixiyhz(i, j, k + 1) + dt/rho_ixiyhz*pd3

                end do
            end do
        end do
        !$omp end parallel do

        ! Particle velocity set d (i+-1/2, j+-1/2, k+-1/2)
        !$omp parallel do private(i, j, k, &
            !$omp pdxxx, pdyxy, pdzxz, &
            !$omp pdxxy, pdyyy, pdzyz, &
            !$omp pdxxz, pdyyz, pdzzz, &
            !$omp pdzxx, pdzxy, pdzyy, &
            !$omp pd1, pd2, pd3, pd4, pd5, pd6, tpaz, tpbz, tpcz, eta, &
            !$omp rho_hxiyiz, rho_ixhyiz, rho_ixiyhz, rho_hxhyhz, &
            !$omp ratiox, ratioy, ratioz, dax, day, daz, &
            !$omp a, b, ax, ay, az, bx, by, bz, kx, ky, kz) collapse(3) schedule(auto)
        do k = max(1, nz1 - 1), nz2 - 1
            do j = ny1 - 1, ny2 - 1
                do i = nx1 - 1, nx2 - 1

                    eta = eta_zz_h(k + 1)
                    tpaz = max(0.0, eta_max - eta)/(topo_hxhy(i + 1, j + 1) + depth_max)*slopex_hxhy(i + 1, j + 1)
                    tpbz = max(0.0, eta_max - eta)/(topo_hxhy(i + 1, j + 1) + depth_max)*slopey_hxhy(i + 1, j + 1)
                    tpcz = eta_max/(topo_hxhy(i + 1, j + 1) + depth_max)

                    rho_hxhyhz = 0.125*sum(rho(i:i + 1, j:j + 1, k:k + 1))

                    pdxxx = pdxw_stencil(stressxx_ixhyhz, i, j + 1, k + 1)/dx
                    pdxxy = pdxw_stencil(stressxy_ixhyhz, i, j + 1, k + 1)/dx
                    pdxxz = pdxw_stencil(stressxz_ixhyhz, i, j + 1, k + 1)/dx
                    pdyxy = pdyw_stencil(stressxy_hxiyhz, i + 1, j, k + 1)/dy
                    pdyyy = pdyw_stencil(stressyy_hxiyhz, i + 1, j, k + 1)/dy
                    pdyyz = pdyw_stencil(stressyz_hxiyhz, i + 1, j, k + 1)/dy

                    pdzxx = pdzw_stencil(stressxx_hxhyiz, i + 1, j + 1, k)/dz*eta_dz_scaling_h(k + 1)
                    pdzxy = pdzw_stencil(stressxy_hxhyiz, i + 1, j + 1, k)/dz*eta_dz_scaling_h(k + 1)
                    pdzxz = pdzw_stencil(stressxz_hxhyiz, i + 1, j + 1, k)/dz*eta_dz_scaling_h(k + 1)
                    pdzyy = pdzw_stencil(stressyy_hxhyiz, i + 1, j + 1, k)/dz*eta_dz_scaling_h(k + 1)
                    pdzyz = pdzw_stencil(stressyz_hxhyiz, i + 1, j + 1, k)/dz*eta_dz_scaling_h(k + 1)
                    pdzzz = pdzw_stencil(stresszz_hxhyiz, i + 1, j + 1, k)/dz*eta_dz_scaling_h(k + 1)

                    ! MPML
                    if (i + 1 <= 1 .or. i + 1 >= nx + 1 &
                            .or. j + 1 <= 1 .or. j + 1 >= ny + 1 &
                            .or. k + 1 >= nz + 1) then

                        ratiox = 0.0
                        ratioy = 0.0
                        ratioz = 0.0

                        if (i + 1 <= 1) then
                            ratiox = dampratio_left(j + 1, k + 1)
                        else if (i + 1 >= nx + 1) then
                            ratiox = dampratio_right(j + 1, k + 1)
                        end if
                        if (j + 1 <= 1) then
                            ratioy = dampratio_front(i + 1, k + 1)
                        else if (j + 1 >= ny + 1) then
                            ratioy = dampratio_back(i + 1, k + 1)
                        end if
                        if (k + 1 >= nz + 1) then
                            ratioz = dampratio_bottom(i + 1, j + 1)
                        end if

                        dax = hdax(i + 1, j + 1, k + 1) + ratioy*hday(i + 1, j + 1, k + 1) + ratioz*hdaz(i + 1, j + 1, k + 1)
                        a = dax/kappaxh(i + 1, j + 1, k + 1)
                        b = dax/kappaxh(i + 1, j + 1, k + 1) + alphaxh(i + 1, j + 1, k + 1)
                        damp_coef(dt, a, b, ax, bx)
                        kx = kappaxh(i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdxxx_hxhyhz, pdxxx, ax, bx, kx, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdxxy_hxhyhz, pdxxy, ax, bx, kx, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdxxz_hxhyhz, pdxxz, ax, bx, kx, i + 1, j + 1, k + 1)

                        day = hday(i + 1, j + 1, k + 1) + ratiox*hdax(i + 1, j + 1, k + 1) + ratioz*hdaz(i + 1, j + 1, k + 1)
                        a = day/kappayh(i + 1, j + 1, k + 1)
                        b = day/kappayh(i + 1, j + 1, k + 1) + alphayh(i + 1, j + 1, k + 1)
                        damp_coef(dt, a, b, ay, by)
                        ky = kappayh(i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdyxy_hxhyhz, pdyxy, ay, by, ky, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdyyy_hxhyhz, pdyyy, ay, by, ky, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdyyz_hxhyhz, pdyyz, ay, by, ky, i + 1, j + 1, k + 1)

                        daz = hdaz(i + 1, j + 1, k + 1) + ratiox*hdax(i + 1, j + 1, k + 1) + ratioy*hday(i + 1, j + 1, k + 1)
                        a = daz/kappazh(i + 1, j + 1, k + 1)
                        b = daz/kappazh(i + 1, j + 1, k + 1) + alphazh(i + 1, j + 1, k + 1)
                        damp_coef(dt, a, b, az, bz)
                        kz = kappazh(i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdzxx_hxhyhz, pdzxx, az, bz, kz, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdzxy_hxhyhz, pdzxy, az, bz, kz, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdzxz_hxhyhz, pdzxz, az, bz, kz, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdzyy_hxhyhz, pdzyy, az, bz, kz, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdzyz_hxhyhz, pdzyz, az, bz, kz, i + 1, j + 1, k + 1)
                        damp_wavefield(memory_pdzzz_hxhyhz, pdzzz, az, bz, kz, i + 1, j + 1, k + 1)

                    end if

                    pd1 = pdxxx + tpaz*pdzxx &
                        + pdyxy + tpbz*pdzxy &
                        + tpcz*pdzxz
                    pd2 = pdxxy + tpaz*pdzxy &
                        + pdyyy + tpbz*pdzyy &
                        + tpcz*pdzyz
                    pd3 = pdxxz + tpaz*pdzxz &
                        + pdyyz + tpbz*pdzyz &
                        + tpcz*pdzzz

                    vx_hxhyhz(i + 1, j + 1, k + 1) = vx_hxhyhz(i + 1, j + 1, k + 1) + dt/rho_hxhyhz*pd1
                    vy_hxhyhz(i + 1, j + 1, k + 1) = vy_hxhyhz(i + 1, j + 1, k + 1) + dt/rho_hxhyhz*pd2
                    vz_hxhyhz(i + 1, j + 1, k + 1) = vz_hxhyhz(i + 1, j + 1, k + 1) + dt/rho_hxhyhz*pd3

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
                dz_s = eta_dz_i(sgmtr%srcr(k)%gz)*(topo_max + depth_max)/eta_max
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
                                    ! hx-iy-iz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        vx_hxiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vx_hxiyiz(sgx + irx, sgy + iry, sgz + irz) + sin(polar)*cos(azimuth)*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        vy_hxiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vy_hxiyiz(sgx + irx, sgy + iry, sgz + irz) + sin(polar)*sin(azimuth)*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        vz_hxiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vz_hxiyiz(sgx + irx, sgy + iry, sgz + irz) + cos(polar)*amp &
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
                                    ! ix-hy-iz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        vx_ixhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vx_ixhyiz(sgx + irx, sgy + iry, sgz + irz) + sin(polar)*cos(azimuth)*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        vy_ixhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vy_ixhyiz(sgx + irx, sgy + iry, sgz + irz) + sin(polar)*sin(azimuth)*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        vz_ixhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vz_ixhyiz(sgx + irx, sgy + iry, sgz + irz) + cos(polar)*amp &
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
                                    ! ix-iy-hz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        vx_ixiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vx_ixiyhz(sgx + irx, sgy + iry, sgz + irz) + sin(polar)*cos(azimuth)*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        vy_ixiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vy_ixiyhz(sgx + irx, sgy + iry, sgz + irz) + sin(polar)*sin(azimuth)*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        vz_ixiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vz_ixiyhz(sgx + irx, sgy + iry, sgz + irz) + cos(polar)*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                        sgx = sgmtr%srcr(k)%hx
                        sgy = sgmtr%srcr(k)%hy
                        sgz = sgmtr%srcr(k)%hz
                        !$omp parallel do private(irx, iry, irz) collapse(3) schedule(auto)
                        do irz = -nkw, nkw
                            do iry = -nkw, nkw
                                do irx = -nkw, nkw
                                    ! hx-hy-hz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        vx_hxhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vx_hxhyhz(sgx + irx, sgy + iry, sgz + irz) + sin(polar)*cos(azimuth)*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        vy_hxhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vy_hxhyhz(sgx + irx, sgy + iry, sgz + irz) + sin(polar)*sin(azimuth)*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        vz_hxhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            vz_hxhyhz(sgx + irx, sgy + iry, sgz + irz) + cos(polar)*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
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
                                    ! ix-iy-iz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxx_ixiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxx_ixiyiz(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressyy_ixiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyy_ixiyiz(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stresszz_ixiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stresszz_ixiyiz(sgx + irx, sgy + iry, sgz + irz) - amp &
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
                                    ! hx-hy-iz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxx_hxhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxx_hxhyiz(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressyy_hxhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyy_hxhyiz(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stresszz_hxhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stresszz_hxhyiz(sgx + irx, sgy + iry, sgz + irz) - amp &
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
                                    ! hx-iy-hz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxx_hxiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxx_hxiyhz(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stressyy_hxiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyy_hxiyhz(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stresszz_hxiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stresszz_hxiyhz(sgx + irx, sgy + iry, sgz + irz) - amp &
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
                                    ! ix-hy-hz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxx_ixhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxx_ixhyhz(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stressyy_ixhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyy_ixhyhz(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stresszz_ixhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stresszz_ixhyhz(sgx + irx, sgy + iry, sgz + irz) - amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                    end if
                                end do
                            end do
                        end do
                        !$omp end parallel do

                    case ('mt')
                        ! Moment tensor
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
                                    ! ix-iy-iz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxx_ixiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxx_ixiyiz(sgx + irx, sgy + iry, sgz + irz) - m11*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressyy_ixiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyy_ixiyiz(sgx + irx, sgy + iry, sgz + irz) - m22*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stresszz_ixiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stresszz_ixiyiz(sgx + irx, sgy + iry, sgz + irz) - m33*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressxy_ixiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxy_ixiyiz(sgx + irx, sgy + iry, sgz + irz) - m12*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressxz_ixiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxz_ixiyiz(sgx + irx, sgy + iry, sgz + irz) - m13*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressyz_ixiyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyz_ixiyiz(sgx + irx, sgy + iry, sgz + irz) - m23*amp &
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
                                    ! hx-hy-iz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxx_hxhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxx_hxhyiz(sgx + irx, sgy + iry, sgz + irz) - m11*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressyy_hxhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyy_hxhyiz(sgx + irx, sgy + iry, sgz + irz) - m22*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stresszz_hxhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stresszz_hxhyiz(sgx + irx, sgy + iry, sgz + irz) - m33*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressxy_hxhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxy_hxhyiz(sgx + irx, sgy + iry, sgz + irz) - m12*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressxz_hxhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxz_hxhyiz(sgx + irx, sgy + iry, sgz + irz) - m13*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_iz(irz)
                                        stressyz_hxhyiz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyz_hxhyiz(sgx + irx, sgy + iry, sgz + irz) - m23*amp &
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
                                    ! hx-iy-hz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxx_hxiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxx_hxiyhz(sgx + irx, sgy + iry, sgz + irz) - m11*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stressyy_hxiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyy_hxiyhz(sgx + irx, sgy + iry, sgz + irz) - m22*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stresszz_hxiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stresszz_hxiyhz(sgx + irx, sgy + iry, sgz + irz) - m33*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stressxy_hxiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxy_hxiyhz(sgx + irx, sgy + iry, sgz + irz) - m12*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stressxz_hxiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxz_hxiyhz(sgx + irx, sgy + iry, sgz + irz) - m13*amp &
                                            *sgmtr%srcr(k)%interp_hx(irx) &
                                            *sgmtr%srcr(k)%interp_iy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stressyz_hxiyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyz_hxiyhz(sgx + irx, sgy + iry, sgz + irz) - m23*amp &
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
                                    ! ix-hy-hz
                                    if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                        stressxx_ixhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxx_ixhyhz(sgx + irx, sgy + iry, sgz + irz) - m11*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stressyy_ixhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyy_ixhyhz(sgx + irx, sgy + iry, sgz + irz) - m22*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stresszz_ixhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stresszz_ixhyhz(sgx + irx, sgy + iry, sgz + irz) - m33*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stressxy_ixhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxy_ixhyhz(sgx + irx, sgy + iry, sgz + irz) - m12*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stressxz_ixhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressxz_ixhyhz(sgx + irx, sgy + iry, sgz + irz) - m13*amp &
                                            *sgmtr%srcr(k)%interp_ix(irx) &
                                            *sgmtr%srcr(k)%interp_hy(iry) &
                                            *sgmtr%srcr(k)%interp_hz(irz)
                                        stressyz_ixhyhz(sgx + irx, sgy + iry, sgz + irz) = &
                                            stressyz_ixhyhz(sgx + irx, sgy + iry, sgz + irz) - m23*amp &
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
