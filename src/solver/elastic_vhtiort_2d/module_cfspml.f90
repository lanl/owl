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


module elastic_vhtiort_2d_cfspml

    use libflit
    use elastic_vhtiort_2d_vars

    implicit none

    ! PML coefficients
    ! Strangely, the values of these parameters are not
    ! consistent with those reported in existing references
    ! e.g., nk=3, nd=3, na=1, R0=5.0e-9
    integer, parameter :: npower_k = 2
    integer, parameter :: npower_d = 2
    integer, parameter :: npower_a = 1
    real :: eigenderiv_max = -0.01

    ! Scaling factor1
    real :: alphamax
    integer :: npower
    real :: R0
    real :: kmax

    real, allocatable, dimension(:, :) :: axii
    real, allocatable, dimension(:, :) :: azii
    real, allocatable, dimension(:, :) :: bxii
    real, allocatable, dimension(:, :) :: bzii
    real, allocatable, dimension(:, :) :: kxii
    real, allocatable, dimension(:, :) :: kzii

    real, allocatable, dimension(:, :) :: axhh
    real, allocatable, dimension(:, :) :: azhh
    real, allocatable, dimension(:, :) :: bxhh
    real, allocatable, dimension(:, :) :: bzhh
    real, allocatable, dimension(:, :) :: kxhh
    real, allocatable, dimension(:, :) :: kzhh

    real, allocatable, dimension(:, :) :: axhi
    real, allocatable, dimension(:, :) :: azhi
    real, allocatable, dimension(:, :) :: bxhi
    real, allocatable, dimension(:, :) :: bzhi
    real, allocatable, dimension(:, :) :: kxhi
    real, allocatable, dimension(:, :) :: kzhi

    real, allocatable, dimension(:, :) :: axih
    real, allocatable, dimension(:, :) :: azih
    real, allocatable, dimension(:, :) :: bxih
    real, allocatable, dimension(:, :) :: bzih
    real, allocatable, dimension(:, :) :: kxih
    real, allocatable, dimension(:, :) :: kzih

    real, allocatable, dimension(:) :: dampratio_left
    real, allocatable, dimension(:) :: dampratio_right
    real, allocatable, dimension(:) :: dampratio_top
    real, allocatable, dimension(:) :: dampratio_bottom

contains

    !
    !> Calculate the eigenvalue for M-PML
    !
    function mpml_eigenvalue(c11, c13, c15, c33, c35, c55, rho, whichmode, k1, k3)

        ! arguments
        real, intent(in) :: c11, c13, c15, c33, c35, c55, rho, k1, k3
        integer, intent(in) :: whichmode
        complex :: mpml_eigenvalue

        ! Get the eigvalue
        select case (whichmode)

            case (1)
                mpml_eigenvalue = -(sqrt(cmplx(-((((c11 + c55)*k1**2 + 2*(c15 + c35)*k1*k3 + (c33 + c55)*k3**2)*rho + &
                    sqrt((((c11 + c55)*k1**2 + 2*(c15 + c35)*k1*k3 + (c33 + c55)*k3**2)**2 + &
                    4*((c15**2 - c11*c55)*k1**4 + 2*(c13*c15 - c11*c35)*k1**3*k3 + &
                    (c13**2 - c11*c33 - 2*c15*c35 + 2*c13*c55)*k1**2*k3**2 + &
                    2*(-(c15*c33) + c13*c35)*k1*k3**3 + (c35**2 - c33*c55)*k3**4))*rho**2) &
                    )/rho**2)))/sqrt(2.0))

            case (2)
                mpml_eigenvalue = -(sqrt(cmplx((-(((c11 + c55)*k1**2 + 2*(c15 + c35)*k1*k3 + (c33 + c55)*k3**2)*rho) + &
                    sqrt((((c11 + c55)*k1**2 + 2*(c15 + c35)*k1*k3 + (c33 + c55)*k3**2)**2 + &
                    4*((c15**2 - c11*c55)*k1**4 + 2*(c13*c15 - c11*c35)*k1**3*k3 + &
                    (c13**2 - c11*c33 - 2*c15*c35 + 2*c13*c55)*k1**2*k3**2 + &
                    2*(-(c15*c33) + c13*c35)*k1*k3**3 + (c35**2 - c33*c55)*k3**4))*rho**2))/ &
                    rho**2))/sqrt(2.0))

        end select

    end function mpml_eigenvalue

    !
    !> Calculate the eigenvalue derivative for M-PML
    !
    function mpml_eigenvalue_derivative(c11, c13, c15, c33, c35, c55, rho, direction, eigval, ratio, k1, k3)

        ! Parameters
        real, intent(in) :: c11, c13, c15, c33, c35, c55, rho
        character(len=*), intent(in) :: direction
        real, intent(in) :: ratio, k1, k3
        complex, intent(in) :: eigval
        complex :: mpml_eigenvalue_derivative

        ! Calcualte the eigevalue derivative
        select case (direction)
            case ('x')
                mpml_eigenvalue_derivative = &
                    (2*c15*c35*k1**2*k3**2 + 3*c15*c33*k1*k3**3 - 2*c35**2*k3**4 &
                    + 2*c33*c55*k3**4 - 2*c15**2*k1**4*ratio + &
                    2*c15*c35*k1**2*k3**2*ratio + c15*c33*k1*k3**3*ratio &
                    - c13**2*k1**2*k3**2*(1 + ratio) - &
                    c13*k1*k3*(c15*k1**2*(1 + 3*ratio) &
                    + k3*(2*c55*k1*(1 + ratio) + c35*k3*(3 + ratio))) + &
                    eigval**2*(k3*(3*c15*k1*(1 + ratio) + 3*c35*k1*(1 + ratio) &
                    + c33*k3*(2 + ratio)) + c55*(k3**2*(2 + ratio) + k1**2*(1 + 2*ratio)))* &
                    rho + 2*eigval**4*(1 + ratio)*rho**2 + c11*k1**2* &
                    (2*c55*k1**2*ratio + c33*k3**2*(1 + ratio) &
                    + c35*k1*k3*(1 + 3*ratio) + eigval**2*(1 + 2*ratio)*rho))/ &
                    (2*((c15**2 - c11*c55)*k1**4 + 2*(c13*c15 - c11*c35)*k1**3*k3 &
                    + (c13**2 - c11*c33 - 2*c15*c35 + 2*c13*c55)*k1**2*k3**2 + &
                    2*(-(c15*c33) + c13*c35)*k1*k3**3 + (c35**2 - c33*c55)*k3**4) - &
                    3*eigval**2*((c11 + c55)*k1**2 + 2*(c15 + c35)*k1*k3 &
                    + (c33 + c55)*k3**2)*rho - 4*eigval**4*rho**2)

            case ('z')
                mpml_eigenvalue_derivative = &
                    (-2*c15**2*k1**4 + k3**2*(-2*(c35**2 &
                    - c33*c55)*k3**2*ratio - c13**2*k1**2*(1 + ratio) - &
                    c13*k1*(2*c55*k1*(1 + ratio) + c35*k3*(1 + 3*ratio))) + &
                    eigval**2*(k3*(3*c35*k1*(1 + ratio) + c33*k3*(1 + 2*ratio)) &
                    + c55*(k1**2*(2 + ratio) + k3**2*(1 + 2*ratio)))*rho + &
                    2*eigval**4*(1 + ratio)*rho**2 + c15*k1*k3*(-(c13*k1**2*(3 + ratio)) &
                    + k3*(2*c35*k1*(1 + ratio) + c33*k3*(1 + 3*ratio)) + &
                    3*eigval**2*(1 + ratio)*rho) + c11*k1**2*(2*c55*k1**2 &
                    + k3*(c33*k3*(1 + ratio) + c35*k1*(3 + ratio)) + &
                    eigval**2*(2 + ratio)*rho))/ &
                    (2*((c15**2 - c11*c55)*k1**4 + 2*(c13*c15 - c11*c35)*k1**3*k3 &
                    + (c13**2 - c11*c33 - 2*c15*c35 + 2*c13*c55)*k1**2*k3**2 + &
                    2*(-(c15*c33) + c13*c35)*k1*k3**3 + (c35**2 - c33*c55)*k3**4) - &
                    3*eigval**2*((c11 + c55)*k1**2 + 2*(c15 + c35)*k1*k3 &
                    + (c33 + c55)*k3**2)*rho - 4*eigval**4*rho**2)

        end select

    end function mpml_eigenvalue_derivative

    !
    !> Get the x/z or z/x damping ratio for M-PML
    !
    !  Refs: Method adopted from Meza-Fajardo and Papageorgiou (2008)
    !        with some extension to deal with genenral anisotropic media
    !
    function damp_ratio(c11, c13, c15, c33, c35, c55, rho, direction) result(ratio)

        ! arguments
        real :: c11, c13, c15, c33, c35, c55, rho
        character(len=*) :: direction

        real :: kx, kz, dtheta, derig, ratio, theta
        complex :: eigval1, eigval2, derig1, derig2
        integer :: i
        real :: alpha, dratio

        ! Devide angle
        dtheta = const_pi/180.0

        ! Critical value
        alpha = eigenderiv_max
        dratio = 0.002

        ! Calculate the ratio
        ratio = 0.0
        do i = 1, 180

            ! theta ranging from 0 to pi: theta = 0 is z-axis
            theta = (i + 0.5)*dtheta
            kx = sin(theta)
            kz = cos(theta)

            ! Eigenvalue
            eigval1 = mpml_eigenvalue(c11, c13, c15, c33, c35, c55, rho, 1, kx, kz)
            eigval2 = mpml_eigenvalue(c11, c13, c15, c33, c35, c55, rho, 2, kx, kz)

            ! Derivative of eigenvalue
            derig1 = mpml_eigenvalue_derivative(c11*1.0e-9, c13*1.0e-9, c15*1.0e-9, &
                c33*1.0e-9, c35*1.0e-9, c55*1.0e-9, rho*1.0e-9, direction, eigval1, ratio, kx, kz)
            derig2 = mpml_eigenvalue_derivative(c11*1.0e-9, c13*1.0e-9, c15*1.0e-9, &
                c33*1.0e-9, c35*1.0e-9, c55*1.0e-9, rho*1.0e-9, direction, eigval2, ratio, kx, kz)
            derig = max(real(derig1), real(derig2))

            do while (derig >= alpha .and. ratio < 1.0)

                ! Increase ratio until the derivatives smaller than alpha
                ratio = ratio + dratio

                derig1 = mpml_eigenvalue_derivative(c11*1.0e-9, c13*1.0e-9, c15*1.0e-9, &
                    c33*1.0e-9, c35*1.0e-9, c55*1.0e-9, rho*1.0e-9, direction, eigval1, ratio, kx, kz)
                derig2 = mpml_eigenvalue_derivative(c11*1.0e-9, c13*1.0e-9, c15*1.0e-9, &
                    c33*1.0e-9, c35*1.0e-9, c55*1.0e-9, rho*1.0e-9, direction, eigval2, ratio, kx, kz)
                derig = max(real(derig1), real(derig2))

            end do

        end do

    end function damp_ratio

    !
    !> Compute the MPML damping ratios
    !
    subroutine compute_damping_profile_ratio

        integer :: i, j, bindex

        call alloc_array(dampratio_left, [-pml + 1 - 1, nz + pml + 1])
        call alloc_array(dampratio_right, [-pml + 1 - 1, nz + pml + 1])
        call alloc_array(dampratio_top, [-pml + 1 - 1, nx + pml + 1])
        call alloc_array(dampratio_bottom, [-pml + 1 - 1, nx + pml + 1])

        ! Left boundary
        bindex = 1
        !$omp parallel do private(j)
        do j = -pml + 1 - 1, nz + pml + 1
            dampratio_left(j) = damp_ratio(c11(bindex, j), c13(bindex, j), 0.0, c33(bindex, j), &
                0.0, c55(bindex, j), rho(bindex, j), 'x')
        end do
        !$omp end parallel do

        ! Right boundary
        bindex = nx
        !$omp parallel do private(j)
        do j = -pml + 1 - 1, nz + pml + 1
            dampratio_right(j) = damp_ratio(c11(bindex, j), c13(bindex, j), 0.0, c33(bindex, j), &
                0.0, c55(bindex, j), rho(bindex, j), 'x')
        end do
        !$omp end parallel do

        ! Top boundary
        bindex = 1
        !$omp parallel do private(i)
        do i = -pml + 1 - 1, nx + pml + 1
            dampratio_top(i) = damp_ratio(c11(i, bindex), c13(i, bindex), 0.0, c33(i, bindex), &
                0.0, c55(i, bindex), rho(i, bindex), 'z')
        end do
        !$omp end parallel do

        ! Bottom boundary
        bindex = nz
        !$omp parallel do private(i)
        do i = -pml + 1 - 1, nx + pml + 1
            dampratio_bottom(i) = damp_ratio(c11(i, bindex), c13(i, bindex), 0.0, c33(i, bindex), &
                0.0, c55(i, bindex), rho(i, bindex), 'z')
        end do
        !$omp end parallel do

        call warn(date_time_compact()//' Maximum MPML damping raitos = ' &
            //num2str(maxval(dampratio_left), '(f5.3)')//' ' &
            //num2str(maxval(dampratio_right), '(f5.3)')//' ' &
            //num2str(maxval(dampratio_top), '(f5.3)')//' ' &
            //num2str(maxval(dampratio_bottom), '(f5.3)'))

    end subroutine

    !
    !> Compute MPML damping coefficients
    !
    subroutine damp_coef(vpx, vpz, ax, az, bx, bz, kx, kz, xdist, zdist, dx, dz, ratiox, ratioz)

        real, intent(in) :: vpx, vpz, dx, dz, ratiox, ratioz
        real, intent(inout) :: ax, az, bx, bz, kx, kz
        real, intent(in) :: xdist, zdist

        real :: dampx, dampz
        real :: alphax, alphaz
        real :: nd
        real :: wx, wz

        ax = 0.0
        az = 0.0
        bx = 0.0
        bz = 0.0
        kx = 1.0
        kz = 1.0

        nd = xdist
        dampx = vpx*(npower_d + 1.0)*log(1.0/R0)/(2.0*pml*dx)*nd**npower_d
        nd = zdist
        dampz = vpz*(npower_d + 1.0)*log(1.0/R0)/(2.0*pml*dz)*nd**npower_d

        wx = dampx + ratioz*dampz
        wz = dampz + ratiox*dampx
        dampx = wx
        dampz = wz

        nd = xdist
        ! Diminishing alpha (Roden and Gedney, 2000); residual alpha to ensure long-time stability
        alphax = max(alphamax*(1.0 - nd**npower_a), 0.25*alphamax)
        kx = 1.0 + (kmax - 1.0)*nd**npower_k
        ax = -2.0*abs(dt)*dampx/kx/(2.0 + abs(dt)*(alphax + dampx/kx))
        bx = (2.0 - abs(dt)*(alphax + dampx/kx))/(2.0 + abs(dt)*(alphax + dampx/kx))

        nd = zdist
        ! Diminishing alpha (Roden and Gedney, 2000); residual alpha to ensure long-time stability
        alphaz = max(alphamax*(1.0 - nd**npower_a), 0.25*alphamax)
        kz = 1.0 + (kmax - 1.0)*nd**npower_k
        az = -2.0*abs(dt)*dampz/kz/(2.0 + abs(dt)*(alphaz + dampz/kz))
        bz = (2.0 - abs(dt)*(alphaz + dampz/kz))/(2.0 + abs(dt)*(alphaz + dampz/kz))

    end subroutine

    !
    !> Compute MPML damping coefficients
    !
    subroutine compute_cfspml_damping_coef

        integer :: i, j
        real :: ax, az, bx, bz, kx, kz
        real :: xdisti, xdisth, zdisti, zdisth
        real :: ratiox, ratioz

        ! Compute CFS-MPML coefficients
        alphamax = 1.0*maxval(sgmtr%srcr(:)%f0)*const_pi
        R0 = 1.0e-5
        kmax = 1.0

        if (aniso_param == 'iso') then
            call alloc_array(dampratio_left, [-pml + 1 - 1, nz + pml + 1])
            call alloc_array(dampratio_right, [-pml + 1 - 1, nz + pml + 1])
            call alloc_array(dampratio_top, [-pml + 1 - 1, nx + pml + 1])
            call alloc_array(dampratio_bottom, [-pml + 1 - 1, nx + pml + 1])
            dampratio_left = 0
            dampratio_right = 0
            dampratio_top = 0
            dampratio_bottom = 0
        else
            call compute_damping_profile_ratio
        end if

        ! Allocate memory for coefficient arrays
        call alloc_array(axii, [1, nx, 1, nz], pad=pml)
        call alloc_array(axih, [1, nx, 1, nz], pad=pml)
        call alloc_array(axhi, [1, nx, 1, nz], pad=pml)
        call alloc_array(axhh, [1, nx, 1, nz], pad=pml)

        call alloc_array(bxii, [1, nx, 1, nz], pad=pml)
        call alloc_array(bxih, [1, nx, 1, nz], pad=pml)
        call alloc_array(bxhi, [1, nx, 1, nz], pad=pml)
        call alloc_array(bxhh, [1, nx, 1, nz], pad=pml)

        call alloc_array(kxii, [1, nx, 1, nz], pad=pml)
        call alloc_array(kxih, [1, nx, 1, nz], pad=pml)
        call alloc_array(kxhi, [1, nx, 1, nz], pad=pml)
        call alloc_array(kxhh, [1, nx, 1, nz], pad=pml)

        call alloc_array(azii, [1, nx, 1, nz], pad=pml)
        call alloc_array(azih, [1, nx, 1, nz], pad=pml)
        call alloc_array(azhi, [1, nx, 1, nz], pad=pml)
        call alloc_array(azhh, [1, nx, 1, nz], pad=pml)

        call alloc_array(bzii, [1, nx, 1, nz], pad=pml)
        call alloc_array(bzih, [1, nx, 1, nz], pad=pml)
        call alloc_array(bzhi, [1, nx, 1, nz], pad=pml)
        call alloc_array(bzhh, [1, nx, 1, nz], pad=pml)

        call alloc_array(kzii, [1, nx, 1, nz], pad=pml)
        call alloc_array(kzih, [1, nx, 1, nz], pad=pml)
        call alloc_array(kzhi, [1, nx, 1, nz], pad=pml)
        call alloc_array(kzhh, [1, nx, 1, nz], pad=pml)

        kxii = 1.0
        kxih = 1.0
        kxhi = 1.0
        kxhh = 1.0

        kzii = 1.0
        kzih = 1.0
        kzhi = 1.0
        kzhh = 1.0

        if (yn_free_surface) then
            ! Free-surface modeling uses depth-varying mesh

            !$omp parallel do private(i, j, ax, bx, kx, az, bz, kz, &
                !$omp xdisti, xdisth, zdisti, zdisth, ratiox, ratioz) collapse(2) schedule(auto)
            do j = -pml + 1, nz + pml
                do i = -pml + 1, nx + pml

                    xdisti = 0.0d0
                    xdisth = 0.0d0
                    zdisti = 0.0d0
                    zdisth = 0.0d0

                    if (i <= 1) then
                        xdisti = abs((i - 1 + 0.0d0)/pml)
                        xdisth = abs((i - 1 - 0.5d0)/pml)
                    else if (i >= nx) then
                        xdisti = abs((i - nx + 0.0d0)/pml)
                        xdisth = abs((i - nx - 0.5d0)/pml)
                    end if

                    if (j <= 1) then
                        zdisti = abs((j - 1 + 0.0d0)/pml)
                        zdisth = abs((j - 1 - 0.5d0)/pml)
                    else if (j >= nz) then
                        zdisti = abs((j - nz + 0.0d0)/pml)
                        zdisth = abs((j - nz - 0.5d0)/pml)
                    end if

                    ! Compute only for boundaries
                    if (.not. (i > 1 .and. i < nx .and. j > 1 .and. j < nz)) then

                        ratiox = 0.0d0
                        ratioz = 0.0d0

                        ! integer-integer
                        if (i <= 1) then
                            ratiox = dampratio_left(j)
                        else if (i >= nx) then
                            ratiox = dampratio_right(j)
                        end if
                        if (j <= 1) then
                            ratioz = dampratio_top(i)
                        else if (j >= nz) then
                            ratioz = dampratio_bottom(i)
                        end if

                        ! integer-integer
                        call damp_coef(pmlvp, pmlvp, ax, az, bx, bz, kx, kz, xdisti, zdisti, dx, dz_i(j), ratiox, ratioz)
                        axii(i, j) = ax
                        bxii(i, j) = bx
                        kxii(i, j) = kx
                        azii(i, j) = az
                        bzii(i, j) = bz
                        kzii(i, j) = kz

                        ! integer-half
                        call damp_coef(pmlvp, pmlvp, ax, az, bx, bz, kx, kz, xdisti, zdisth, dx, dz_i(j), ratiox, ratioz)
                        axih(i, j) = ax
                        bxih(i, j) = bx
                        kxih(i, j) = kx
                        azih(i, j) = az
                        bzih(i, j) = bz
                        kzih(i, j) = kz

                        ! half-integer
                        call damp_coef(pmlvp, pmlvp, ax, az, bx, bz, kx, kz, xdisth, zdisti, dx, dz_i(j), ratiox, ratioz)
                        axhi(i, j) = ax
                        bxhi(i, j) = bx
                        kxhi(i, j) = kx
                        azhi(i, j) = az
                        bzhi(i, j) = bz
                        kzhi(i, j) = kz

                        ! half-half
                        call damp_coef(pmlvp, pmlvp, ax, az, bx, bz, kx, kz, xdisth, zdisth, dx, dz_i(j), ratiox, ratioz)
                        axhh(i, j) = ax
                        bxhh(i, j) = bx
                        kxhh(i, j) = kx
                        azhh(i, j) = az
                        bzhh(i, j) = bz
                        kzhh(i, j) = kz

                    end if

                end do
            end do
            !$omp end parallel do

        else

            !$omp parallel do private(i, j, ax, bx, kx, az, bz, kz, &
                !$omp xdisti, xdisth, zdisti, zdisth, ratiox, ratioz) collapse(2) schedule(auto)
            do j = -pml + 1, nz + pml
                do i = -pml + 1, nx + pml

                    xdisti = 0.0d0
                    xdisth = 0.0d0
                    zdisti = 0.0d0
                    zdisth = 0.0d0

                    if (i <= 1) then
                        xdisti = abs((i - 1 + 0.0d0)/pml)
                        xdisth = abs((i - 1 - 0.5d0)/pml)
                    else if (i >= nx) then
                        xdisti = abs((i - nx + 0.0d0)/pml)
                        xdisth = abs((i - nx - 0.5d0)/pml)
                    end if

                    if (j <= 1) then
                        zdisti = abs((j - 1 + 0.0d0)/pml)
                        zdisth = abs((j - 1 - 0.5d0)/pml)
                    else if (j >= nz) then
                        zdisti = abs((j - nz + 0.0d0)/pml)
                        zdisth = abs((j - nz - 0.5d0)/pml)
                    end if

                    ! Compute only for boundaries
                    if (.not. (i > 1 .and. i < nx .and. j > 1 .and. j < nz)) then

                        ratiox = 0.0d0
                        ratioz = 0.0d0

                        ! integer-integer
                        if (i <= 1) then
                            ratiox = dampratio_left(j)
                        else if (i >= nx) then
                            ratiox = dampratio_right(j)
                        end if
                        if (j <= 1) then
                            ratioz = dampratio_top(i)
                        else if (j >= nz) then
                            ratioz = dampratio_bottom(i)
                        end if

                        ! integer-integer
                        call damp_coef(pmlvp, pmlvp, ax, az, bx, bz, kx, kz, xdisti, zdisti, dx, dz, ratiox, ratioz)
                        axii(i, j) = ax
                        bxii(i, j) = bx
                        kxii(i, j) = kx
                        azii(i, j) = az
                        bzii(i, j) = bz
                        kzii(i, j) = kz

                        ! integer-half
                        call damp_coef(pmlvp, pmlvp, ax, az, bx, bz, kx, kz, xdisti, zdisth, dx, dz, ratiox, ratioz)
                        axih(i, j) = ax
                        bxih(i, j) = bx
                        kxih(i, j) = kx
                        azih(i, j) = az
                        bzih(i, j) = bz
                        kzih(i, j) = kz

                        ! half-integer
                        call damp_coef(pmlvp, pmlvp, ax, az, bx, bz, kx, kz, xdisth, zdisti, dx, dz, ratiox, ratioz)
                        axhi(i, j) = ax
                        bxhi(i, j) = bx
                        kxhi(i, j) = kx
                        azhi(i, j) = az
                        bzhi(i, j) = bz
                        kzhi(i, j) = kz

                        ! half-half
                        call damp_coef(pmlvp, pmlvp, ax, az, bx, bz, kx, kz, xdisth, zdisth, dx, dz, ratiox, ratioz)
                        axhh(i, j) = ax
                        bxhh(i, j) = bx
                        kxhh(i, j) = kx
                        azhh(i, j) = az
                        bzhh(i, j) = bz
                        kzhh(i, j) = kz

                    end if

                end do
            end do
            !$omp end parallel do

        end if

    end subroutine

end module
