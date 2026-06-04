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

module mod_anisotropy

    use libflit

    implicit none

    interface rotate_cij
        module procedure :: cij_rotate_2d
        module procedure :: cij_rotate_3d
    end interface

    interface thomsen_to_cij
        module procedure :: thomsen_to_cij2d
        module procedure :: thomsen_to_cij3d
    end interface

    interface alkhalifah_tsvankin_to_cij
        module procedure :: alkhalifah_tsvankin_to_cij2d
        module procedure :: alkhalifah_tsvankin_to_cij3d
    end interface

    interface cij_to_sij
        module procedure :: cij_to_sij_2d
        module procedure :: cij_to_sij_3d
    end interface

    private
    public :: rotate_cij
    public :: thomsen_to_cij
    public :: alkhalifah_tsvankin_to_cij
    public :: cij_to_sij
    public :: min_max_phase_velocity_2d
    public :: min_max_phase_velocity_3d

contains

    !
    !> Convert elasticity to compliance for 2D elastic media
    !
    subroutine cij_to_sij_2d(c11, c13, c15, c33, c35, c55, &
            s11, s13, s15, s33, s35, s55)

        real, intent(in) :: c11, c13, c15, c33, c35, c55
        real, intent(inout) ::s11, s13, s15, s33, s35, s55
        real :: temp1, temp2

        ! Analytical expressions for the elasticity matrix inverse
        temp1 = c15**2*c33 - 2*c13*c15*c35 + c11*c35**2 + c13**2*c55 - c11*c33*c55
        temp2 = c15**2*c33 - 2*c13*c15*c35 + c13**2*c55 + c11*(c35**2 - c33*c55)

        s11 = (c35**2 - c33*c55)/temp1
        s13 = (-(c15*c35) + c13*c55)/temp2
        s15 = (c15*c33 - c13*c35)/temp1
        s33 = (c15**2 - c11*c55)/temp1
        s35 = (-(c13*c15) + c11*c35)/temp2
        s55 = (c13**2 - c11*c33)/temp1

    end subroutine

    !
    !> Convert elasticity to compliance for 3D elastic media
    !
    subroutine cij_to_sij_3d(c11, c12, c13, c14, c15, c16, &
            c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, &
            c55, c56, c66, s11, s12, s13, s14, s15, s16, &
            s22, s23, s24, s25, s26, s33, s34, s35, s36, s44, s45, s46, &
            s55, s56, s66)

        real, intent(in) :: c11, c12, c13, c14, c15, c16, &
            c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, &
            c55, c56, c66
        real, intent(out) :: s11, s12, s13, s14, s15, s16, &
            s22, s23, s24, s25, s26, s33, s34, s35, s36, s44, s45, s46, &
            s55, s56, s66

        real, dimension(1:6, 1:6) :: s

        s = inv(reshape([ &
            c11, c12, c13, c14, c15, c16, &
            c12, c22, c23, c24, c25, c26, &
            c13, c23, c33, c34, c35, c36, &
            c14, c24, c34, c44, c45, c46, &
            c15, c25, c35, c45, c55, c56, &
            c16, c26, c36, c46, c56, c66], shape=[6, 6]))
        s11 = s(1, 1)
        s12 = s(1, 2)
        s13 = s(1, 3)
        s14 = s(1, 4)
        s15 = s(1, 5)
        s16 = s(1, 6)
        s22 = s(2, 2)
        s23 = s(2, 3)
        s24 = s(2, 4)
        s25 = s(2, 5)
        s26 = s(2, 6)
        s33 = s(3, 3)
        s34 = s(3, 4)
        s35 = s(3, 5)
        s36 = s(3, 6)
        s44 = s(4, 4)
        s45 = s(4, 5)
        s46 = s(4, 6)
        s55 = s(5, 5)
        s56 = s(5, 6)
        s66 = s(6, 6)

    end subroutine

    pure function bond_matrix(r) result(b)

        real, dimension(1:3, 1:3), intent(in) :: r
        real, allocatable, dimension(:, :) :: b

        allocate (b(1:6, 1:6))
        b(1, :) = [ &
            r(1, 1)**2, &
            r(1, 2)**2, &
            r(1, 3)**2, &
            2.0*r(1, 2)*r(1, 3), &
            2.0*r(1, 1)*r(1, 3), &
            2.0*r(1, 1)*r(1, 2)]
        b(2, :) = [ &
            r(2, 1)**2, &
            r(2, 2)**2, &
            r(2, 3)**2, &
            2.0*r(2, 2)*r(2, 3), &
            2.0*r(2, 1)*r(2, 3), &
            2.0*r(2, 1)*r(2, 2)]
        b(3, :) = [ &
            r(3, 1)**2, &
            r(3, 2)**2, &
            r(3, 3)**2, &
            2.0*r(3, 2)*r(3, 3), &
            2.0*r(3, 1)*r(3, 3), &
            2.0*r(3, 1)*r(3, 2)]
        b(4, :) = [ &
            r(2, 1)*r(3, 1), &
            r(2, 2)*r(3, 2), &
            r(2, 3)*r(3, 3), &
            r(2, 2)*r(3, 3) + r(2, 3)*r(3, 2), &
            r(2, 1)*r(3, 3) + r(2, 3)*r(3, 1), &
            r(2, 1)*r(3, 2) + r(2, 2)*r(3, 1)]
        b(5, :) = [ &
            r(1, 1)*r(3, 1), &
            r(1, 2)*r(3, 2), &
            r(1, 3)*r(3, 3), &
            r(1, 2)*r(3, 3) + r(1, 3)*r(3, 2), &
            r(1, 1)*r(3, 3) + r(1, 3)*r(3, 1), &
            r(1, 1)*r(3, 2) + r(1, 2)*r(3, 1)]
        b(6, :) = [ &
            r(1, 1)*r(2, 1), &
            r(1, 2)*r(2, 2), &
            r(1, 3)*r(2, 3), &
            r(1, 2)*r(2, 3) + r(1, 3)*r(2, 2), &
            r(1, 1)*r(2, 3) + r(1, 3)*r(2, 1), &
            r(1, 1)*r(2, 2) + r(1, 2)*r(2, 1)]

    end function

    !
    !> Rotate Cij matrix counter-clockwise around y-axis by an angle thetay
    !
    elemental subroutine cij_rotate_2d(c11, c13, c15, c33, c35, c55, thetay)

        real, intent(inout) :: c11, c13, c15, c33, c35, c55
        real, intent(in) :: thetay

        real, allocatable, dimension(:, :) :: c, a
        real :: c12, c14, c16
        real :: c22, c23, c24, c25, c26
        real :: c34, c36
        real :: c44, c45, c46
        real :: c56
        real :: c66

        c12 = 0.0
        c14 = 0.0
        c16 = 0.0
        c22 = 0.0
        c23 = 0.0
        c24 = 0.0
        c25 = 0.0
        c26 = 0.0
        c34 = 0.0
        c36 = 0.0
        c44 = 0.0
        c45 = 0.0
        c46 = 0.0
        c56 = 0.0
        c66 = 0.0

        ! Rotation matrix; For 2D case, the rotation order does not matter
        a = bond_matrix(rotation_matrix(thetay, 'y'))

        ! Elasticity matrix
        c = reshape([ &
            c11, c12, c13, c14, c15, c16, &
            c12, c22, c23, c24, c25, c26, &
            c13, c23, c33, c34, c35, c36, &
            c14, c24, c34, c44, c45, c46, &
            c15, c25, c35, c45, c55, c56, &
            c16, c26, c36, c46, c56, c66], [6, 6])

        ! Rotation
        c = matmul(a, matmul(c, transpose(a)))

        ! Take relevant elements for 2D scenario
        c11 = c(1, 1)
        c13 = c(1, 3)
        c15 = c(1, 5)
        c33 = c(3, 3)
        c35 = c(3, 5)
        c55 = c(5, 5)

    end subroutine

    !
    !> Rotate Cij matrix
    !
    !> thetax Counter-clockwise rotation angle around x-axis
    !> thetay Counter-clockwise rotation angle around y-axis
    !> thetaz Counter-clockwise rotation angle around z-axis
    !> order Rotation order, =xyz (default), xzy, yxz, yzx, zxy, zyx
    !>
    !> For VTI medium, only two angles are needed (theta_y, theta_z)
    !> For orthorhombic medium, three may be needed
    !
    elemental subroutine cij_rotate_3d(c11, c12, c13, c14, c15, c16, &
            c22, c23, c24, c25, c26, c33, c34, c35, c36, &
            c44, c45, c46, c55, c56, c66, thetax, thetay, thetaz, &
            order)

        real, intent(inout) :: c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66
        real, intent(in) :: thetax, thetay, thetaz
        character(len=*), intent(in), optional :: order

        real, allocatable, dimension(:, :) :: c, a
        character(len=3) :: rotation_order

        if (present(order)) then
            rotation_order = order
        else
            rotation_order = 'xyz'
        end if

        ! Rotation matrix
        ! Note that to rotate a VTI medium, the first rotation cannot be around z;
        ! otherwise, the azimuth angle does not act at all since VTI is symmetric in the x-y plane
        a = bond_matrix(rotation_matrix([thetax, thetay, thetaz], rotation_order))

        ! Elasticity matrix
        c = reshape([ &
            c11, c12, c13, c14, c15, c16, &
            c12, c22, c23, c24, c25, c26, &
            c13, c23, c33, c34, c35, c36, &
            c14, c24, c34, c44, c45, c46, &
            c15, c25, c35, c45, c55, c56, &
            c16, c26, c36, c46, c56, c66], [6, 6])

        ! Rotation
        c = matmul(a, matmul(c, transpose(a)))

        ! Take relevant elements for 2D
        c11 = c(1, 1)
        c12 = c(1, 2)
        c13 = c(1, 3)
        c14 = c(1, 4)
        c15 = c(1, 5)
        c16 = c(1, 6)
        c22 = c(2, 2)
        c23 = c(2, 3)
        c24 = c(2, 4)
        c25 = c(2, 5)
        c26 = c(2, 6)
        c33 = c(3, 3)
        c34 = c(3, 4)
        c35 = c(3, 5)
        c36 = c(3, 6)
        c44 = c(4, 4)
        c45 = c(4, 5)
        c46 = c(4, 6)
        c55 = c(5, 5)
        c56 = c(5, 6)
        c66 = c(6, 6)

    end subroutine

    !
    !> Convert Thomsen parameters to 2D TTI Cij
    !
    elemental subroutine thomsen_to_cij2d( &
            vp, vs, rho, eps, del, the, c11, c13, c15, c33, c35, c55)

        real, intent(in) :: vp, vs, rho, eps, del, the
        real, intent(out) :: c11, c13, c15, c33, c35, c55

        real :: tmp

        c33 = vp**2*rho
        c55 = vs**2*rho
        c11 = c33*(1 + 2*eps)
        tmp = 2*c33*(c33 - c55)*del + (c33 - c55)**2
        c13 = sqrt(abs(tmp)) - c55
        c15 = 0.0
        c35 = 0.0

        if (the /= 0) then
            ! Counter-clockwise rotation
            call cij_rotate_2d(c11, c13, c15, c33, c35, c55, -the)
        end if

    end subroutine

    !
    !> Convert Thomsen parameters to 3D TTI Cij
    !
    elemental subroutine thomsen_to_cij3d( &
            vp, vs, rho, eps, del, gam, the, phi, &
            c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66)

        real, intent(in) :: vp, vs, rho, eps, del, gam, the, phi
        real, intent(out) :: c11, c12, c13, c14, c15, c16, &
            c22, c23, c24, c25, c26, &
            c33, c34, c35, c36, c44, c45, c46, &
            c55, c56, c66

        real :: tmp

        c33 = vp**2*rho
        c44 = vs**2*rho
        c11 = c33*(1 + 2*eps)
        c66 = c44*(1 + 2*gam)
        c12 = c11 - 2*c66

        tmp = 2*c33*(c33 - c44)*del + (c33 - c44)**2
        c13 = sqrt(abs(tmp)) - c44

        c22 = c11
        c23 = c13
        c55 = c44

        c14 = 0.0
        c15 = 0.0
        c16 = 0.0
        c24 = 0.0
        c25 = 0.0
        c26 = 0.0
        c34 = 0.0
        c35 = 0.0
        c36 = 0.0
        c45 = 0.0
        c46 = 0.0
        c56 = 0.0

        if ((the /= 0 .or. phi /= 0) .and. .not. (eps == 0 .and. del == 0 .and. gam == 0)) then
            ! Counter-clockwise rotation
            call cij_rotate_3d(c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
                c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, 0.0, -the, -phi, 'xyz')
        end if

    end subroutine

    !
    !> Convert Alkhalifah-Tsvankin parameters to 2D TTI Cij
    !
    elemental subroutine alkhalifah_tsvankin_to_cij2d( &
            vp, vs, rho, eps, eta, the, c11, c13, c15, c33, c35, c55)

        real, intent(in) :: vp, vs, rho, eps, eta, the
        real, intent(out) :: c11, c13, c15, c33, c35, c55

        c11 = vp**2*rho
        c55 = vs**2*rho
        c33 = c11/(1 + 2*eps)
        c13 = rho*(sqrt(abs((vp**2/(1 + 2*eps) - vs**2)*(vp**2/(1 + 2*eta) - vs**2))) - vs**2)
        c15 = 0.0
        c35 = 0.0

        if (the /= 0) then
            ! Counter-clockwise rotation
            call cij_rotate_2d(c11, c13, c15, c33, c35, c55, -the)
        end if

    end subroutine

    !
    !> Convert Alkhalifah-Tsvankin parameters to 3D TTI Cij
    !
    elemental subroutine alkhalifah_tsvankin_to_cij3d( &
            vp, vs, rho, eps, eta, gam, the, phi, &
            c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66)

        real, intent(in) :: vp, vs, rho, eps, eta, gam, the, phi
        real, intent(out) :: c11, c12, c13, c14, c15, c16, &
            c22, c23, c24, c25, c26, &
            c33, c34, c35, c36, c44, c45, c46, &
            c55, c56, c66

        c11 = vp**2*rho
        c55 = vs**2*rho
        c33 = c11/(1 + 2*eps)
        c66 = c44*(1 + 2*gam)

        c12 = c11 - 2*c66
        c13 = rho*(sqrt(abs((vp**2/(1 + 2*eps) - vs**2)*(vp**2/(1 + 2*eta) - vs**2))) - vs**2)

        c22 = c11
        c23 = c13
        c55 = c44

        c14 = 0.0
        c15 = 0.0
        c16 = 0.0
        c24 = 0.0
        c25 = 0.0
        c26 = 0.0
        c34 = 0.0
        c35 = 0.0
        c36 = 0.0
        c45 = 0.0
        c46 = 0.0
        c56 = 0.0

        if ((the /= 0 .or. phi /= 0) .and. .not. (eps == 0 .and. eta == 0 .and. gam == 0)) then
            ! Counter-clockwise rotation
            call cij_rotate_3d(c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
                c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, 0.0, -the, -phi, 'xyz')
        end if

    end subroutine


    !
    !> Calculate phase velocities in a 2D anisotropic medium along (k1, k3)
    !
    subroutine christoffel2(c11, c13, c15, c33, c35, c55, rho, k1, k3, qv)

        ! Parameters
        real, intent(in) :: c11, c13, c15, c33, c35, c55, rho, k1, k3
        real, dimension(:), intent(out) :: qv

        real :: a11, a12, a22

        ! k1 = sin(theta), k3 = cos(theta)
        ! Therefore, theta is the angle between wavenumber vector and z-axis
        a11 = c11*k1**2 + 2*c15*k1*k3 + c55*k3**2
        a12 = c15*k1**2 + (c13 + c55)*k1*k3 + c35*k3**2
        a22 = c55*k1**2 + 2*c35*k1*k3 + c33*k3**2

        ! Solve for the eigenvalues
        qv = zeros(2)
        qv(1) = 0.5*(a11 + sqrt(4.0*a12**2 + (a11 - a22)**2) + a22)
        qv(2) = 0.5*(a11 - sqrt(4.0*a12**2 + (a11 - a22)**2) + a22)

        qv = sqrt(qv/rho)

    end subroutine

    !
    !> Solve Kelvin-Christoffel equation
    !
    subroutine christoffel3(c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho, k1, k2, k3, qv)

        real, intent(in) :: c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho
        real, intent(in) :: k1, k2, k3
        real, dimension(:), intent(inout) :: qv

        real :: a11, a12, a13, a22, a23, a33
        real, allocatable, dimension(:, :) :: a, vr
        real, allocatable, dimension(:) :: wr

        ! Elements in Christoffel matrix
        a11 = c11*k1**2 + 2*c16*k1*k2 + c66*k2**2 + 2*c15*k1*k3 + 2*c56*k2*k3 + c55*k3**2
        a12 = c16*k1**2 + (c12 + c66)*k1*k2 + c26*k2**2 + (c14 + c56)*k1*k3 + (c25 + c46)*k2*k3 + c45*k3**2
        a13 = c15*k1**2 + (c14 + c56)*k1*k2 + c46*k2**2 + (c13 + c55)*k1*k3 + (c36 + c45)*k2*k3 + c35*k3**2
        a22 = c66*k1**2 + 2*c26*k1*k2 + c22*k2**2 + 2*c46*k1*k3 + 2*c24*k2*k3 + c44*k3**2
        a23 = c56*k1**2 + (c25 + c46)*k1*k2 + c24*k2**2 + (c36 + c45)*k1*k3 + (c23 + c44)*k2*k3 + c34*k3**2
        a33 = c55*k1**2 + 2*c45*k1*k2 + c44*k2**2 + 2*c35*k1*k3 + 2*c34*k2*k3 + c33*k3**2

        a = reshape([a11, a12, a13, a12, a22, a23, a13, a23, a33], [3, 3])
        wr = zeros(3)
        vr = zeros(3, 3)

        ! Solve the eigenvalue problem
        call eigen_symm3x3(a, wr, vr)

        ! Return phase velocity
        qv = sqrt(wr/rho)

    end subroutine

    !
    !> Compute qP- and qS-wave phase velocities in the range of polar angle [0, pi]
    !> for a 2D anisotropic medium
    !
    subroutine min_max_phase_velocity_2d(c11, c13, c15, c33, c35, c55, rho, n, qp, qs)

        ! Parameters
        real, intent(in) :: c11, c13, c15, c33, c35, c55, rho
        integer, intent(in) :: n
        real, intent(inout) :: qp, qs

        integer :: i
        real :: theta, dtheta, k1, k3
        real, allocatable, dimension(:, :) :: qv

        qv = zeros(n, 2)

        dtheta = const_pi/n

        ! Iteration through different directions
        do i = 1, n

            theta = (i - 1)*dtheta
            k1 = sin(theta)
            k3 = cos(theta)

            ! Take the minimum and maximum phase velocities
            ! as the qS- and qP-wave phase velocities, respectively
            call christoffel2(c11, c13, c15, c33, c35, c55, rho, k1, k3, qv(i, :))

        end do

        qp = maxval(qv)
        qs = minval(qv)

    end subroutine

    !
    !> Compute qP- and qS-wave phase velocities in the range of polar angle [0, pi]
    !> and azimuth angle [0, pi] for a 3D anisotropic medium
    !
    subroutine min_max_phase_velocity_3d(c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho, n, m, qp, qs)

        ! Parameters
        real, intent(in) :: c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho
        integer, intent(in) :: n, m
        real, intent(inout) :: qp, qs

        integer :: i, j
        real :: theta, dtheta, phi, dphi, k1, k2, k3
        real, allocatable, dimension(:, :, :) :: qv

        qv = zeros(n, m, 3)

        dtheta = const_pi/n
        dphi = const_pi/m

        ! Iteration through different directions
        do j = 1, m
            do i = 1, n

                theta = (i - 1)*dtheta
                phi = (j - 1)*dphi

                k1 = sin(theta)*cos(phi)
                k2 = sin(theta)*sin(phi)
                k3 = cos(theta)

                ! Take the minimum and maximum phase velocities
                ! as the qS- and qP-wave phase velocities, respectively
                call christoffel3(c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, &
                    c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho, k1, k2, k3, &
                    qv(i, j, :))

            end do
        end do

        qp = maxval(qv)
        qs = minval(qv)

    end subroutine

end module
