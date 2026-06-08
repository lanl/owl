

program test

    use libflit
    use mod_anisotropy, only: rotate_cij, thomsen_to_cij

    implicit none

    integer :: i, j

    real, allocatable, dimension(:, :) :: vp, vs, eps, del, the, c11, c13, c15, c33, c35, c55, rho
    integer :: n1, n2

    real :: f0 = 12.5
    real :: sx = 50.0
    real :: sz = 2.5
    real :: rz = 0.0

    call make_directory('./model')

    n1 = 200
    n2 = 300

    call getpar_float('f0', f0, 0.0, required=.true.)
    call getpar_float('sx', sx, 400.0)
    call getpar_float('sz', sz, 0.0, required=.true.)
    call getpar_float('rz', rz, 0.0)

    vp = zeros(n1, n2)
    vs = zeros(n1, n2)
    eps = zeros(n1, n2)
    del = zeros(n1, n2)
    the = zeros(n1, n2)
    rho = zeros(n1, n2) + 2000.0

    !$omp parallel do private(i, j)
    do j = 1, n2
        do i = 1, n1

            if (j <= nint(n2/2.0)) then
                vp(i, j) = 2500.0
                vs(i, j) = 1800.0
                eps(i, j) = 0.3
                del(i, j) = -0.2
                the(i, j) = -real(40.0*const_deg2rad)
            else
                vp(i, j) = 3000.0
                vs(i, j) = 2000.0
                eps(i, j) = 0.0
                del(i, j) = 0.0
                the(i, j) = 0.0
            end if

        end do
    end do
    !$omp end parallel do

    c11 = zeros_like(vp)
    c13 = zeros_like(vp)
    c15 = zeros_like(vp)
    c33 = zeros_like(vp)
    c35 = zeros_like(vp)
    c55 = zeros_like(vp)

    call thomsen_to_cij(vp, vs, rho, eps, del, -the, c11, c13, c15, c33, c35, c55)

    call output_array(vp, './model/vp.bin')
    call output_array(vs, './model/vs.bin')
    call output_array(eps, './model/eps.bin')
    call output_array(del, './model/del.bin')
    call output_array(the, './model/the.bin')

    call output_array(c11, './model/c11.bin')
    call output_array(c13, './model/c13.bin')
    call output_array(c15, './model/c15.bin')
    call output_array(c33, './model/c33.bin')
    call output_array(c35, './model/c35.bin')
    call output_array(c55, './model/c55.bin')
    call output_array(rho, './model/rho.bin')

    ! Geometry
    call make_directory('./geometry')
    open(3, file='./geometry/geometry.txt')
    write(3, *) 'shot_1_geometry.txt'
    close(3)

    open(3, file='./geometry/shot_1_geometry.txt')
    write(3, *) 1
    write(3, *)
    write(3, *) 1
    write(3, *) sx, 0.0, sz
    write(3, *) 'explosion'
    write(3, *) 'ricker', f0, 1e6, 0.0
    write(3, *) 0, 0

    write(3, *)

    write(3, *) n2
    do i = 1, n2
        write(3, '(3es, es)') (i - 1)*10.0, 0.0, rz, 1.0
    end do

    close(3)

    open(3, file='./geometry/rec.txt')
    do i = 1, n2
        write(3, '(3es, es)') (i - 1)*10.0, 0.0, rz, 1.0
    end do
    close(3)

    !contains
    !
    !    elemental subroutine cij_rotate_2d(c11, c13, c15, c33, c35, c55, thetay)
    !
    !        real, intent(inout) :: c11, c13, c15, c33, c35, c55
    !        real, intent(in) :: thetay
    !
    !        real, dimension(1:6, 1:6) :: c, a
    !        real :: c12, c14, c16
    !        real :: c22, c23, c24, c25, c26
    !        real :: c34, c36
    !        real :: c44, c45, c46
    !        real :: c56
    !        real :: c66
    !
    !        c12 = 0.0
    !        c14 = 0.0
    !        c16 = 0.0
    !        c22 = 0.0
    !        c23 = 0.0
    !        c24 = 0.0
    !        c25 = 0.0
    !        c26 = 0.0
    !        c34 = 0.0
    !        c36 = 0.0
    !        c44 = 0.0
    !        c45 = 0.0
    !        c46 = 0.0
    !        c56 = 0.0
    !        c66 = 0.0
    !
    !        ! Rotation matrix; For 2D case, the rotation order does not matter
    !        a = bond_matrix(rotation_matrix(thetay, 'y'))
    !
    !        ! Elasticity matrix
    !        c = reshape([ &
        !            c11, c12, c13, c14, c15, c16, &
        !            c12, c22, c23, c24, c25, c26, &
        !            c13, c23, c33, c34, c35, c36, &
        !            c14, c24, c34, c44, c45, c46, &
        !            c15, c25, c35, c45, c55, c56, &
        !            c16, c26, c36, c46, c56, c66], [6, 6])
    !
    !        ! Rotation
    !        c = matmul(a, matmul(c, transpose(a)))
    !
    !        ! Take relevant elements for 2D scenario
    !        c11 = c(1, 1)
    !        c13 = c(1, 3)
    !        c15 = c(1, 5)
    !        c33 = c(3, 3)
    !        c35 = c(3, 5)
    !        c55 = c(5, 5)
    !
    !    end subroutine cij_rotate_2d
    !
    !    pure function bond_matrix(r) result(b)
    !
    !        real, dimension(1:3, 1:3), intent(in) :: r
    !        real, allocatable, dimension(:, :) :: b
    !
    !        allocate (b(1:6, 1:6))
    !        b(1, :) = [ &
        !            r(1, 1)**2, &
        !            r(1, 2)**2, &
        !            r(1, 3)**2, &
        !            2.0*r(1, 2)*r(1, 3), &
        !            2.0*r(1, 1)*r(1, 3), &
        !            2.0*r(1, 1)*r(1, 2)]
    !        b(2, :) = [ &
        !            r(2, 1)**2, &
        !            r(2, 2)**2, &
        !            r(2, 3)**2, &
        !            2.0*r(2, 2)*r(2, 3), &
        !            2.0*r(2, 1)*r(2, 3), &
        !            2.0*r(2, 1)*r(2, 2)]
    !        b(3, :) = [ &
        !            r(3, 1)**2, &
        !            r(3, 2)**2, &
        !            r(3, 3)**2, &
        !            2.0*r(3, 2)*r(3, 3), &
        !            2.0*r(3, 1)*r(3, 3), &
        !            2.0*r(3, 1)*r(3, 2)]
    !        b(4, :) = [ &
        !            r(2, 1)*r(3, 1), &
        !            r(2, 2)*r(3, 2), &
        !            r(2, 3)*r(3, 3), &
        !            r(2, 2)*r(3, 3) + r(2, 3)*r(3, 2), &
        !            r(2, 1)*r(3, 3) + r(2, 3)*r(3, 1), &
        !            r(2, 1)*r(3, 2) + r(2, 2)*r(3, 1)]
    !        b(5, :) = [ &
        !            r(1, 1)*r(3, 1), &
        !            r(1, 2)*r(3, 2), &
        !            r(1, 3)*r(3, 3), &
        !            r(1, 2)*r(3, 3) + r(1, 3)*r(3, 2), &
        !            r(1, 1)*r(3, 3) + r(1, 3)*r(3, 1), &
        !            r(1, 1)*r(3, 2) + r(1, 2)*r(3, 1)]
    !        b(6, :) = [ &
        !            r(1, 1)*r(2, 1), &
        !            r(1, 2)*r(2, 2), &
        !            r(1, 3)*r(2, 3), &
        !            r(1, 2)*r(2, 3) + r(1, 3)*r(2, 2), &
        !            r(1, 1)*r(2, 3) + r(1, 3)*r(2, 1), &
        !            r(1, 1)*r(2, 2) + r(1, 2)*r(2, 1)]
    !
    !    end function bond_matrix
    !
    !    elemental subroutine thomsen_to_cij2d( &
        !            vp, vs, rho, eps, del, the, c11, c13, c15, c33, c35, c55)
    !
    !        real, intent(in) :: vp, vs, rho, eps, del, the
    !        real, intent(out) :: c11, c13, c15, c33, c35, c55
    !
    !        real :: tmp
    !
    !        c33 = vp**2*rho
    !        c55 = vs**2*rho
    !        c11 = c33*(1 + 2*eps)
    !        tmp = 2*c33*(c33 - c55)*del + (c33 - c55)**2
    !        c13 = sqrt(abs(tmp)) - c55
    !        c15 = 0.0
    !        c35 = 0.0
    !
    !        if (the /= 0) then
    !            ! Clockwise rotation
    !            call cij_rotate_2d(c11, c13, c15, c33, c35, c55, -the)
    !        end if
    !
    !    end subroutine thomsen_to_cij2d

end program test
