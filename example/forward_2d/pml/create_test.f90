

program test

    use libflit

    implicit none

    integer :: i

    real, allocatable, dimension(:, :) :: c11, c13, c15, c33, c35, c55, rho
    integer :: n1, n2

    real :: f0 = 12.5
    real :: sx = 50.0
    real :: sz = 2.5
    real :: rz = 0.0

    n1 = 300
    n2 = 300

    call getpar_float('f0', f0, 0.0, required=.true.)
    call getpar_float('sx', sx, 400.0)
    call getpar_float('sz', sz, 0.0, required=.true.)
    call getpar_float('rz', rz, 0.0)

    call make_directory('./model')

    rho = zeros(n1, n2) + 2000.0
    c11 = zeros(n1, n2) + 4e9
    c13 = zeros(n1, n2) + 7.5e9
    c15 = zeros(n1, n2) + 0
    c33 = zeros(n1, n2) + 20e9
    c35 = zeros(n1, n2) + 0
    c55 = zeros(n1, n2) + 2e9

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

end program
