
program test

    use libflit

    implicit none

    real, allocatable, dimension(:, :) :: w

    ! Geometry
    call make_directory('./geometry')
    open(3, file='./geometry/geometry.txt')
    write(3, *) 'shot_1_geometry.txt'
    close(3)

    open(33, file='./geometry/shot_1_geometry.txt')
    write(33, *) 1
    write(33, *)
    write(33, *) 1
    write(33, *) 300.0, 0.0, 1000.0

    write(33, *) 'explosion'
    write(33, *) 'gaussian', 20.0, 1e5, 0.0
    write(33, *) 0, 0
    write(33, *)

    write(33, *) 1
    write(33, *) 1700.0, 0.0, 1000.0, 1.0
    close(33)

    ! Model
    call make_directory('./model')

    w = zeros(201, 201) + 2500.0
    call output_array(w, './model/vp_init.bin')

end program test
