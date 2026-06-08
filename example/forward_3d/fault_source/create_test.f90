
program test

    use libflit
    use librgm

    implicit none

    type(rgm3_curved) :: p
    integer :: n1, n2, n3
    real, allocatable, dimension(:, :, :) :: vp, vs, rho
    real, allocatable, dimension(:, :) :: xyzt0
    integer :: np
    integer :: i, j
    real :: drx, dry

    call make_directory('./model')

    n1 = 81
    n2 = 201
    n3 = 201

    p%n1 = n1
    p%n2 = n2
    p%n3 = n3

    p%nl = 20
    p%refl_shape = 'cauchy'
    p%ng = 3
    p%seed = 1234
    p%refl_mu2 = [0.0, 200.0]
    p%refl_mu3 = [0.0, 200.0]
    p%refl_sigma2 = [80.0, 150.0]
    p%refl_sigma3 = [40.0, 150.0]
    p%refl_height = [0, 50]
    p%refl_shape_top = 'perlin'
    p%refl_height_top = [0, 10]
    p%lwv = 0.5
    p%nf = 0

    call p%generate

    vp = rescale(p%vp, [2000.0, 4500.0])
    vs = rescale(p%vp, [1500.0, 2500.0])
    rho = p%rho

    call output_array(vp , './model/vp.bin')
    call output_array(vs , './model/vs.bin')
    call output_array(rho, './model/rho.bin')

    ! Geometry
    call make_directory('./geometry')
    open (3, file='./geometry/geometry.txt')
    write (3, *) 'shot_1_geometry.txt'
    close (3)

    np = count_nonempty_lines('./fault_source_xyzt0.txt')
    xyzt0 = load('./fault_source_xyzt0.txt', np, 4, ascii=.true.)

    drx = 200.0
    dry = 200.0

    open (3, file='./geometry/shot_1_geometry.txt')
    write (3, *) 1
    write (3, *)
    write (3, *) np
    do i = 1, np
        write (3, '(3es)') xyzt0(i, 1:3)*1.0e3
        write (3, '(a, 6es)') 'mt', 0.4330127, 0.4330127, -0.8660254, -0.4330127, 0.35355339, -0.35355339
        write (3, '(a, 3es)') 'gaussian', 1.0, 1.0e6, xyzt0(i, 4)
        write (3, '(2es)') 0, 0
    end do
    write (3, *)
    write (3, *) 51*51
    do i = 1, 51
        do j = 1, 51
            write (3, '(3es, es)') (i - 1.0)*drx, (j - 1.0)*dry, 0.0, 1.0
        end do
    end do
    close (3)

end program
