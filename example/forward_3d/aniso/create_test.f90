
program test

    use libflit

    integer :: n1, n2, n3
    real, allocatable, dimension(:, :, :) :: c11, c12, c13, c22, c23, c33, c44, c55, c66
    real, allocatable, dimension(:, :, :) :: vp, vs, rho, eps, del, gam, the

    call make_directory('./model')

    n1 = 50
    n2 = 60
    n3 = 70

    ! Cij
    c11 = zeros(n1, n2, n3) + 6e9
    c33 = zeros(n1, n2, n3) + 20e9
    c44 = zeros(n1, n2, n3) + 2e9
    c13 = zeros(n1, n2, n3) + 7.5e9

    c22 = c11
    c55 = c44
    c66 = zeros(n1, n2, n3) + 3e9

    c12 = zeros(n1, n2, n3) - 2e9
    c23 = c13

    rho = zeros(n1, n2, n3) + 1.0e3

    call output_array(c11, './model/c11.bin')
    call output_array(c12, './model/c12.bin')
    call output_array(c13, './model/c13.bin')
    call output_array(c22, './model/c22.bin')
    call output_array(c23, './model/c23.bin')
    call output_array(c33, './model/c33.bin')
    call output_array(c44, './model/c44.bin')
    call output_array(c55, './model/c55.bin')
    call output_array(c66, './model/c66.bin')
    call output_array(rho, './model/rho.bin')


    ! Thomsen
    vp = zeros(n1, n2, n3) + 3000.0
    vs = vp/sqrt(3.0)
    rho = zeros(n1, n2, n3) + 1.0e3
    eps = zeros(n1, n2, n3) + 1.0
    del = zeros(n1, n2, n3) + 1.7
    gam = zeros(n1, n2, n3) + 0.4
    the = zeros(n1, n2, n3) + const_pi_half

    call output_array(vp , './model/vp.bin')
    call output_array(vs , './model/vs.bin')
    call output_array(rho, './model/rho.bin')
    call output_array(eps, './model/eps.bin')
    call output_array(del, './model/del.bin')
    call output_array(gam, './model/gam.bin')
    call output_array(the, './model/the.bin')

    ! Geometry
    call make_directory('./geometry')
    open (3, file='./geometry/geometry.txt')
    write (3, *) 'shot_1_geometry.txt'
    close (3)

    open (3, file='./geometry/shot_1_geometry.txt')
    write (3, *) 1
    write (3, *)
    write (3, *) 1
    write (3, *) 200.0, 200.0, 100.0
    write(3, *) 'explosion'
    write (3, *) 'ricker', 20, 1e6, 0
    write (3, *) 0, 0
    write (3, *)
    write (3, *) 1
    write (3, '(3es, es)') 10.0, 10.0, 200.0, 1.0
    close (3)

end program
