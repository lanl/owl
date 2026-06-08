
program test

    use libflit
    use librgm

    type(rgm2_curved) :: p
    real ,allocatable, dimension(:) :: tp
    type(fractal_noise_1d) :: xx

    real, allocatable, dimension(:, :) :: vp, vs, rho
    integer :: n1, n2

    real :: f0 = 12.5
    real :: sx = 50.0
    real :: sz = 2.5
    real :: rz = 0.0
    integer :: pml = 15

    n1 = 50
    n2 = 300

    call make_directory('./model')

    call getpar_float('f0', f0, 0.0, required=.true.)
    call getpar_float('sx', sx, 400.0)
    call getpar_float('sz', sz, 0.0, required=.true.)

    xx%n1 = n2 + 2*pml + 2
    xx%seed = 12343331
    xx%seed = 1231212

    tp = xx%generate()
    tp = gauss_filt(tp, 4.0)
    tp = rescale(tp, [0.0, 200.0])

    open(3, file='./model/ftopo.txt')
    do i = -pml, n2 + pml + 1
        write(3, *) (i - 1)*10.0, tp(i + pml + 1)
    end do
    close(3)

    ! Model 1
    p%n1 = n1
    p%n2 = n2
    p%refl_shape = 'gaussian'
    p%refl_shape_top = 'perlin'
    p%refl_smooth_top = 2
    p%refl_sigma2 = [90.0, 120.0]
    p%seed = 1235
    p%ng = 4
    p%refl_sigma2 = [30.0, 50.0]
    p%lwv = 0.3
    p%lwh = 0.4
    p%refl_height = [0, 80]
    p%refl_height_top = [0, 2]
    p%nl = 15
    p%disp = [5.0, 10.0]
    p%yn_elastic = .true.
    p%nf = 5
    call p%generate

    vp = rescale(p%vp, [3000.0, 4500.0])
    vs = rescale(p%vs, [1732.0, 3000.0])
    rho = p%rho

    call output_array(vp, './model/vp.bin')
    call output_array(vs, './model/vs.bin')
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
        write(3, '(3es, es)') (i - 1)*10.0, 0.0, rz + 0, 1.0
    end do

    close(3)

    open(3, file='./geometry/rec.txt')
    do i = 1, n2
        write(3, '(3es, es)') (i - 1)*10.0, 0.0, rz + 0, 1.0
    end do
    close(3)

end program
