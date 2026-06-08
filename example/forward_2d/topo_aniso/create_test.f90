

program test

    use libflit
    use librgm
    use mod_anisotropy, only: thomsen_to_cij

    implicit none

    type(rgm2_curved) :: p
    real ,allocatable, dimension(:) :: tp
    type(fractal_noise_1d) :: xx
    integer :: i, j

    real, allocatable, dimension(:, :) :: vp, vs, eps, del, the, c11, c13, c15, c33, c35, c55, rho
    real, allocatable, dimension(:, :) :: yc
    real, allocatable, dimension(:, :, :) :: w
    integer :: n1, n2

    real :: f0 = 12.5
    real :: sx = 50.0
    real :: sz = 2.5
    real :: rz = 0.0
    integer :: pml = 15

    n1 = 50
    n2 = 300

    call getpar_float('f0', f0, 0.0, required=.true.)
    call getpar_float('sx', sx, 400.0)
    call getpar_float('sz', sz, 0.0, required=.true.)
    call getpar_float('rz', rz, 0.0)

    xx%n1 = n2 + 2*pml + 2
    xx%seed = 567

    tp = xx%generate()
    tp = gauss_filt(tp, 2.0)
    tp = rescale(tp, [0.0, 100.0])

    open(3, file='ftopo.txt')
    do i = -pml, n2 + pml + 1
        write(3, *) (i - 1)*10.0, tp(i + pml + 1)
    end do
    close(3)

    call make_directory('./model')

    ! Model 1
    p%n1 = n1
    p%n2 = n2
    p%refl_shape = 'cauchy'
    p%refl_shape_top = 'perlin'
    p%refl_smooth_top = 2
    p%refl_sigma2 = [90.0, 120.0]
    p%seed = 223
    p%ng = 2
    p%refl_sigma2 = [30.0, 50.0]
    p%lwv = 0.3
    p%lwh = 0.4
    p%refl_height = [0, 40]
    p%refl_height_top = [0, 2]
    p%nl = 15
    p%disp = [2.0, 5.0]
    p%yn_elastic = .true.
    p%nf = 10
    call p%generate

    vp = rescale(p%vp, [2500.0, 4500.0])
    vs = rescale(p%vp, [1600.0, 3000.0])
    eps = rescale(1.0/p%vp, [0.2, 0.3])
    del = rescale((1.0/p%vp)**2, [-0.1, 0.2])
    rho = rescale(p%vp**2, [2000.0, 3500.0])

    p%nf = 0
    call p%generate

    call local_dip(p%vp, the)
    the = gauss_filt(the, [3.0, 3.0])

    call output_array(vp, './model/vp.bin')
    call output_array(vs, './model/vs.bin')
    call output_array(eps, './model/eps.bin')
    call output_array(del, './model/del.bin')
    call output_array(the, './model/the.bin')


    c11 = zeros_like(vp)
    c13 = zeros_like(vp)
    c15 = zeros_like(vp)
    c33 = zeros_like(vp)
    c35 = zeros_like(vp)
    c55 = zeros_like(vp)

    ! note that specfem uses a different theta convention
    call thomsen_to_cij(vp, vs, rho, eps, del, -the, c11, c13, c15, c33, c35, c55)

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


    w = zeros(n1, n2, 6)
    yc = zeros(n1, n2)
    !$omp parallel do private(j)
    do j = 1, n2
        yc(:, j) = linspace(-tp(j + pml + 1), 400.0, n1)
        w(:, j, 1) = ginterp(linspace(-maxval(tp), 400.0, n1), vp(:, j), yc(:, j), 'linear')
        w(:, j, 2) = ginterp(linspace(-maxval(tp), 400.0, n1), vs(:, j), yc(:, j), 'linear')
        w(:, j, 3) = ginterp(linspace(-maxval(tp), 400.0, n1), rho(:, j), yc(:, j), 'linear')
        w(:, j, 4) = ginterp(linspace(-maxval(tp), 400.0, n1), eps(:, j), yc(:, j), 'linear')
        w(:, j, 5) = ginterp(linspace(-maxval(tp), 400.0, n1), del(:, j), yc(:, j), 'linear')
        w(:, j, 6) = ginterp(linspace(-maxval(tp), 400.0, n1), the(:, j), yc(:, j), 'linear')
    end do
    !$omp end parallel do
    open(3, file='./model/model.txt')
    do j = 1, n2
        do i = 1, n1
            write(3, '(8es)') (j - 1)*10.0, yc(i, j), w(i, j, :)
        end do
    end do
    close(3)

end program test
