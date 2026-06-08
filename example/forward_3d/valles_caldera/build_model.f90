
program test

    use libflit
    use librgm

    implicit none

    integer :: n1, n2, n3, i, j, k
    real :: d1, d2, d3
    real, allocatable, dimension(:, :, :) :: vp, rho
    real :: hmax
    real, allocatable, dimension(:, :) :: topo, fl
    type(rgm3_curved) :: p

    n1 = 121
    n2 = 640
    n3 = 818

    d1 = 50
    d2 = 50
    d3 = 50

    p%n1 = n1
    p%n2 = n2
    p%n3 = n3

    p%nl = 30
    p%refl_shape = 'gaussian'
    p%ng = 1
    p%seed = 1234
    p%refl_mu2 = [320.0, 320.0]
    p%refl_mu3 = [360.0, 360.0]
    p%refl_sigma2 = [150.0, 150.0]*2
    p%refl_sigma3 = [150.0, 150.0]*2
    p%refl_height = [0, 121]
    p%refl_height_top = [0, 2]
    p%nf = 10
    p%delta_dip = 0
    p%disp = [2.0, 7.0]
    p%yn_salt = .true.
    p%nsalt = 1
    p%salt_top_z = [0.0, 0.01]
    p%salt_radius = [180, 210]/2.5*2

    call p%generate()

    p%vp = rescale(p%vp, [3000.0, 6000.0])
    call output_array(p%vp, './model/vp.bin')

    p%vp = rescale(p%vp, [2000.0, 3500.0])
    call output_array(p%vp, './model/vs.bin')

    !	rho = load('./vp.bin', n1, n2, n3)
    !	rho = rescale(rho**2, [2000.0, 4000.0])
    !	call output_array(rho, './rho.bin')

    topo = load('./model/dem_utm.bin', 666, 844)
    topo = topo(14:14 + n2 - 1, 13:13 + n3 - 1)

    print *, shape(topo)

    topo = flip(topo, [1])
    call output_array(topo, './model/topo.bin')


    fl = topo
    vp = ones(n1, n2, n3)
    rho = load('./model/vp.bin', n1, n2, n3)
    hmax = maxval(topo)
    !$omp parallel do private(i, j, k) collapse(3)
    do k = 1, n3
        do j = 1, n2
            do i = 1, n1
                if ((i - 1)*d1 <= hmax - topo(j, k)) then
                    vp(i, j, k) = 0.0
                    rho(i, j, k) = nan()
                end if
            end do
        end do
    end do
    !$omp end parallel do
    call output_array(vp, './model/mask.bin')
    call output_array(rho, './model/vp_masked.bin')

    topo = pad(topo, [15, 15, 15, 15])
    topo = topo - minval(topo)

    open(3, file='./model/topo.txt')
    do i = 1, size(topo, 2)
        do j = 1, size(topo, 1)
            write(3, '(3es)') (i - 1 - 15)*d3, (j - 1 - 15)*d2, topo(j, i)
        end do
    end do

    call make_directory('./geometry')
    open(3, file='./geometry/geometry.txt')
    write(3, '(a)') 'shot_1_geometry.txt'
    close(3)

    open (3, file='./geometry/shot_1_geometry.txt')
    write (3, *) 1
    write (3, *)
    write (3, *) 1
    write (3, *) 150*50.0, 550*50.0, 500.0
    ! strike, dip, rake = 135, 60, -90
    write (3, *) 'mt 0.43 -0.43 0.0 -0.25 0.55 -0.55 '
    write (3, *) 'ricker', 2.5, 1e9, 0
    write (3, *) 0, 0
    write (3, *)
    write (3, *) 100
    do i = 1, 10
        do j = 1, 10
            write (3, '(3es, es)') (i - 1)*4000.0 + 2000.0, (j - 1)*3200.0 + 1600.0, 0.0, 1.0
        end do
    end do
    close (3)

end program
