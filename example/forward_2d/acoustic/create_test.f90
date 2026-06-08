

program test

    use libflit
    use librgm

    type(rgm2_curved) :: p
    integer :: n1, n2
    real :: d1, d2

    real :: f0 = 25
    real :: sx = 400.0
    real :: sz = 500.0
    real :: rz = 0.0

    n1 = 100
    n2 = 200
    d1 = 10.0
    d2 = 10.0

    ! Velocity and density models
    call make_directory('./model')

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

    p%vp = rescale(p%vp, [2000.0, 4500.0])
    p%rho = rescale(p%rho, [1500.0, 3500.0])

    call output_array(p%vp, './model/vp.bin')
    call output_array(p%rho, './model/rho.bin')

    ! Optionally, these models can be saved to hdf5
    call fh5_open('./model/model.h5', fh5_fid, mode='w')
    call fh5_write_attr(fh5_fid, '/', 'nx', n2)
    call fh5_write_attr(fh5_fid, '/', 'nz', n1)
    call fh5_write_attr(fh5_fid, '/', 'dx (m)', d2)
    call fh5_write_attr(fh5_fid, '/', 'dz (m)', d1)
    call fh5_write_attr(fh5_fid, '/', 'origin x (m)', 0.0)
    call fh5_write_attr(fh5_fid, '/', 'origin z (m)', 0.0)
    call fh5_write(fh5_fid, '/vp', p%vp)
    call fh5_write_attr(fh5_fid, '/vp', 'units', 'm/s')
    call fh5_write(fh5_fid, '/rho', p%rho)
    call fh5_write_attr(fh5_fid, '/rho', 'units', 'kg/m^3')
    call fh5_close(fh5_fid)

    ! Source-receiver geometry
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

    write(3, *) p%n2
    do i = 1, p%n2
        write(3, '(3es, es)') (i - 1)*10.0, 0.0, rz, 1.0
    end do


    close(3)

end program test
