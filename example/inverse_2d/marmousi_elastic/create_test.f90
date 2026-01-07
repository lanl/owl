
program test

    use libflit

    real, allocatable, dimension(:, :) :: w, ww, vp, vs
    integer :: nsf

    w = zeros(221, 801)
    call input_array(w, './model/marmousi_vp.bin')

    w = interp_to(w, [111, 401])

    nsf = 5

    vp = w
    vs = rescale(vp, [680.0, 2700.0])

    eps = rescale(median_filt(1.0/vp), [0.2, 0.3])
    del = rescale(median_filt(1.0/vs), [0.1, 0.15])

    the = zeros_like(vp)
    the = gauss_filt(deriv(vp, dim=1), [1.0, 1.0])
    call gstdip(the, the)
    the = gauss_filt(median_filt(the, [4, 4]), [2.0, 2.0])


    call output_array(vp, './model/vp.bin')
    call output_array(vs, './model/vs.bin')

    call output_array(eps, './model/eps.bin')
    call output_array(del, './model/del.bin')
    call output_array(the, './model/the.bin')

    ww = vp
    vp = 1.0/gauss_filt(1.0/vp, [8.0, 8.0])
    vp = 1.0/gauss_filt(1.0/vp, [8.0, 8.0])
    vp(1:nsf, :) = ww(1:nsf, :)
    call output_array(vp, './model/vp_init.bin')

    ww = vs
    vs = 1.0/gauss_filt(1.0/vs, [8.0, 8.0])
    vs = 1.0/gauss_filt(1.0/vs, [8.0, 8.0])
    vs(1:nsf, :) = ww(1:nsf, :)
    call output_array(vs, './model/vs_init.bin')

    w = 1.0
    w(1:nsf, :) = 0.0
    call output_array(w, './model/mask.bin')

    ! Geometry
    call make_directory('./geometry')

    open(3, file='./geometry/geometry.txt')
    do i = 1, 80

        write(3, *) 'shot_'//num2str(i)//'_geometry.txt'

        open(33, file='./geometry/shot_'//num2str(i)//'_geometry.txt')
        write(33, *) i
        write(33, *)
        write(33, *) 1
        write(33, *) (i - 1)*100.0 + 50.0, 0.0, 20.0
        !        write(33, *) 'force', 0.0, 0.0
        write(33, *) 'explosion'
        write(33, *) 'ricker', 5.0, 100.0, 0.0
        write(33, *) 0, 0
        write(33, *)
        write(33, *) 401
        do j = 1, 401
            write(33, *) (j - 1)*20.0 + 10.0, 0.0, 20.0, 1.0
        end do
        close(33)

    end do
    close(3)

end program test


