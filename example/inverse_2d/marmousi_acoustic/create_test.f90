
program test

    use libflit

    implicit none

    real, allocatable, dimension(:, :) :: w, ww

    call make_directory('./model')

    w = zeros(221, 801)
    call input_array(w, './marmousi_vp.bin')

    w = interp_to(w, [111, 401])

    call output_array(w, './model/vp.bin')

    ww = w

    w = 1.0/gauss_filt(1.0/w, [10.0, 10.0])
    w = 1.0/gauss_filt(1.0/w, [10.0, 10.0])
    w(1:5, :) = ww(1:5, :)
    call output_array(w, './model/vp_init.bin')

    w = 1.0
    w(1:5, :) = 0.0
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
        write(33, *) 'explosion'
        write(33, *) 'gaussian_deriv', 10.0, 1e6, 0.0
        write(33, *) 0, 0
        write(33, *)
        write(33, *) 401
        do j = 1, 401
            write(33, *) (j - 1)*20.0 + 10.0, 0.0, 20.0, 1.0
        end do
        close(33)

    end do
    close(3)

end program


