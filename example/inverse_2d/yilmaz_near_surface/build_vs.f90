
program test

    use libflit

    implicit none

    real, allocatable, dimension(:, :) :: vp
    integer :: n1, n2

    n1 = 21
    n2 = 97

    vp = load('./fatt/iteration_100/model/updated_vp.bin', n1, n2)

    vp = interp(vp, shape(vp), [1.0, 1.0], [0.0, 0.0], &
        [61, 291], [0.333333, 0.333333], [0.0, 0.0])
    call output_array(vp, './model/vp_fatt.bin')

    vp = rescale(vp, [200.0, 4000.0])
    call output_array(vp, './model/vs_fatt.bin')

end program

