
program test

    use libflit
    use librgm

    implicit none

    integer :: i, j, ishotx, ishoty, l
    real :: minv
    real, allocatable, dimension(:, :, :) :: v, m

    type(rgm3_curved) :: p
    integer :: n1, n2, n3
    real :: f0 = 15.0
    real :: sz = 0.0
    real :: rz = 0.0

    n1 = 81
    n2 = 201
    n3 = 221

    p%n1 = n1
    p%n2 = n2
    p%n3 = n3

    p%seed = 5678
    p%refl_shape = 'perlin'
    p%refl_shape_top = 'perlin'
    p%refl_smooth = 3
    p%refl_smooth_top = 0
    p%lwv = 0.5
    p%lwh = 0.6
    p%refl_height = [0, 80]
    p%refl_height_top = [0, 2]
    p%unconf = 1
    p%unconf_z = [0.1, 0.2]
    p%unconf_nl = 16
    p%unconf_refl_height = [0, 0]
    p%nl = 15
    p%disp = [6.0, 15.0]
    p%yn_elastic = .true.
    p%delta_v = 750
    p%nf = 7

    call p%generate

    p%vp = rescale(p%vp, [2200.0, 4500.0])
    minv = minval(p%vp)

    v = p%vp
    m = ones_like(p%vp)

    call make_directory('./model')

    p%vp (1:5, :, :) = 2100.0
    m(1:5, :, :) = 0.0

    call output_array(p%vp, './model/vp.bin')
    call output_array(m, './model/mask.bin')

    p%vp = gauss_filt(v, [4.0, 13.0, 13.0])
    p%vp(1:5, :, :) = 2100.0
    call output_array(p%vp, './model/vp_init.bin')

    ! Geometry
    call make_directory('./geometry')
    open (3, file='./geometry/geometry.txt')
    do i = 1, 110
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'
    end do

    close (3)

    l = 1
    do ishotx = 1, 11
        do ishoty = 1, 10

            open (3, file='./geometry/shot_'//num2str(l)//'_geometry.txt')
            write (3, *) l
            write (3, *)
            write (3, *) 1
            write (3, *) 200.0 + (ishotx - 1)*400.0, 200.0 + (ishoty - 1)*400.0, 0.0, sz
            write (3, *) 'explosion'
            write (3, *) 'ricker', f0, 1e6, 0.0
            write (3, *) 0, 0

            write (3, *)

            write (3, *) floor(n2/4.0) * floor(n3/4.0)
            do i = 1, floor(n3/4.0)
                do j = 1, floor(n2/4.0)
                    write (3, '(3es, es)') 40 + (i - 1)*80.0,  40 + (j - 1)*80.0, rz, 1.0
                end do
            end do

            close (3)

            l = l + 1

        end do
    end do

end program
