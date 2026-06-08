
program test

    use libflit
    use librgm

    implicit none

    type(rgm2_curved) :: p
    real, allocatable, dimension(:) :: tp
    type(fractal_noise_1d) :: xx
    real, allocatable, dimension(:, :) :: vp, vs, rho, m
    real :: minvp, minvs, minrho
    integer :: n1, n2, i, j, ishot

    real :: f0 = 7.5
    real :: sz = 0.0
    real :: rz = 0.0
    integer :: pml = 15

    call make_directory('./model')

    n1 = 101
    n2 = 301

    ! Topography
    xx%n1 = n2
    xx%seed = 1122
    xx%octaves = 7
    xx%periods1 = 7

    tp = xx%generate()
    call pad_array(tp, [pml + 1, pml + 1])

    tp = gauss_filt(tp, 1.0)
    tp = rescale(tp, [0.0, 300.0])

    open (3, file='./model/ftopo.txt')
    do i = -pml, n2 + pml + 1
        write (3, *) (i - 1)*10.0, tp(i)
    end do
    close (3)

    ! Mask
    m = ones(n1, n2)
    do j = 1, n2
        do i = 1, n1
            if (i - 1.0 < (maxval(tp) - tp(j))/10.0) then
                m(i, j) = 0.0
            end if
        end do
    end do
    call output_array(m, './model/mask.bin')

    ! Models
    p%n1 = n1
    p%n2 = n2
    p%refl_shape = 'perlin'
    p%refl_shape_top = 'perlin'
    p%refl_smooth = 3
    p%refl_smooth_top = 0
    p%seed = 2341
    p%lwv = 0.5
    p%lwh = 0.5
    p%refl_height = [0, 80]
    p%refl_height_top = [0, 2]
    p%nl = 15
    p%disp = [20.0, 40.0]
    p%yn_elastic = .true.
    p%nf = 6
    call p%generate

    p%vp = rescale(p%vp, [2500.0, 4500.0])
    p%vs = rescale(p%vs, [1600.0, 2500.0])
    p%rho = 1.0/gauss_filt(1.0/p%rho, [6.0, 20.0])

    vp = 1.0/gauss_filt(1.0/p%vp, [6.0, 20.0])
    vs = 1.0/gauss_filt(1.0/p%vs, [6.0, 20.0])
    rho = p%rho

    minvp = minval(p%vp, mask=(m ==1))
    minvs = minval(p%vs, mask=(m ==1))
    minrho = minval(p%rho, mask=(m ==1))

    where (m == 0)
        p%vp = minvp
        p%vs = minvs
        p%rho = minrho
        vp = minvp
        vs = minvs
        rho = minrho
    end where

    call output_array(p%vp/p%vs, './model/vpvs_ratio.bin')

    call output_array(p%vp, './model/vp.bin')
    call output_array(p%vs, './model/vs.bin')
    call output_array(p%rho, './model/rho.bin')

    call output_array(vp, './model/vp_init.bin')
    call output_array(vs, './model/vs_init.bin')
    call output_array(rho, './model/rho_init.bin')

    ! Geometry
    call make_directory('./geometry')
    open (3, file='./geometry/geometry.txt')
    do i = 1, 60
        write (3, *) 'shot_'//num2str(i)//'_geometry.txt'
    end do

    close (3)

    do ishot = 1, 60
        open (3, file='./geometry/shot_'//num2str(ishot)//'_geometry.txt')
        write (3, *) ishot
        write (3, *)
        write (3, *) 1
        write (3, *) 20.0 + (ishot - 1)*50.0, 0.0, sz
        write (3, *) 'explosion'
        write (3, *) 'ricker', f0, 1e6, 0.0
        write (3, *) 0, 0

        write (3, *)

        write (3, *) p%n2
        do i = 1, p%n2
            write (3, '(3es, es)') (i - 1)*10.0, 0.0, rz, 1.0
        end do

        close (3)
    end do

end program
