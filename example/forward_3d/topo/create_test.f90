
program test

    use libflit
    use librgm

    implicit none

    type(fractal_noise_2d) :: xx
    integer :: nrx, nry
    integer :: n1, n2, n3, pml, i, j
    real, allocatable, dimension(:, :, :) :: vp, vs, rho
    real, allocatable, dimension(:, :) :: tp

    real :: f0 = 12.5
    real :: sz = 0.0
    real :: rz = 0.0

    call make_directory('./model')

    n1 = 40
    n2 = 35
    n3 = 200
    pml = 15

    xx%seed = 112233
    xx%n1 = n2
    xx%n2 = n3
    tp = xx%generate()
    tp = gauss_filt(tp, [4.0, 4.0])
    tp = tp - mean(tp)
    tp = taper(tp, len=[10, 10, 10, 10], method=['blackman', 'blackman', 'blackman', 'blackman'])
    tp = pad(tp, [pml + 1, pml + 1, pml + 1, pml + 1])
    tp = rescale(tp, [0.0, 100.0])

    call output_array(tp(pml + 1 + 1:pml + 1 + n2, pml + 1 + 1:pml + 1 + n3), './topo.bin')

    open (3, file='ftopo.txt')
    do j = -pml, n3 + pml + 1
        do i = -pml, n2 + pml + 1
            write (3, *) (j - 1)*10.0, (i - 1)*10.0, tp(i + pml + 1, j + pml + 1)
        end do
    end do
    close (3)

    vp = zeros(n1, n2, n3) + 3000.0
    vs = zeros(n1, n2, n3) + 3000.0/sqrt(3.0)
    rho = zeros(n1, n2, n3) + 2000.0

    call output_array(vp, './model/vp.bin')
    call output_array(vs, './model/vs.bin')
    call output_array(rho, './model/rho.bin')

    call make_directory('./geometry')
    open (3, file='./geometry/geometry.txt')
    write (3, *) 'shot_1_geometry.txt'
    close (3)

    nrx = floor(n3/2.0)
    nry = floor(n2/2.0)

    open (3, file='./geometry/shot_1_geometry.txt')
    write (3, *) 1
    write (3, *)
    write (3, *) 1
    write (3, *) 200.0, 100.0, sz
    ! Explosion source
    write (3, *) 'explosion'
    ! Force vector source
    !    write(3, *) 'force', 45.0, 45.0
    ! MT source
    !    write(3, '(a, 6es)') 'mt', 1.0, 0.5, 0.5, 0.2, 0.1, -0.2
    write (3, *) 'ricker', f0, 1e6, 0
    write (3, *) 0, 0
    write (3, *)

    write (3, *) nrx*nry

    do j = 1, nry
        do i = 1, nrx
            write (3, '(3es, es)') (i - 1)*20 + 10.0, (j - 1)*20 + 10.0, rz, 1.0
        end do
    end do

    close (3)

    open(3, file='./geometry/rec.txt')
    do j = 1, nry
        do i = 1, nrx
            write (3, '(3es, es)') (i - 1)*20 + 10.0, (j - 1)*20 + 10.0, rz
        end do
    end do
    close(3)


    open(3, file='./geometry/rec_subset.txt')
    do j = 1, nry
        if (j == 6 .or. j == 17) then
            do i = 1, nrx
                write (3, '(3es, es)') (i - 1)*20 + 10.0, (j - 1)*20 + 10.0, rz
            end do
        end if
    end do
    close(3)

    open(3, file='./geometry/src.txt')
    write(3, '(3es)') 200.0, 100.0
    close(3)

end program
