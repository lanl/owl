
program test

    use libflit
    use mod_su

    implicit none

    integer :: i, j
    type(su) :: w
    real, allocatable, dimension(:, :) :: vp
    real, allocatable, dimension(:) :: h
    integer :: nt
    real :: f0

    ! Initial velocity model
    !    raw_file = vp_raw.bin
    !    dtype = float32
    !    endian = little
    !    order = column-based
    !    array_shape = (nz, nx)
    !    nx = 3820
    !    nz = 1200
    !    dx = 125000.0
    !    dz = 5.0

    call make_directory('./model')
    vp = load('./vp_raw.bin', 1200, 3820)
    vp = interp_to(vp, [480, 3820])
    call output_array(vp, './model/vp.bin')

    ! Geometry
    call make_directory('./geometry')

    nt = count_nonempty_lines('./Wavelet.txt')
    h = load('./Wavelet.txt', nt, ascii=.true.)
    call output_array(h, './geometry/chevron_wavelet.bin')

    f0 = 10.0

    open (6, file='../geometry/chevron_wavelet.txt')
    do i = 1, nt
        write (6, '(2es)') (i - 1)*0.6666666e-3, h(i)
    end do
    close(6)

    open (6, file='../geometry/chevron_geometry.txt', status='replace')

    do i = 1, 1600

        write(6, *) 'shot_'//num2str(i)//'_geometry.txt'

        call w%load('./csg/shot_'//tidy(num2str(i, '(i)'))//'_seismogram_p.su')

        open(11, file='../geometry/shot_'//num2str(i)//'_geometry.txt')

        write (11, *) i
        write (11, *)

        write (11, *) 1
        write (11, '(3es)') real(w%trace(1)%header%SourceX*1.0e-4), 0.0, 15.0
        write (11, *) 'explosion'
        write (11, '(x, a, 3es)') 'custom', f0, 6.0e-3, 0.0175
        write (11, *) './geometry/chevron_wavelet.txt'
        write (11, *) 0, 0
        write (11, *)

        write (11, *) w%nr
        do j = 1, w%nr
            write (11, '(4es)') real(w%trace(j)%header%GroupX*1.0e-4), 0.0, 15.0, 1.0
        end do

        close (11)

        print *, i

    end do
    close (6)

end program
