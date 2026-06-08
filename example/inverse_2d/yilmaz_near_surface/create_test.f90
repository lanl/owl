
program test

    use libflit
    use mod_su

    implicit none

    type(su) :: w
    real, allocatable, dimension(:) :: d, t
    integer :: i, j
    real, allocatable, dimension(:, :) :: data

    call fh5_open('./data_and_pick.h5', fh5_fid, mode='r')

    call make_directory('./data_raw')
    call make_directory('./data_processed')

    do i = 1, 49

        call fh5_read_attr(fh5_fid, '/', 'nt', w%nt)
        call fh5_read_attr(fh5_fid, '/', 'dt', w%dt)
        call fh5_read(fh5_fid, '/Shot '//num2str(i)//'/Data', data)
        !        print *, minval(data), maxval(data), shape(data)
        w%nr = size(data, 2)

        call w%init(nt=w%nt, nr=w%nr)

        call w%from_array(data)
        call w%output('./data_raw/shot_' // num2str(i) // '_seismogram_z.su')

        ! Correct time delay for the 13th shot
        if (i == 13) then
            d = zeros(w%nt)
            do j = 1, w%nr
                d(1:12) = 0.0
                d(13:2000) = w%trace(j)%data(1:2000-13 + 1)
                w%trace(j)%data = d
            end do
        end if

        ! Processing data
        call w%resamp(ddt=0.001, nnt=501)
        call w%freqfilt(f=[10.0, 20.0, 120.0, 130.0], a=[0.0, 1.0, 1.0, 0.0])
        call w%output('./data_processed/shot_' // num2str(i) // '_seismogram_z.su')

        call output_array(w%to_array(), './data_processed/shot_' // num2str(i) // '.bin')

        print *, i

    end do

    ! Traveltime picks
    call make_directory('./pick')

    do i = 1, 49
        t = zeros(48)
        call fh5_read(fh5_fid, '/Shot '//num2str(i)//'/Pick', t)
        call output_array(t, './pick/shot_'//num2str(i)//'_traveltime_p.bin')
        open(3, file='./pick/time_'//num2str(i)//'.txt')
        do j = 1, 48
            write(3, *) t((i - 1)*48 + j), j - 1
        end do
        print *, i
    end do

    call fh5_close(fh5_fid)

    ! FATT geometry
    call make_directory('./geometry')

    open(3, file='./geometry/geometry.txt')
    do i = 1, 49

        write(3, *) 'shot_'//num2str(i)//'_geometry.txt'

        open(33, file='./geometry/shot_'//num2str(i)//'_geometry.txt')
        write(33, *) i
        write(33, *)
        write(33, *) 1
        write(33, *) (i - 1)*2.0, 0.0, 0.0, 0.0
        write(33, *)
        write(33, *) 48
        do j = 1, 48
            write(33, *) 1.0 + (j - 1)*2.0, 0.0, 0.0, 1.0
        end do
        close(33)

    end do
    close(3)

    ! FWI geometry
    call make_directory('./geometry_fwi')

    open(3, file='./geometry_fwi/geometry.txt')
    do i = 1, 49

        write(3, *) 'shot_'//num2str(i)//'_geometry.txt'

        open(33, file='./geometry_fwi/shot_'//num2str(i)//'_geometry.txt')
        write(33, *) i
        write(33, *)
        write(33, *) 1
        write(33, *) (i - 1)*2.0, 0.0, 0.0
        write(33, *) 'force', 0.0, 0.0
        write(33, *) 'custom', 60.0, -1e6, 0.0
        write(33, *) './geometry_fwi/wavelet.txt'
        write(33, *) 0, 0
        write(33, *)
        write(33, *) 48
        do j = 1, 48
            write(33, *) 1.0 + (j - 1)*2.0, 0.0, 0.0, 1.0
        end do
        close(33)

    end do
    close(3)

end program
