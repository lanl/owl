
program test

    use libflit

    implicit none

    real, allocatable, dimension(:, :, :) :: w
    real, allocatable, dimension(:, :) :: sxyz
    integer :: nx, ny, nz
    integer :: i, j, k, l

    call make_directory('./model')

    sxyz = zeros(20000, 3)

    nx = 101
    ny = 101
    nz = 101

    l = 1
    do k = 3, nz, 3
        do j = 3, ny, 3
            do i = 3, nx, 3

                if ((i == 3 .or. i == 97) .and. j >= 3 .and. j <= 98 .and. k >= 3 .and. k <= 98) then
                    sxyz(l, :) = [i, j, k]
                    l = l + 1
                end if

                if ((j == 3 .or. j == 97) .and. i >= 3 .and. i <= 98 .and. k >= 3 .and. k <= 98) then
                    sxyz(l, :) = [i, j, k]
                    l = l + 1
                end if

                if ((k == 3 .or. k == 97) .and. j >= 3 .and. j <= 98 .and. i >= 3 .and. i <= 98) then
                    sxyz(l, :) = [i, j, k]
                    l = l + 1
                end if

            end do
        end do
    end do

    l = l - 1
    sxyz = sxyz(1:l, :)

    sxyz = unique(sxyz, cols=[1, 2, 3])

    ! Geometry
    call make_directory('./geometry')
    open(3, file='./geometry/geometry.txt')
    write(3, *) 'shot_1_geometry.txt'
    close(3)

    open(33, file='./geometry/shot_1_geometry.txt')
    write(33, *) 1
    write(33, *)
    write(33, *) 1
    write(33, *) 500.0, 500.0, 500.0

    write(33, '(a, 6es)') 'mt', 1.0, -1.0, 0.5, -1.0, 1.0, -1.0
    write(33, *) 'ricker', 15.0, 1e3, 0.0
    write(33, *) 0, 0
    write(33, *)
    write(33, *) size(sxyz, 1)
    do i = 1, size(sxyz, 1)
        write(33, *) (sxyz(i, 1) - 1)*10.0, (sxyz(i, 2) - 1)*10.0, (sxyz(i, 3) - 1)*10.0, 1.0
    end do
    close(33)

    call output_array([1.0, -1.0, 0.5, -1.0, 1.0, -1.0], './model/mt.bin')
    call output_array([0.5, 0.5, 1.0, 0.0, 0.0, 0.0], './model/mt_init.bin')

    ! Model
    w = zeros(nz, ny, nx) + 2500.0
    call output_array(w, './model/vp.bin')

    w = zeros(nz, ny, nx) + 1500.0
    call output_array(w, './model/vs.bin')

end program
