
program test

    use libflit
    use mod_su

    type(su) :: a, b
    real, allocatable, dimension(:, :) :: c, d
    integer :: ir

    character(len=256) :: dir1, dir2, label1, label2

    ir = 60

    ! x component

    call getpar_string('dir1', dir1, '', required=.true.)
    call getpar_string('dir2', dir2, '', required=.true.)

    call getpar_string('label1', label1, '', required=.true.)
    call getpar_string('label2', label2, '', required=.true.)

    call a%load(tidy(dir1)//'/shot_1_seismogram_x.su')
    call b%load(tidy(dir2)//'/shot_1_seismogram_x.su')

    c = a%to_array()
    print *, 'specfem', minval(c(:, ir)), maxval(c(:, ir))

    d = -b%to_array()
    print *, 'owl', minval(d(:, ir)), maxval(d(:, ir))

    call output_array(c, './x_'//tidy(label1)//'.bin')
    call output_array(d, './x_'//tidy(label2)//'.bin')


    ! z component
    call a%load(tidy(dir1)//'/shot_1_seismogram_z.su')
    call b%load(tidy(dir2)//'/shot_1_seismogram_z.su')

    c = a%to_array()
    print *, 'specfem', minval(c(:, ir)), maxval(c(:, ir))

    d = b%to_array()
    print *, 'owl', minval(d(:, ir)), maxval(d(:, ir))

    call output_array(c, './z_'//tidy(label1)//'.bin')
    call output_array(d, './z_'//tidy(label2)//'.bin')


end program test
