!
! © 2025. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!


module mod_utility

    use libflit
    use mod_parameters
    use mod_model
    use mod_source_receiver
    use mod_grid

    implicit none

contains

    !
    ! Iteration directory of model
    !
    function dir_iter_model(iter) result(dir)

        integer :: iter
        character(len=:), allocatable :: dir
        character(len=1024) :: tmpdir

        tmpdir = tidy(dir_working)//'/iteration_'//num2str(iter)//'/model'

        allocate (character(len=len_trim(tmpdir)) :: dir)
        dir = tidy(tmpdir)

    end function dir_iter_model

    !
    ! Iteration directory of synthetic data
    !
    function dir_iter_synthetic(iter) result(dir)

        integer :: iter
        character(len=:), allocatable :: dir
        character(len=1024) :: tmpdir

        if (put_synthetic_in_scratch .and. .not. yn_enforce_update) then
            ! If the forward modeling is for computing misfit only
            ! and the enforce update is turned off, i.e., allowing step size searches
            ! Then put the forward modeling data to scratch/synthetic (a temporary directory)
            tmpdir = tidy(dir_scratch)//'/synthetic'
        else
            ! Otherwise, put the data to current iteration synthetic directory
            tmpdir = tidy(dir_working)//'/iteration_'//num2str(iter)//'/synthetic'
        end if

        allocate (character(len=len_trim(tmpdir)) :: dir)
        dir = tidy(tmpdir)

    end function dir_iter_synthetic

    !
    ! Iteration directory of muted synthetic data
    !
    function dir_iter_synthetic_processed(iter) result(dir)

        integer :: iter
        character(len=:), allocatable :: dir
        character(len=1024) :: tmpdir

        if (put_synthetic_in_scratch .and. .not. yn_enforce_update) then
            ! If the forward modeling is for computing misfit only
            ! and the enforce update is turned off, i.e., allowing step size searches
            ! Then put the forward modeling data to scratch/synthetic (a temporary directory)
            tmpdir = tidy(dir_scratch)//'/synthetic_processed'
        else
            ! Otherwise, put the data to current iteration synthetic directory
            tmpdir = tidy(dir_working)//'/iteration_'//num2str(iter)//'/synthetic_processed'
        end if

        allocate (character(len=len_trim(tmpdir)) :: dir)
        dir = tidy(tmpdir)

    end function dir_iter_synthetic_processed

    !
    ! Iteration directory of adjoint source
    !
    function dir_iter_adjoint_source(iter) result(dir)

        integer :: iter
        character(len=:), allocatable :: dir
        character(len=1024) :: tmpdir

        tmpdir = tidy(dir_working)//'/iteration_'//num2str(iter)//'/adjoint_source'

        allocate (character(len=len_trim(tmpdir)) :: dir)
        dir = tidy(tmpdir)

    end function dir_iter_adjoint_source

    !
    ! Iteration directory of encoded recored data
    !
    function dir_iter_record(iter) result(dir)

        integer :: iter
        character(len=:), allocatable :: dir
        character(len=1024) :: tmpdir

        ! If source encoding, then for each iteration, re-encoding is needed
        ! as the encoding matrix changes over iterations
        !        if (yn_source_encoding) then
        tmpdir = tidy(dir_working)//'/iteration_'//num2str(iter)//'/record_encoded'
        !        else
        tmpdir = tidy(dir_record)
        !        end if

        allocate (character(len=len_trim(tmpdir)) :: dir)
        dir = tidy(tmpdir)

    end function dir_iter_record

    !
    ! Iteration directory of base seismograms
    !
    function dir_iter_base(iter) result(dir)

        integer :: iter
        character(len=:), allocatable :: dir
        character(len=1024) :: tmpdir

        !        if (yn_source_encoding) then
        tmpdir = tidy(dir_working)//'/iteration_'//num2str(iter)//'/base_encoded'
        !        else
        tmpdir = tidy(dir_base)
        !        end if

        allocate (character(len=len_trim(tmpdir)) :: dir)
        dir = tidy(tmpdir)

    end function dir_iter_base

    !
    ! Make necessary directories in each iteration
    !
    subroutine make_iter_dir

        character(len=1024) :: dir_iteration

        dir_iteration = tidy(dir_working)//'/iteration_'//num2str(iter)

        if (rankid == 0) then

            call make_directory(tidy(dir_iteration)//'/synthetic')
            call make_directory(tidy(dir_iteration)//'/adjoint_source')
            !            if (yn_source_encoding) then
            !                call make_directory(tidy(dir_iteration)//'/record_encoded')
            !                call make_directory(tidy(dir_iteration)//'/base_encoded')
            !            end if
            call make_directory(tidy(dir_iteration)//'/model')

        end if

        call mpibarrier

    end subroutine make_iter_dir

    !
    ! Print error info
    !
    subroutine print_misfit

        integer :: i, tmpi, tmpr
        character(len=1) :: nl = char(10)
        real :: misfit0
        integer :: m

        if (resume_from_iter == 1 .and. iter == 0) then

            misfit0 = 0.0

            ! Allocate memory for misfit arrays
            !            call alloc_array(shot_misfit, [1, nss, 0, niter_max])
            !            call alloc_array(step_misfit, [1, nss])
            call alloc_array(shot_misfit, [1, ns, 0, niter_max])
            call alloc_array(step_misfit, [1, ns])
            call alloc_array(data_misfit, [0, niter_max + 1])

            ! Create misfit history ascii file if new inversion
            open (33, file=tidy(file_data_misfit), status='replace')
            close (33)

        end if

        if (resume_from_iter > 1 .and. iter == 0) then

            ! Allocate memory for misfit arrays
            !            call alloc_array(shot_misfit, [1, nss, 0, niter_max])
            !            call alloc_array(step_misfit, [1, nss])
            call alloc_array(shot_misfit, [1, ns, 0, niter_max])
            call alloc_array(step_misfit, [1, ns])
            call alloc_array(data_misfit, [0, niter_max + 1])

            ! Read existing misfit history if continued inversion
            open (33, file=tidy(file_data_misfit), status='old')
            do i = 0, resume_from_iter - 1
                read (33, *) tmpi, data_misfit(i), tmpr
            end do
            close (33)

            open (33, file=tidy(file_shot_misfit), status='old', access='stream', form='unformatted')
            do i = 0, resume_from_iter - 1
                read (33) shot_misfit(:, i)
            end do
            close (33)

            misfit0 = data_misfit(0)

        end if

        ! Only in rank 0
        if (rankid == 0 .and. iter >= 1) then

            ! Exit when encountering NaN
            if (.not. (data_misfit(iter) <= float_huge) .or. isnan(data_misfit(iter))) then
                call warn(date_time_compact()//' Error: Misfit is NaN. Exiting. ')
                call mpistop
            end if

            if (iter == 1 .and. misfit0 == 0) then
                misfit0 = data_misfit(0)
                open (33, file=tidy(file_data_misfit), position='append')
                write (33, '(i4,es18.7,es18.7)') 0, misfit0, 1.0
                close (33)
            end if

            ! Misfit information
            call warn(date_time_compact()//' >>>>>>>>>> Data misfit: ' &
                //num2str(data_misfit(iter), '(es18.7)'))
            if (data_misfit(0) == 0) then
                call warn(date_time_compact()//' >>>>>>>>>> Normalized data misfit: ' &
                    //num2str(0.0, '(es18.7)'))
                stop
            else
                call warn(date_time_compact()//' >>>>>>>>>> Normalized data misfit: ' &
                    //num2str(data_misfit(iter)/data_misfit(0), '(es18.7)'))
            end if

            ! Save misfit history to ASCII file
            open (33, file=tidy(file_data_misfit), status='old', form='formatted', access='direct', recl=41)
            write (33, '(i4,es18.7,es18.7,a)', rec=iter + 1) &
                iter, &
                data_misfit(iter), &
                data_misfit(iter)/data_misfit(0), nl
            close (33)

            ! For an inversion that resumes from iter > 1, truncate
            m = count_nonempty_lines(file_data_misfit)
            if (m > iter + 1) then
                open (33, file=tidy(file_data_misfit), status='replace')
                do i = 0, iter
                    write (33, '(i4,es18.7,es18.7,a)') i, data_misfit(i), data_misfit(i)/data_misfit(0)
                end do
                close (33)
            end if

            ! Output shot misfit
            shot_misfit(:, iter) = step_misfit
            call output_array(shot_misfit(:, 0:iter), tidy(file_shot_misfit))

        end if

        step_misfit = 0.0

        call mpibarrier

        if (iter >= 3) then
            if (data_misfit(iter) == data_misfit(iter - 1) .and. data_misfit(iter - 1) == data_misfit(iter - 2)) then

                if (yn_flat_stop) then
                    ! If flat-stop is required, then stop

                    if (rankid == 0) then
                        call warn(' ')
                        call warn(date_time_compact()//' The inversion has three successive iterations with identical data misfits. ')
                        call warn(date_time_compact()//' Inversion has stopped. Exiting. ')
                        call warn(' ')
                    end if
                    call mpibarrier
                    call mpiend

                else

                    ! Otherwise, when the program arrives here, it means the inversion cannot perform effective misfit reduction
                    ! anymore. Most likely, it is in some local minimum and cannot find an effective gradient.
                    ! In this case, in the next iteration, we allow some relaxation for misfits.
                    trigger_jumpout = .true.

                    if (rankid == 0) then
                        call warn(date_time_compact()//' Misfit relaxation is enabled. ')
                        call warn(date_time_compact()//' From next iteration on, misfit may increase. ')
                    end if

                end if

            end if
        end if

    end subroutine print_misfit

    !
    ! Print step number information
    !
    subroutine print_step_info(step_number, step, data_misfit)

        integer, intent(in) :: step_number
        real, intent(in) :: step, data_misfit

        if (rankid == 0) then
            if (.not. (data_misfit <= float_huge)) then
                call warn(date_time_compact()//' Error: Misfit is NaN. Exiting. ')
                call mpistop
            end if
            call warn(date_time_compact()//' Trial step number: '//num2str(step_number))
            call warn(date_time_compact()//' Trial step size: '//num2str(step, '(es18.7)'))
            call warn(date_time_compact()//' Data misfit: '//num2str(data_misfit, '(es18.7)'))
        end if

    end subroutine print_step_info

    !
    ! Divide shots into different groups
    !
    subroutine divide_shots

        if (rankid == 0) then
            call assert(ngroup <= nrank, ' <divide_shots> Error: ngroup must <= nrank. Exiting. ')
            call assert(ngroup <= ns, ' <divide_shots> Error: ngroup must <= ns. Exiting. ')
            call assert(ngroup*rank1_group*rank2_group*rank3_group == nrank, &
                ' <divide_shots> Error: ngroup * rankx * ranky * rankz must = nrank. Exiting. ')
        end if
        call mpibarrier

        call alloc_array(shot_in_group, [0, ngroup - 1, 1, 2])
        call cut(1, ns, ngroup, shot_in_group)

    end subroutine divide_shots

    !
    ! Check stability and dispersion
    !
    subroutine check_dt_f0(dt, dtstable, f0, f0clean)

        real, intent(in) :: dtstable, f0, f0clean
        real, intent(inout) :: dt

        if (dtstable == 0) then
            if (rankid_group == 0) then
                call warn(date_time_compact()//' Error: Stable dt = 0 ')
            end if
            stop
        end if

        if (dt > dtstable) then
            if (rankid_group == 0) then
                call warn(date_time_compact() &
                    //' Warning: dt = '//num2str(dt, '(es)') &
                    //' > stable dt = '//num2str(dtstable, '(es)')//' s')
            end if
            dt = nice(0.95*dtstable, 0.25)
            if (rankid_group == 0) then
                call warn(date_time_compact() &
                    //' Warning: Set dt = '//num2str(dt, '(es)')//' s')
            end if
        end if

        if (dt == 0) then
            dt = nice(0.95*dtstable, 0.25)
            if (rankid_group == 0) then
                call warn(date_time_compact() &
                    //' Warning: dt = 0, set dt = '//num2str(dt, '(es)')//' s')
            end if
        end if

        if (f0 > f0clean) then
            if (rankid_group == 0) then
                call warn(date_time_compact() &
                    //' Warning: f0 = '//num2str(f0, '(es)') &
                    //' > clean f0 = '//num2str(f0clean, '(es)')//' Hz')
            end if
        end if

    end subroutine check_dt_f0


    !
    !> Set adaptive range
    !
    subroutine set_adaptive_range(geom)

        type(source_receiver_geometry), intent(inout) :: geom

        ! then set 1 or 0 based on current offset min and max
        where (geom%recr(:)%aoff < offset_min)
            geom%recr(:)%weight = 0.0
        end where
        where (geom%recr(:)%aoff > offset_max)
            geom%recr(:)%weight = 0.0
        end where

        ! range selection
        if (yn_adpx) then
            shot_xbeg = max(xmin, min( &
                minval(geom%recr(:)%x, mask=geom%recr(:)%weight /= 0), &
                minval(geom%srcr(:)%x)) - adp_extrax)
            shot_xend = min(xmax, max( &
                maxval(geom%recr(:)%x, mask=geom%recr(:)%weight /= 0), &
                maxval(geom%srcr(:)%x)) + adp_extrax)
            shot_xbeg = floor((shot_xbeg - ox)/dx)*dx + ox
            shot_xend = ceiling((shot_xend - ox)/dx)*dx + ox
            where (geom%recr(:)%x < shot_xbeg) geom%recr(:)%weight = 0.0
            where (geom%recr(:)%x > shot_xend) geom%recr(:)%weight = 0.0
        else
            shot_xbeg = ox
            shot_xend = ox + (nx - 1)*dx
        end if

        if (yn_adpy) then
            shot_ybeg = max(ymin, min( &
                minval(geom%recr(:)%y, mask=geom%recr(:)%weight /= 0), &
                minval(geom%srcr(:)%y)) - adp_extray)
            shot_yend = min(ymax, max( &
                maxval(geom%recr(:)%y, mask=geom%recr(:)%weight /= 0), &
                maxval(geom%srcr(:)%y)) + adp_extray)
            shot_ybeg = floor((shot_ybeg - oy)/dy)*dy + oy
            shot_yend = ceiling((shot_yend - oy)/dy)*dy + oy
            where (geom%recr(:)%y < shot_ybeg) geom%recr(:)%weight = 0.0
            where (geom%recr(:)%y > shot_yend) geom%recr(:)%weight = 0.0
        else
            shot_ybeg = oy
            shot_yend = oy + (ny - 1)*dy
        end if

        if (yn_adpz) then
            shot_zbeg = max(zmin, min( &
                minval(geom%recr(:)%z, mask=geom%recr(:)%weight /= 0), &
                minval(geom%srcr(:)%z)) - adp_extraz)
            shot_zend = min(zmax, max( &
                maxval(geom%recr(:)%z, mask=geom%recr(:)%weight /= 0), &
                maxval(geom%srcr(:)%z)) + adp_extraz)
            shot_zbeg = floor((shot_zbeg - oz)/dz)*dz + oz
            shot_zend = ceiling((shot_zend - oz)/dz)*dz + oz
            where (geom%recr(:)%z < shot_zbeg) geom%recr(:)%weight = 0.0
            where (geom%recr(:)%z > shot_zend) geom%recr(:)%weight = 0.0
        else
            shot_zbeg = oz
            shot_zend = oz + (nz - 1)*dz
        end if

        ! Range for each shot
        shot_nxbeg = clip(int((shot_xbeg - ox)/dx) + 1, 1, nx)
        shot_nxend = clip(int((shot_xend - ox)/dx) + 1, 1, nx)
        shot_nybeg = clip(int((shot_ybeg - oy)/dy) + 1, 1, ny)
        shot_nyend = clip(int((shot_yend - oy)/dy) + 1, 1, ny)
        shot_nzbeg = clip(int((shot_zbeg - oz)/dz) + 1, 1, nz)
        shot_nzend = clip(int((shot_zend - oz)/dz) + 1, 1, nz)

        shot_nx = shot_nxend - shot_nxbeg + 1
        shot_ny = shot_nyend - shot_nybeg + 1
        shot_nz = shot_nzend - shot_nzbeg + 1

        if (verbose .and. rankid_group == 0) then
            call warn(date_time_compact()//' Shot '//num2str(set_srcid(ishot))//' dimensions:')
            call warn(date_time_compact()//'      xmin, xmax = ' &
                //num2str(shot_xbeg, '(es)')//', '//num2str(shot_xend, '(es)'))
            if (shot_ny > 1) then
                call warn(date_time_compact()//'      ymin, ymax = ' &
                    //num2str(shot_ybeg, '(es)')//', '//num2str(shot_yend, '(es)'))
            end if
            call warn(date_time_compact()//'      zmin, zmax = ' &
                //num2str(shot_zbeg, '(es)')//', '//num2str(shot_zend, '(es)'))
            if (shot_ny > 1) then
                call warn(date_time_compact()//'      nx, ny, nz = ' &
                    //num2str(shot_nx)//', '//num2str(shot_ny)//', '//num2str(shot_nz))
            else
                call warn(date_time_compact()//'      nx, nz = ' &
                    //num2str(shot_nx)//', '//num2str(shot_nz))
            end if
            call warn(date_time_compact()//'      nxmin, nxmax = ' &
                //num2str(shot_nxbeg)//', '//num2str(shot_nxend))
            if (shot_ny > 1) then
                call warn(date_time_compact()//'      nymin, nymax = ' &
                    //num2str(shot_nybeg)//', '//num2str(shot_nyend))
            end if
            call warn(date_time_compact()//'      nzmin, nzmax = ' &
                //num2str(shot_nzbeg)//', '//num2str(shot_nzend))
        end if

    end subroutine

end module
