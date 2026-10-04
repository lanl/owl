!
! © 2025-2026. Triad National Security, LLC. All rights reserved.
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


module elastic_vhtiort_3d_boundary_saving

    use libflit
    use elastic_vhtiort_3d_vars

    implicit none

    integer :: xbwbeg1
    integer :: xbwbeg2
    integer :: xbwend1
    integer :: xbwend2
    integer :: xbwbegy
    integer :: xbwendy
    integer :: xbwbegz
    integer :: xbwendz

    integer :: ybwbeg1
    integer :: ybwbeg2
    integer :: ybwend1
    integer :: ybwend2
    integer :: ybwbegx
    integer :: ybwendx
    integer :: ybwbegz
    integer :: ybwendz

    integer :: zbwbegx
    integer :: zbwendx
    integer :: zbwbegy
    integer :: zbwendy
    integer :: zbwbeg1
    integer :: zbwbeg2
    integer :: zbwend1
    integer :: zbwend2

    integer :: bwiounit
    integer :: bwrecl

contains

    !
    !> Clip the range [lo, hi] to [n1, n2]; an empty range becomes [n1, n1 - 1],
    !> which is inside the bounds of the block arrays
    !
    subroutine clip_range(lo, hi, n1, n2)

        integer, intent(inout) :: lo, hi
        integer, intent(in) :: n1, n2

        lo = max(lo, n1)
        hi = min(hi, n2)
        if (hi < lo) then
            lo = n1
            hi = n1 - 1
        end if

    end subroutine clip_range

    !
    !> Save final step wavefields for elastic media
    !
    subroutine output_final_step_wavefield

        integer :: funit

        open (newunit=funit, file=tidy(dir_working)//'/shot_' &
            //num2str(sgmtr%id) &
            //'_final_step_wavefield.bin.'//num2str(rankid_group), &
            form='unformatted', access='stream', status='replace', action='write')

        write (funit) vx, vy, vz, stressxx, stressyy, stresszz, stressyz, stressxz, stressxy

        close (funit)

    end subroutine output_final_step_wavefield

    !
    !> Input final step wavefields for elastic media
    !
    subroutine input_final_step_wavefield

        integer :: funit

        open (newunit=funit, file=tidy(dir_working)//'/shot_' &
            //num2str(sgmtr%id) &
            //'_final_step_wavefield.bin.'//num2str(rankid_group), &
            form='unformatted', access='stream', status='old', action='read')

        read (funit) vx, vy, vz, stressxx, stressyy, stresszz, stressyz, stressxz, stressxy

        close (funit, status='delete')

        ! The final step was saved after add_source changed the owned points, so its
        ! halos are stale; make the restored state consistent across blocks
        call commute_array_group(vx, fdhalf)
        call commute_array_group(vy, fdhalf)
        call commute_array_group(vz, fdhalf)
        call commute_array_group(stressxx, fdhalf)
        call commute_array_group(stressyy, fdhalf)
        call commute_array_group(stresszz, fdhalf)
        call commute_array_group(stressyz, fdhalf)
        call commute_array_group(stressxz, fdhalf)
        call commute_array_group(stressxy, fdhalf)

    end subroutine input_final_step_wavefield

    !
    !> Prepare boundary saving
    !
    subroutine prepare_boundary_saving

        ! The layers to save along x are x = 2 - fdhalf to 1 and x = nx + 1 to nx + fdhalf,
        ! over y and z padded by fdhalf; likewise along y and z. Each range is clipped
        ! to this block, and an empty range has size zero (see clip_range).
        ! x
        xbwbeg1 = 2 - fdhalf
        xbwbeg2 = 1
        xbwend1 = nx + 1
        xbwend2 = nx + fdhalf
        xbwbegy = 1 - fdhalf
        xbwendy = ny + fdhalf
        xbwbegz = 1 - fdhalf
        xbwendz = nz + fdhalf

        ! y
        ybwbegx = 1 - fdhalf
        ybwendx = nx + fdhalf
        ybwbeg1 = 2 - fdhalf
        ybwbeg2 = 1
        ybwend1 = ny + 1
        ybwend2 = ny + fdhalf
        ybwbegz = 1 - fdhalf
        ybwendz = nz + fdhalf

        ! z
        zbwbegx = 1 - fdhalf
        zbwendx = nx + fdhalf
        zbwbegy = 1 - fdhalf
        zbwendy = ny + fdhalf
        zbwbeg1 = 2 - fdhalf
        zbwbeg2 = 1
        zbwend1 = nz + 1
        zbwend2 = nz + fdhalf

        call clip_range(xbwbeg1, xbwbeg2, nx1, nx2)
        call clip_range(xbwend1, xbwend2, nx1, nx2)
        call clip_range(xbwbegy, xbwendy, ny1, ny2)
        call clip_range(xbwbegz, xbwendz, nz1, nz2)
        call clip_range(ybwbegx, ybwendx, nx1, nx2)
        call clip_range(ybwbeg1, ybwbeg2, ny1, ny2)
        call clip_range(ybwend1, ybwend2, ny1, ny2)
        call clip_range(ybwbegz, ybwendz, nz1, nz2)
        call clip_range(zbwbegx, zbwendx, nx1, nx2)
        call clip_range(zbwbegy, zbwendy, ny1, ny2)
        call clip_range(zbwbeg1, zbwbeg2, nz1, nz2)
        call clip_range(zbwend1, zbwend2, nz1, nz2)

        bwrecl = 0
        ! x
        bwrecl = bwrecl + (xbwbeg2 - xbwbeg1 + 1)*(xbwendy - xbwbegy + 1)*(xbwendz - xbwbegz + 1)
        bwrecl = bwrecl + (xbwend2 - xbwend1 + 1)*(xbwendy - xbwbegy + 1)*(xbwendz - xbwbegz + 1)
        ! y
        bwrecl = bwrecl + (ybwendx - ybwbegx + 1)*(ybwbeg2 - ybwbeg1 + 1)*(ybwendz - ybwbegz + 1)
        bwrecl = bwrecl + (ybwendx - ybwbegx + 1)*(ybwend2 - ybwend1 + 1)*(ybwendz - ybwbegz + 1)
        ! z
        bwrecl = bwrecl + (zbwendx - zbwbegx + 1)*(zbwendy - zbwbegy + 1)*(zbwbeg2 - zbwbeg1 + 1)
        bwrecl = bwrecl + (zbwendx - zbwbegx + 1)*(zbwendy - zbwbegy + 1)*(zbwend2 - zbwend1 + 1)

        ! For elastic wavefield reconstruction, must save all three particle velocity wavefields
        bwrecl = bwrecl*3

    end subroutine prepare_boundary_saving

    !
    !> Save boundary wavefield
    !
    subroutine save_boundary_wavefield(t)

        integer, intent(in) :: t

        if (bwrecl == 0) then
            return
        end if

        write (bwiounit, rec=t) &
            ! left boundary
        vx(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vx(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2)

    end subroutine save_boundary_wavefield

    !
    !> Inject boundary wavefield as boundary condition:
    !>        wavefield=R, not wavefield=wavefield+R
    !
    subroutine inject_boundary_wavefield(t)

        integer, intent(in) :: t

        if (bwrecl == 0) then
            return
        end if

        read (bwiounit, rec=t) &
            ! left boundary
        vx(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vx(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2)

    end subroutine inject_boundary_wavefield

    !
    !> Open boundary-saving file
    !
    subroutine open_boundary_saving

        if (bwrecl == 0) then
            return
        end if

        open (newunit=bwiounit, file=tidy(dir_working)//'/shot_' &
            //num2str(sgmtr%id) &
            //'_boundary_wavefield.bin.'//num2str(rankid_group), &
            form='unformatted', access='direct', recl=4*bwrecl)

    end subroutine open_boundary_saving

    !
    !> Close boundary saving file
    !
    subroutine close_boundary_saving(delete)

        logical, intent(in), optional :: delete

        if (bwrecl == 0) then
            return
        end if

        if (present(delete)) then
            if (delete) then
                close (bwiounit, status='delete')
            else
                close (bwiounit)
            end if
        else
            close (bwiounit)
        end if

    end subroutine

end module
