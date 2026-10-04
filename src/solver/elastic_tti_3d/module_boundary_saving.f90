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


module elastic_tti_3d_boundary_saving

    use libflit
    use elastic_tti_3d_vars

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

        write (funit) &
            vx_hxiyiz, vy_hxiyiz, vz_hxiyiz, &
            vx_ixhyiz, vy_ixhyiz, vz_ixhyiz, &
            vx_ixiyhz, vy_ixiyhz, vz_ixiyhz, &
            vx_hxhyhz, vy_hxhyhz, vz_hxhyhz, &
            stressxx_ixiyiz, stressyy_ixiyiz, stresszz_ixiyiz, &
            stressxy_ixiyiz, stressxz_ixiyiz, stressyz_ixiyiz, &
            stressxx_hxhyiz, stressyy_hxhyiz, stresszz_hxhyiz, &
            stressxy_hxhyiz, stressxz_hxhyiz, stressyz_hxhyiz, &
            stressxx_hxiyhz, stressyy_hxiyhz, stresszz_hxiyhz, &
            stressxy_hxiyhz, stressxz_hxiyhz, stressyz_hxiyhz, &
            stressxx_ixhyhz, stressyy_ixhyhz, stresszz_ixhyhz, &
            stressxy_ixhyhz, stressxz_ixhyhz, stressyz_ixhyhz

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

        read (funit) &
            vx_hxiyiz, vy_hxiyiz, vz_hxiyiz, &
            vx_ixhyiz, vy_ixhyiz, vz_ixhyiz, &
            vx_ixiyhz, vy_ixiyhz, vz_ixiyhz, &
            vx_hxhyhz, vy_hxhyhz, vz_hxhyhz, &
            stressxx_ixiyiz, stressyy_ixiyiz, stresszz_ixiyiz, &
            stressxy_ixiyiz, stressxz_ixiyiz, stressyz_ixiyiz, &
            stressxx_hxhyiz, stressyy_hxhyiz, stresszz_hxhyiz, &
            stressxy_hxhyiz, stressxz_hxhyiz, stressyz_hxhyiz, &
            stressxx_hxiyhz, stressyy_hxiyhz, stresszz_hxiyhz, &
            stressxy_hxiyhz, stressxz_hxiyhz, stressyz_hxiyhz, &
            stressxx_ixhyhz, stressyy_ixhyhz, stresszz_ixhyhz, &
            stressxy_ixhyhz, stressxz_ixhyhz, stressyz_ixhyhz

        close (funit, status='delete')

        ! The final step was saved after add_source changed the owned points, so its
        ! halos are stale; make the restored state consistent across blocks
        call commute_array_group(vx_hxiyiz, fdhalf)
        call commute_array_group(vy_hxiyiz, fdhalf)
        call commute_array_group(vz_hxiyiz, fdhalf)
        call commute_array_group(vx_ixhyiz, fdhalf)
        call commute_array_group(vy_ixhyiz, fdhalf)
        call commute_array_group(vz_ixhyiz, fdhalf)
        call commute_array_group(vx_ixiyhz, fdhalf)
        call commute_array_group(vy_ixiyhz, fdhalf)
        call commute_array_group(vz_ixiyhz, fdhalf)
        call commute_array_group(vx_hxhyhz, fdhalf)
        call commute_array_group(vy_hxhyhz, fdhalf)
        call commute_array_group(vz_hxhyhz, fdhalf)
        call commute_array_group(stressxx_ixiyiz, fdhalf)
        call commute_array_group(stressyy_ixiyiz, fdhalf)
        call commute_array_group(stresszz_ixiyiz, fdhalf)
        call commute_array_group(stressxy_ixiyiz, fdhalf)
        call commute_array_group(stressxz_ixiyiz, fdhalf)
        call commute_array_group(stressyz_ixiyiz, fdhalf)
        call commute_array_group(stressxx_hxhyiz, fdhalf)
        call commute_array_group(stressyy_hxhyiz, fdhalf)
        call commute_array_group(stresszz_hxhyiz, fdhalf)
        call commute_array_group(stressxy_hxhyiz, fdhalf)
        call commute_array_group(stressxz_hxhyiz, fdhalf)
        call commute_array_group(stressyz_hxhyiz, fdhalf)
        call commute_array_group(stressxx_hxiyhz, fdhalf)
        call commute_array_group(stressyy_hxiyhz, fdhalf)
        call commute_array_group(stresszz_hxiyhz, fdhalf)
        call commute_array_group(stressxy_hxiyhz, fdhalf)
        call commute_array_group(stressxz_hxiyhz, fdhalf)
        call commute_array_group(stressyz_hxiyhz, fdhalf)
        call commute_array_group(stressxx_ixhyhz, fdhalf)
        call commute_array_group(stressyy_ixhyhz, fdhalf)
        call commute_array_group(stresszz_ixhyhz, fdhalf)
        call commute_array_group(stressxy_ixhyhz, fdhalf)
        call commute_array_group(stressxz_ixhyhz, fdhalf)
        call commute_array_group(stressyz_ixhyhz, fdhalf)

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
        bwrecl = bwrecl*3*4

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
        vx_hxiyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx_hxiyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vx_hxiyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx_hxiyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy_hxiyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy_hxiyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy_hxiyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy_hxiyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz_hxiyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz_hxiyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz_hxiyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz_hxiyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vx_ixhyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx_ixhyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vx_ixhyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx_ixhyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy_ixhyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy_ixhyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy_ixhyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy_ixhyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz_ixhyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz_ixhyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz_ixhyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz_ixhyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vx_ixiyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx_ixiyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vx_ixiyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx_ixiyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy_ixiyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy_ixiyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy_ixiyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy_ixiyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz_ixiyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz_ixiyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz_ixiyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz_ixiyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vx_hxhyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx_hxhyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vx_hxhyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx_hxhyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy_hxhyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy_hxhyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy_hxhyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy_hxhyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz_hxhyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz_hxhyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz_hxhyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz_hxhyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2)

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
        vx_hxiyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx_hxiyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vx_hxiyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx_hxiyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy_hxiyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy_hxiyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy_hxiyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy_hxiyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz_hxiyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz_hxiyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz_hxiyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz_hxiyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_hxiyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vx_ixhyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx_ixhyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vx_ixhyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx_ixhyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy_ixhyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy_ixhyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy_ixhyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy_ixhyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz_ixhyiz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz_ixhyiz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz_ixhyiz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz_ixhyiz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_ixhyiz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vx_ixiyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx_ixiyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vx_ixiyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx_ixiyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy_ixiyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy_ixiyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy_ixiyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy_ixiyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz_ixiyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz_ixiyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz_ixiyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz_ixiyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_ixiyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vx_hxhyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx_hxhyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vx_hxhyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx_hxhyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy_hxhyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy_hxhyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy_hxhyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy_hxhyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz_hxhyhz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz_hxhyhz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz_hxhyhz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz_hxhyhz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_hxhyhz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2)

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
