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


submodule(elastic_tti_3d) elastic_tti_3d_fwi_gradient

    use libflit
    use elastic_tti_3d_vars
    use elastic_tti_3d_boundary_saving
    use elastic_tti_3d_wavefield
    use mod_anisotropy

    implicit none

#define interior_region nx1_interior:nx2_interior, ny1_interior:ny2_interior, nz1_interior:nz2_interior

#include '../../lib/macro_thomsen_3d.f90'
#include '../../lib/macro_alkhalifah_tsvankin_3d.f90'

contains

    module subroutine compute_adjoint(this)

        class(wave_solver_elastic_tti_3d), intent(inout) :: this

        integer :: l, ir, irx, iry, irz, rgx, rgy, rgz
        type(grid3) :: grd
        integer :: i, j, k, t
        real, allocatable, dimension(:, :, :) :: grad
        real, allocatable, dimension(:, :) :: seis_vx, seis_vy, seis_vz
        real :: amp1, amp2, amp3, amp4
        real :: wmin1, wmin2, wmin3
        real :: wmax1, wmax2, wmax3
        logical :: wnan

        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield
        call alloc_adjoint_wavefield

        yn_energy_precond = this%energy_precond
        kernel_v = this%kernel_v
        kernel_a = this%kernel_a

        ! Adjoint source
        if (yn_compx) then
            call this%seis_vxr%load(tidy(dir_adjoint)//'/shot_'//num2str(sgmtr%id)//'_seismogram_x.su')
            call this%seis_vxr%zero_foreign_rank_traces_group
            call this%seis_vxr%resamp(nnt=nt, ddt=dt)
            call this%seis_vxr%collect_group
            seis_vx = this%seis_vxr%to_array()
        else
            seis_vx = zeros(nt, sgmtr%nr)
        end if

        if (yn_compy) then
            call this%seis_vyr%load(tidy(dir_adjoint)//'/shot_'//num2str(sgmtr%id)//'_seismogram_y.su')
            call this%seis_vyr%zero_foreign_rank_traces_group
            call this%seis_vyr%resamp(nnt=nt, ddt=dt)
            call this%seis_vyr%collect_group
            seis_vy = this%seis_vyr%to_array()
        else
            seis_vy = zeros(nt, sgmtr%nr)
        end if

        if (yn_compz) then
            call this%seis_vzr%load(tidy(dir_adjoint)//'/shot_'//num2str(sgmtr%id)//'_seismogram_z.su')
            call this%seis_vzr%zero_foreign_rank_traces_group
            call this%seis_vzr%resamp(nnt=nt, ddt=dt)
            call this%seis_vzr%collect_group
            seis_vz = this%seis_vzr%to_array()
        else
            seis_vz = zeros(nt, sgmtr%nr)
        end if

        ! Elastic compliances
        call alloc_array(s11, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s12, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s13, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s14, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s15, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s16, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s22, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s23, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s24, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s25, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s26, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s33, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s34, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s35, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s36, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s44, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s45, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s46, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s55, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s56, [1, nx, 1, ny, 1, nz], pad=pml + 1)
        call alloc_array(s66, [1, nx, 1, ny, 1, nz], pad=pml + 1)

        !$omp parallel do private(i, j, k) collapse(3)
        do k = nz1 - 1, nz2 + 1
            do j = ny1 - 1, ny2 + 1
                do i = nx1 - 1, nx2 + 1
                    call cij_to_sij( &
                        c11(i, j, k), c12(i, j, k), c13(i, j, k), &
                        c14(i, j, k), c15(i, j, k), c16(i, j, k), &
                        c22(i, j, k), c23(i, j, k), c24(i, j, k), &
                        c25(i, j, k), c26(i, j, k), c33(i, j, k), &
                        c34(i, j, k), c35(i, j, k), c36(i, j, k), &
                        c44(i, j, k), c45(i, j, k), c46(i, j, k), &
                        c55(i, j, k), c56(i, j, k), c66(i, j, k), &
                        s11(i, j, k), s12(i, j, k), s13(i, j, k), &
                        s14(i, j, k), s15(i, j, k), s16(i, j, k), &
                        s22(i, j, k), s23(i, j, k), s24(i, j, k), &
                        s25(i, j, k), s26(i, j, k), s33(i, j, k), &
                        s34(i, j, k), s35(i, j, k), s36(i, j, k), &
                        s44(i, j, k), s45(i, j, k), s46(i, j, k), &
                        s55(i, j, k), s56(i, j, k), s66(i, j, k))
                end do
            end do
        end do
        !$omp end parallel do

        call alloc_array(grad_c11, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c12, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c13, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c14, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c15, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c16, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c22, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c23, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c24, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c25, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c26, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c33, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c34, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c35, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c36, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c44, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c45, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c46, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c55, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c56, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_c66, [1, nx, 1, ny, 1, nz])
        call alloc_array(grad_rho, [1, nx, 1, ny, 1, nz])

        ! Allocate memory for arrays
        grad_c11 = zeros(nx, ny, nz)
        grad_c12 = zeros(nx, ny, nz)
        grad_c13 = zeros(nx, ny, nz)
        grad_c14 = zeros(nx, ny, nz)
        grad_c15 = zeros(nx, ny, nz)
        grad_c16 = zeros(nx, ny, nz)
        grad_c22 = zeros(nx, ny, nz)
        grad_c23 = zeros(nx, ny, nz)
        grad_c24 = zeros(nx, ny, nz)
        grad_c25 = zeros(nx, ny, nz)
        grad_c26 = zeros(nx, ny, nz)
        grad_c33 = zeros(nx, ny, nz)
        grad_c34 = zeros(nx, ny, nz)
        grad_c35 = zeros(nx, ny, nz)
        grad_c36 = zeros(nx, ny, nz)
        grad_c44 = zeros(nx, ny, nz)
        grad_c45 = zeros(nx, ny, nz)
        grad_c46 = zeros(nx, ny, nz)
        grad_c55 = zeros(nx, ny, nz)
        grad_c56 = zeros(nx, ny, nz)
        grad_c66 = zeros(nx, ny, nz)
        grad_rho = zeros(nx, ny, nz)
        grad_mt = zeros(nc_mt)

        energy_src_v = zeros(nx, ny, nz)
        energy_rec_v = zeros(nx, ny, nz)
        energy_src_a = zeros(nx, ny, nz)
        energy_rec_a = zeros(nx, ny, nz)

        if (kernel_v == '') then
            call alloc_array(strainxx, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainyy, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainzz, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainyz, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainxz, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainxy, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainxxr, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainyyr, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainzzr, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainyzr, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainxzr, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(strainxyr, [nx1, nx2, ny1, ny2, nz1, nz2])
        else
            call alloc_array(strainxx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainyy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainzz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainyz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxxr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainyyr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainzzr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainyzr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxzr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxyr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxx_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainyy_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainzz_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainyz_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxz_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxy_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxxr_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainyyr_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainzzr_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainyzr_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxzr_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(strainxyr_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
        end if

        if (kernel_a == '') then
            call alloc_array(src_vx, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(rec_vx, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(src_vy, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(rec_vy, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(src_vz, [nx1, nx2, ny1, ny2, nz1, nz2])
            call alloc_array(rec_vz, [nx1, nx2, ny1, ny2, nz1, nz2])
        else
            call alloc_array(src_vx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(rec_vx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(src_vy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(rec_vy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(src_vz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(rec_vz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(src_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
            call alloc_array(rec_hilbert, [nx1, nx2, ny1, ny2, nz1, nz2], pad=htlen)
        end if

        call alloc_array(prev_stressxx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_stressyy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_stresszz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_stressyz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_stressxz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_stressxy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_vx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_vy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(prev_vz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(snapvx, [1, nx, 1, ny, 1, nz], pad=pml)
        call alloc_array(snapvy, [1, nx, 1, ny, 1, nz], pad=pml)
        call alloc_array(snapvz, [1, nx, 1, ny, 1, nz], pad=pml)

        ! Prepare boundary saving
        call prepare_boundary_saving
        call open_boundary_saving

        ! Compute gradients by cross-correlation
        l = np
        do t = nt, 1, -1

            if (yn_update_medium .and. t >= sgmtr%srcr(1)%hnt) then

                ! -------------- Wavefield reconstruction ----------------------
                ! Store previous stress wavefields for cross-correlation
                !$omp parallel do private(i, j, k) collapse(3)
                do k = nz1, nz2
                    do j = ny1, ny2
                        do i = nx1, nx2
                            prev_stressxx(i, j, k) = &
                                stressxx_ixiyiz(i, j, k) &
                                + 0.25*sum(stressxx_hxhyiz(i:i + 1, j:j + 1, k)) &
                                + 0.25*sum(stressxx_hxiyhz(i:i + 1, j, k:k + 1)) &
                                + 0.25*sum(stressxx_ixhyhz(i, j:j + 1, k:k + 1))
                            prev_stressyy(i, j, k) = &
                                stressyy_ixiyiz(i, j, k) &
                                + 0.25*sum(stressyy_hxhyiz(i:i + 1, j:j + 1, k)) &
                                + 0.25*sum(stressyy_hxiyhz(i:i + 1, j, k:k + 1)) &
                                + 0.25*sum(stressyy_ixhyhz(i, j:j + 1, k:k + 1))
                            prev_stresszz(i, j, k) = &
                                stresszz_ixiyiz(i, j, k) &
                                + 0.25*sum(stresszz_hxhyiz(i:i + 1, j:j + 1, k)) &
                                + 0.25*sum(stresszz_hxiyhz(i:i + 1, j, k:k + 1)) &
                                + 0.25*sum(stresszz_ixhyhz(i, j:j + 1, k:k + 1))
                            prev_stressyz(i, j, k) = &
                                stressyz_ixiyiz(i, j, k) &
                                + 0.25*sum(stressyz_hxhyiz(i:i + 1, j:j + 1, k)) &
                                + 0.25*sum(stressyz_hxiyhz(i:i + 1, j, k:k + 1)) &
                                + 0.25*sum(stressyz_ixhyhz(i, j:j + 1, k:k + 1))
                            prev_stressxz(i, j, k) = &
                                stressxz_ixiyiz(i, j, k) &
                                + 0.25*sum(stressxz_hxhyiz(i:i + 1, j:j + 1, k)) &
                                + 0.25*sum(stressxz_hxiyhz(i:i + 1, j, k:k + 1)) &
                                + 0.25*sum(stressxz_ixhyhz(i, j:j + 1, k:k + 1))
                            prev_stressxy(i, j, k) = &
                                stressxy_ixiyiz(i, j, k) &
                                + 0.25*sum(stressxy_hxhyiz(i:i + 1, j:j + 1, k)) &
                                + 0.25*sum(stressxy_hxiyhz(i:i + 1, j, k:k + 1)) &
                                + 0.25*sum(stressxy_ixhyhz(i, j:j + 1, k:k + 1))
                            prev_vx(i, j, k) = &
                                0.5*sum(vx_hxiyiz(i:i + 1, j, k)) &
                                + 0.5*sum(vx_ixhyiz(i, j:j + 1, k)) &
                                + 0.5*sum(vx_ixiyhz(i, j, k:k + 1)) &
                                + 0.125*sum(vx_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                            prev_vy(i, j, k) = &
                                0.5*sum(vy_hxiyiz(i:i + 1, j, k)) &
                                + 0.5*sum(vy_ixhyiz(i, j:j + 1, k)) &
                                + 0.5*sum(vy_ixiyhz(i, j, k:k + 1)) &
                                + 0.125*sum(vy_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                            prev_vz(i, j, k) = &
                                0.5*sum(vz_hxiyiz(i:i + 1, j, k)) &
                                + 0.5*sum(vz_ixhyiz(i, j:j + 1, k)) &
                                + 0.5*sum(vz_ixiyhz(i, j, k:k + 1)) &
                                + 0.125*sum(vz_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                        end do
                    end do
                end do
                !$omp end parallel do

                ! -------------- Forward wavefield reconstruction -----------------------
                if (yn_free_surface) then
                    call update_wavefield_free_surface(-dt, &
                        stressxx_ixiyiz, stressyy_ixiyiz, stresszz_ixiyiz, &
                        stressxy_ixiyiz, stressxz_ixiyiz, stressyz_ixiyiz, &
                        memory_pdxvx_ixiyiz, &
                        memory_pdxvy_ixiyiz, &
                        memory_pdxvz_ixiyiz, &
                        memory_pdyvx_ixiyiz, &
                        memory_pdyvy_ixiyiz, &
                        memory_pdyvz_ixiyiz, &
                        memory_pdzvx_ixiyiz, &
                        memory_pdzvy_ixiyiz, &
                        memory_pdzvz_ixiyiz, &
                        stressxx_hxhyiz, stressyy_hxhyiz, stresszz_hxhyiz, &
                        stressxy_hxhyiz, stressxz_hxhyiz, stressyz_hxhyiz, &
                        memory_pdxvx_hxhyiz, &
                        memory_pdxvy_hxhyiz, &
                        memory_pdxvz_hxhyiz, &
                        memory_pdyvx_hxhyiz, &
                        memory_pdyvy_hxhyiz, &
                        memory_pdyvz_hxhyiz, &
                        memory_pdzvx_hxhyiz, &
                        memory_pdzvy_hxhyiz, &
                        memory_pdzvz_hxhyiz, &
                        stressxx_hxiyhz, stressyy_hxiyhz, stresszz_hxiyhz, &
                        stressxy_hxiyhz, stressxz_hxiyhz, stressyz_hxiyhz, &
                        memory_pdxvx_hxiyhz, &
                        memory_pdxvy_hxiyhz, &
                        memory_pdxvz_hxiyhz, &
                        memory_pdyvx_hxiyhz, &
                        memory_pdyvy_hxiyhz, &
                        memory_pdyvz_hxiyhz, &
                        memory_pdzvx_hxiyhz, &
                        memory_pdzvy_hxiyhz, &
                        memory_pdzvz_hxiyhz, &
                        stressxx_ixhyhz, stressyy_ixhyhz, stresszz_ixhyhz, &
                        stressxy_ixhyhz, stressxz_ixhyhz, stressyz_ixhyhz, &
                        memory_pdxvx_ixhyhz, &
                        memory_pdxvy_ixhyhz, &
                        memory_pdxvz_ixhyhz, &
                        memory_pdyvx_ixhyhz, &
                        memory_pdyvy_ixhyhz, &
                        memory_pdyvz_ixhyhz, &
                        memory_pdzvx_ixhyhz, &
                        memory_pdzvy_ixhyhz, &
                        memory_pdzvz_ixhyhz, &
                        vx_hxiyiz, vy_hxiyiz, vz_hxiyiz, &
                        memory_pdxxx_hxiyiz, &
                        memory_pdxxy_hxiyiz, &
                        memory_pdxxz_hxiyiz, &
                        memory_pdyxy_hxiyiz, &
                        memory_pdyyy_hxiyiz, &
                        memory_pdyyz_hxiyiz, &
                        memory_pdzxx_hxiyiz, memory_pdzxy_hxiyiz, memory_pdzxz_hxiyiz, &
                        memory_pdzyy_hxiyiz, memory_pdzyz_hxiyiz, &
                        memory_pdzzz_hxiyiz, &
                        vx_ixhyiz, vy_ixhyiz, vz_ixhyiz, &
                        memory_pdxxx_ixhyiz, &
                        memory_pdxxy_ixhyiz, &
                        memory_pdxxz_ixhyiz, &
                        memory_pdyxy_ixhyiz, &
                        memory_pdyyy_ixhyiz, &
                        memory_pdyyz_ixhyiz, &
                        memory_pdzxx_ixhyiz, memory_pdzxy_ixhyiz, memory_pdzxz_ixhyiz, &
                        memory_pdzyy_ixhyiz, memory_pdzyz_ixhyiz, &
                        memory_pdzzz_ixhyiz, &
                        vx_ixiyhz, vy_ixiyhz, vz_ixiyhz, &
                        memory_pdxxx_ixiyhz, &
                        memory_pdxxy_ixiyhz, &
                        memory_pdxxz_ixiyhz, &
                        memory_pdyxy_ixiyhz, &
                        memory_pdyyy_ixiyhz, &
                        memory_pdyyz_ixiyhz, &
                        memory_pdzxx_ixiyhz, memory_pdzxy_ixiyhz, memory_pdzxz_ixiyhz, &
                        memory_pdzyy_ixiyhz, memory_pdzyz_ixiyhz, &
                        memory_pdzzz_ixiyhz, &
                        vx_hxhyhz, vy_hxhyhz, vz_hxhyhz, &
                        memory_pdxxx_hxhyhz, &
                        memory_pdxxy_hxhyhz, &
                        memory_pdxxz_hxhyhz, &
                        memory_pdyxy_hxhyhz, &
                        memory_pdyyy_hxhyhz, &
                        memory_pdyyz_hxhyhz, &
                        memory_pdzxx_hxhyhz, memory_pdzxy_hxhyhz, memory_pdzxz_hxhyhz, &
                        memory_pdzyy_hxhyhz, memory_pdzyz_hxhyhz, &
                        memory_pdzzz_hxhyhz)
                else
                    call update_wavefield(-dt, &
                        stressxx_ixiyiz, stressyy_ixiyiz, stresszz_ixiyiz, &
                        stressxy_ixiyiz, stressxz_ixiyiz, stressyz_ixiyiz, &
                        memory_pdxvx_ixiyiz, &
                        memory_pdxvy_ixiyiz, &
                        memory_pdxvz_ixiyiz, &
                        memory_pdyvx_ixiyiz, &
                        memory_pdyvy_ixiyiz, &
                        memory_pdyvz_ixiyiz, &
                        memory_pdzvx_ixiyiz, &
                        memory_pdzvy_ixiyiz, &
                        memory_pdzvz_ixiyiz, &
                        stressxx_hxhyiz, stressyy_hxhyiz, stresszz_hxhyiz, &
                        stressxy_hxhyiz, stressxz_hxhyiz, stressyz_hxhyiz, &
                        memory_pdxvx_hxhyiz, &
                        memory_pdxvy_hxhyiz, &
                        memory_pdxvz_hxhyiz, &
                        memory_pdyvx_hxhyiz, &
                        memory_pdyvy_hxhyiz, &
                        memory_pdyvz_hxhyiz, &
                        memory_pdzvx_hxhyiz, &
                        memory_pdzvy_hxhyiz, &
                        memory_pdzvz_hxhyiz, &
                        stressxx_hxiyhz, stressyy_hxiyhz, stresszz_hxiyhz, &
                        stressxy_hxiyhz, stressxz_hxiyhz, stressyz_hxiyhz, &
                        memory_pdxvx_hxiyhz, &
                        memory_pdxvy_hxiyhz, &
                        memory_pdxvz_hxiyhz, &
                        memory_pdyvx_hxiyhz, &
                        memory_pdyvy_hxiyhz, &
                        memory_pdyvz_hxiyhz, &
                        memory_pdzvx_hxiyhz, &
                        memory_pdzvy_hxiyhz, &
                        memory_pdzvz_hxiyhz, &
                        stressxx_ixhyhz, stressyy_ixhyhz, stresszz_ixhyhz, &
                        stressxy_ixhyhz, stressxz_ixhyhz, stressyz_ixhyhz, &
                        memory_pdxvx_ixhyhz, &
                        memory_pdxvy_ixhyhz, &
                        memory_pdxvz_ixhyhz, &
                        memory_pdyvx_ixhyhz, &
                        memory_pdyvy_ixhyhz, &
                        memory_pdyvz_ixhyhz, &
                        memory_pdzvx_ixhyhz, &
                        memory_pdzvy_ixhyhz, &
                        memory_pdzvz_ixhyhz, &
                        vx_hxiyiz, vy_hxiyiz, vz_hxiyiz, &
                        memory_pdxxx_hxiyiz, &
                        memory_pdxxy_hxiyiz, &
                        memory_pdxxz_hxiyiz, &
                        memory_pdyxy_hxiyiz, &
                        memory_pdyyy_hxiyiz, &
                        memory_pdyyz_hxiyiz, &
                        memory_pdzxz_hxiyiz, &
                        memory_pdzyz_hxiyiz, &
                        memory_pdzzz_hxiyiz, &
                        vx_ixhyiz, vy_ixhyiz, vz_ixhyiz, &
                        memory_pdxxx_ixhyiz, &
                        memory_pdxxy_ixhyiz, &
                        memory_pdxxz_ixhyiz, &
                        memory_pdyxy_ixhyiz, &
                        memory_pdyyy_ixhyiz, &
                        memory_pdyyz_ixhyiz, &
                        memory_pdzxz_ixhyiz, &
                        memory_pdzyz_ixhyiz, &
                        memory_pdzzz_ixhyiz, &
                        vx_ixiyhz, vy_ixiyhz, vz_ixiyhz, &
                        memory_pdxxx_ixiyhz, &
                        memory_pdxxy_ixiyhz, &
                        memory_pdxxz_ixiyhz, &
                        memory_pdyxy_ixiyhz, &
                        memory_pdyyy_ixiyhz, &
                        memory_pdyyz_ixiyhz, &
                        memory_pdzxz_ixiyhz, &
                        memory_pdzyz_ixiyhz, &
                        memory_pdzzz_ixiyhz, &
                        vx_hxhyhz, vy_hxhyhz, vz_hxhyhz, &
                        memory_pdxxx_hxhyhz, &
                        memory_pdxxy_hxhyhz, &
                        memory_pdxxz_hxhyhz, &
                        memory_pdyxy_hxhyhz, &
                        memory_pdyyy_hxhyhz, &
                        memory_pdyyz_hxhyhz, &
                        memory_pdzxz_hxhyhz, &
                        memory_pdzyz_hxhyhz, &
                        memory_pdzzz_hxhyhz)
                end if

                ! Read final step wavefield
                if (t == nt) then
                    call input_final_step_wavefield
                end if

                ! Read boundary wavefield
                call inject_boundary_wavefield(t)

                ! Record wavefield snapshot if necessary
                if (np /= 0 .and. l <= np) then
                    if (t - 1 == nint(snaps(l)/dt)) then

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

                        if (yn_free_surface) then

                            call alloc_array(snapvx, [1, nx, 1, ny, 1, nz], pad=pml)
                            call alloc_array(snapvy, [1, nx, 1, ny, 1, nz], pad=pml)
                            call alloc_array(snapvz, [1, nx, 1, ny, 1, nz], pad=pml)

                            !$omp parallel do private(i, j, k) collapse(3)
                            do k = nz1, nz2
                                do j = ny1, ny2
                                    do i = nx1, nx2
                                        snapvx(i, j, k) = &
                                            0.5*sum(vx_hxiyiz(i:i + 1, j, k)) &
                                            + 0.5*sum(vx_ixhyiz(i, j:j + 1, k)) &
                                            + 0.5*sum(vx_ixiyhz(i, j, k:k + 1)) &
                                            + 0.125*sum(vx_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                        snapvy(i, j, k) = &
                                            0.5*sum(vy_hxiyiz(i:i + 1, j, k)) &
                                            + 0.5*sum(vy_ixhyiz(i, j:j + 1, k)) &
                                            + 0.5*sum(vy_ixiyhz(i, j, k:k + 1)) &
                                            + 0.125*sum(vy_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                        snapvz(i, j, k) = &
                                            0.5*sum(vz_hxiyiz(i:i + 1, j, k)) &
                                            + 0.5*sum(vz_ixhyiz(i, j:j + 1, k)) &
                                            + 0.5*sum(vz_ixiyhz(i, j, k:k + 1)) &
                                            + 0.125*sum(vz_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                    end do
                                end do
                            end do
                            !$omp end parallel do

                            call reduce_array_group(snapvx)
                            call reduce_array_group(snapvy)
                            call reduce_array_group(snapvz)

                            snapvx = 0.25*snapvx
                            snapvy = 0.25*snapvy
                            snapvz = 0.25*snapvz

                            if (rankid_group == 0) then

                                open(3, file=tidy(dir_snapshot)//'/shot_' &
                                    //num2str(sgmtr%id) &
                                    //'_reconstructed_wavefield_' &
                                    //num2str(l)//'.txt')
                                do i = -pml + 1, nx + pml
                                    do j = -pml + 1, ny + pml
                                        do k = 1, nz + pml
                                            write(3, '(6es)') (i - 1)*dx + ox, (j - 1)*dy + oy, zz_i(i, j, k), &
                                                snapvx(i, j, k), snapvy(i, j, k), snapvz(i, j, k)
                                        end do
                                    end do
                                end do
                                close(3)

                                call map_irregular_to_regular(snapvx, this, [1, this%nx, 1, this%ny, 1, this%nz])
                                call map_irregular_to_regular(snapvy, this, [1, this%nx, 1, this%ny, 1, this%nz])
                                call map_irregular_to_regular(snapvz, this, [1, this%nx, 1, this%ny, 1, this%nz])

                                call output_array(snapvx, tidy(dir_snapshot)//'/shot_' &
                                    //num2str(sgmtr%id) &
                                    //'_reconstructed_wavefield_x_' &
                                    //num2str(l)//'.bin', store=321)
                                call output_array(snapvy, tidy(dir_snapshot)//'/shot_' &
                                    //num2str(sgmtr%id) &
                                    //'_reconstructed_wavefield_y_' &
                                    //num2str(l)//'.bin', store=321)
                                call output_array(snapvz, tidy(dir_snapshot)//'/shot_' &
                                    //num2str(sgmtr%id) &
                                    //'_reconstructed_wavefield_z_' &
                                    //num2str(l)//'.bin', store=321)

                            end if

                        else

                            snapvx = 0.0
                            snapvy = 0.0
                            snapvz = 0.0

                            !$omp parallel do private(i, j, k) collapse(3)
                            do k = nz1, nz2
                                do j = ny1, ny2
                                    do i = nx1, nx2
                                        snapvx(i, j, k) = &
                                            0.5*sum(vx_hxiyiz(i:i + 1, j, k)) &
                                            + 0.5*sum(vx_ixhyiz(i, j:j + 1, k)) &
                                            + 0.5*sum(vx_ixiyhz(i, j, k:k + 1)) &
                                            + 0.125*sum(vx_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                        snapvy(i, j, k) = &
                                            0.5*sum(vy_hxiyiz(i:i + 1, j, k)) &
                                            + 0.5*sum(vy_ixhyiz(i, j:j + 1, k)) &
                                            + 0.5*sum(vy_ixiyhz(i, j, k:k + 1)) &
                                            + 0.125*sum(vy_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                        snapvz(i, j, k) = &
                                            0.5*sum(vz_hxiyiz(i:i + 1, j, k)) &
                                            + 0.5*sum(vz_ixhyiz(i, j:j + 1, k)) &
                                            + 0.5*sum(vz_ixiyhz(i, j, k:k + 1)) &
                                            + 0.125*sum(vz_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                                    end do
                                end do
                            end do
                            !$omp end parallel do

                            call reduce_array_group(snapvx)
                            call reduce_array_group(snapvy)
                            call reduce_array_group(snapvz)

                            snapvx = 0.25*snapvx
                            snapvy = 0.25*snapvy
                            snapvz = 0.25*snapvz

                            ! Output
                            if (rankid_group == 0) then
                                call output_array(snapvx(1:nx, 1:ny, 1:nz), tidy(dir_snapshot)//'/shot_' &
                                    //num2str(sgmtr%id) &
                                    //'_reconstructed_wavefield_x_' &
                                    //num2str(l)//'.bin', store=321)
                                call output_array(snapvy(1:nx, 1:ny, 1:nz), tidy(dir_snapshot)//'/shot_' &
                                    //num2str(sgmtr%id) &
                                    //'_reconstructed_wavefield_y_' &
                                    //num2str(l)//'.bin', store=321)
                                call output_array(snapvz(1:nx, 1:ny, 1:nz), tidy(dir_snapshot)//'/shot_' &
                                    //num2str(sgmtr%id) &
                                    //'_reconstructed_wavefield_z_' &
                                    //num2str(l)//'.bin', store=321)
                            end if

                        end if

                        l = l - 1
                    end if
                end if

            end if

            ! -------------- Adjoint wavefield reverse-time propagation ----------------------
            if (yn_free_surface) then
                call update_wavefield_free_surface(-dt, &
                    stressxxr_ixiyiz, stressyyr_ixiyiz, stresszzr_ixiyiz, &
                    stressxyr_ixiyiz, stressxzr_ixiyiz, stressyzr_ixiyiz, &
                    memory_pdxvxr_ixiyiz, &
                    memory_pdxvyr_ixiyiz, &
                    memory_pdxvzr_ixiyiz, &
                    memory_pdyvxr_ixiyiz, &
                    memory_pdyvyr_ixiyiz, &
                    memory_pdyvzr_ixiyiz, &
                    memory_pdzvxr_ixiyiz, &
                    memory_pdzvyr_ixiyiz, &
                    memory_pdzvzr_ixiyiz, &
                    stressxxr_hxhyiz, stressyyr_hxhyiz, stresszzr_hxhyiz, &
                    stressxyr_hxhyiz, stressxzr_hxhyiz, stressyzr_hxhyiz, &
                    memory_pdxvxr_hxhyiz, &
                    memory_pdxvyr_hxhyiz, &
                    memory_pdxvzr_hxhyiz, &
                    memory_pdyvxr_hxhyiz, &
                    memory_pdyvyr_hxhyiz, &
                    memory_pdyvzr_hxhyiz, &
                    memory_pdzvxr_hxhyiz, &
                    memory_pdzvyr_hxhyiz, &
                    memory_pdzvzr_hxhyiz, &
                    stressxxr_hxiyhz, stressyyr_hxiyhz, stresszzr_hxiyhz, &
                    stressxyr_hxiyhz, stressxzr_hxiyhz, stressyzr_hxiyhz, &
                    memory_pdxvxr_hxiyhz, &
                    memory_pdxvyr_hxiyhz, &
                    memory_pdxvzr_hxiyhz, &
                    memory_pdyvxr_hxiyhz, &
                    memory_pdyvyr_hxiyhz, &
                    memory_pdyvzr_hxiyhz, &
                    memory_pdzvxr_hxiyhz, &
                    memory_pdzvyr_hxiyhz, &
                    memory_pdzvzr_hxiyhz, &
                    stressxxr_ixhyhz, stressyyr_ixhyhz, stresszzr_ixhyhz, &
                    stressxyr_ixhyhz, stressxzr_ixhyhz, stressyzr_ixhyhz, &
                    memory_pdxvxr_ixhyhz, &
                    memory_pdxvyr_ixhyhz, &
                    memory_pdxvzr_ixhyhz, &
                    memory_pdyvxr_ixhyhz, &
                    memory_pdyvyr_ixhyhz, &
                    memory_pdyvzr_ixhyhz, &
                    memory_pdzvxr_ixhyhz, &
                    memory_pdzvyr_ixhyhz, &
                    memory_pdzvzr_ixhyhz, &
                    vxr_hxiyiz, vyr_hxiyiz, vzr_hxiyiz, &
                    memory_pdxxxr_hxiyiz, &
                    memory_pdxxyr_hxiyiz, &
                    memory_pdxxzr_hxiyiz, &
                    memory_pdyxyr_hxiyiz, &
                    memory_pdyyyr_hxiyiz, &
                    memory_pdyyzr_hxiyiz, &
                    memory_pdzxxr_hxiyiz, memory_pdzxyr_hxiyiz, memory_pdzxzr_hxiyiz, &
                    memory_pdzyyr_hxiyiz, memory_pdzyzr_hxiyiz, &
                    memory_pdzzzr_hxiyiz, &
                    vxr_ixhyiz, vyr_ixhyiz, vzr_ixhyiz, &
                    memory_pdxxxr_ixhyiz, &
                    memory_pdxxyr_ixhyiz, &
                    memory_pdxxzr_ixhyiz, &
                    memory_pdyxyr_ixhyiz, &
                    memory_pdyyyr_ixhyiz, &
                    memory_pdyyzr_ixhyiz, &
                    memory_pdzxxr_ixhyiz, memory_pdzxyr_ixhyiz, memory_pdzxzr_ixhyiz, &
                    memory_pdzyyr_ixhyiz, memory_pdzyzr_ixhyiz, &
                    memory_pdzzzr_ixhyiz, &
                    vxr_ixiyhz, vyr_ixiyhz, vzr_ixiyhz, &
                    memory_pdxxxr_ixiyhz, &
                    memory_pdxxyr_ixiyhz, &
                    memory_pdxxzr_ixiyhz, &
                    memory_pdyxyr_ixiyhz, &
                    memory_pdyyyr_ixiyhz, &
                    memory_pdyyzr_ixiyhz, &
                    memory_pdzxxr_ixiyhz, memory_pdzxyr_ixiyhz, memory_pdzxzr_ixiyhz, &
                    memory_pdzyyr_ixiyhz, memory_pdzyzr_ixiyhz, &
                    memory_pdzzzr_ixiyhz, &
                    vxr_hxhyhz, vyr_hxhyhz, vzr_hxhyhz, &
                    memory_pdxxxr_hxhyhz, &
                    memory_pdxxyr_hxhyhz, &
                    memory_pdxxzr_hxhyhz, &
                    memory_pdyxyr_hxhyhz, &
                    memory_pdyyyr_hxhyhz, &
                    memory_pdyyzr_hxhyhz, &
                    memory_pdzxxr_hxhyhz, memory_pdzxyr_hxhyhz, memory_pdzxzr_hxhyhz, &
                    memory_pdzyyr_hxhyhz, memory_pdzyzr_hxhyhz, &
                    memory_pdzzzr_hxhyhz)
            else
                call update_wavefield(-dt, &
                    stressxxr_ixiyiz, stressyyr_ixiyiz, stresszzr_ixiyiz, &
                    stressxyr_ixiyiz, stressxzr_ixiyiz, stressyzr_ixiyiz, &
                    memory_pdxvxr_ixiyiz, &
                    memory_pdxvyr_ixiyiz, &
                    memory_pdxvzr_ixiyiz, &
                    memory_pdyvxr_ixiyiz, &
                    memory_pdyvyr_ixiyiz, &
                    memory_pdyvzr_ixiyiz, &
                    memory_pdzvxr_ixiyiz, &
                    memory_pdzvyr_ixiyiz, &
                    memory_pdzvzr_ixiyiz, &
                    stressxxr_hxhyiz, stressyyr_hxhyiz, stresszzr_hxhyiz, &
                    stressxyr_hxhyiz, stressxzr_hxhyiz, stressyzr_hxhyiz, &
                    memory_pdxvxr_hxhyiz, &
                    memory_pdxvyr_hxhyiz, &
                    memory_pdxvzr_hxhyiz, &
                    memory_pdyvxr_hxhyiz, &
                    memory_pdyvyr_hxhyiz, &
                    memory_pdyvzr_hxhyiz, &
                    memory_pdzvxr_hxhyiz, &
                    memory_pdzvyr_hxhyiz, &
                    memory_pdzvzr_hxhyiz, &
                    stressxxr_hxiyhz, stressyyr_hxiyhz, stresszzr_hxiyhz, &
                    stressxyr_hxiyhz, stressxzr_hxiyhz, stressyzr_hxiyhz, &
                    memory_pdxvxr_hxiyhz, &
                    memory_pdxvyr_hxiyhz, &
                    memory_pdxvzr_hxiyhz, &
                    memory_pdyvxr_hxiyhz, &
                    memory_pdyvyr_hxiyhz, &
                    memory_pdyvzr_hxiyhz, &
                    memory_pdzvxr_hxiyhz, &
                    memory_pdzvyr_hxiyhz, &
                    memory_pdzvzr_hxiyhz, &
                    stressxxr_ixhyhz, stressyyr_ixhyhz, stresszzr_ixhyhz, &
                    stressxyr_ixhyhz, stressxzr_ixhyhz, stressyzr_ixhyhz, &
                    memory_pdxvxr_ixhyhz, &
                    memory_pdxvyr_ixhyhz, &
                    memory_pdxvzr_ixhyhz, &
                    memory_pdyvxr_ixhyhz, &
                    memory_pdyvyr_ixhyhz, &
                    memory_pdyvzr_ixhyhz, &
                    memory_pdzvxr_ixhyhz, &
                    memory_pdzvyr_ixhyhz, &
                    memory_pdzvzr_ixhyhz, &
                    vxr_hxiyiz, vyr_hxiyiz, vzr_hxiyiz, &
                    memory_pdxxxr_hxiyiz, &
                    memory_pdxxyr_hxiyiz, &
                    memory_pdxxzr_hxiyiz, &
                    memory_pdyxyr_hxiyiz, &
                    memory_pdyyyr_hxiyiz, &
                    memory_pdyyzr_hxiyiz, &
                    memory_pdzxzr_hxiyiz, &
                    memory_pdzyzr_hxiyiz, &
                    memory_pdzzzr_hxiyiz, &
                    vxr_ixhyiz, vyr_ixhyiz, vzr_ixhyiz, &
                    memory_pdxxxr_ixhyiz, &
                    memory_pdxxyr_ixhyiz, &
                    memory_pdxxzr_ixhyiz, &
                    memory_pdyxyr_ixhyiz, &
                    memory_pdyyyr_ixhyiz, &
                    memory_pdyyzr_ixhyiz, &
                    memory_pdzxzr_ixhyiz, &
                    memory_pdzyzr_ixhyiz, &
                    memory_pdzzzr_ixhyiz, &
                    vxr_ixiyhz, vyr_ixiyhz, vzr_ixiyhz, &
                    memory_pdxxxr_ixiyhz, &
                    memory_pdxxyr_ixiyhz, &
                    memory_pdxxzr_ixiyhz, &
                    memory_pdyxyr_ixiyhz, &
                    memory_pdyyyr_ixiyhz, &
                    memory_pdyyzr_ixiyhz, &
                    memory_pdzxzr_ixiyhz, &
                    memory_pdzyzr_ixiyhz, &
                    memory_pdzzzr_ixiyhz, &
                    vxr_hxhyhz, vyr_hxhyhz, vzr_hxhyhz, &
                    memory_pdxxxr_hxhyhz, &
                    memory_pdxxyr_hxhyhz, &
                    memory_pdxxzr_hxhyhz, &
                    memory_pdyxyr_hxhyhz, &
                    memory_pdyyyr_hxhyhz, &
                    memory_pdyyzr_hxhyhz, &
                    memory_pdzxzr_hxhyhz, &
                    memory_pdzyzr_hxhyhz, &
                    memory_pdzzzr_hxhyhz)
            end if

            ! Read and add adjoint source
            !$omp parallel private(ir, irx, iry, irz, rgx, rgy, rgz, amp1, amp2, amp3, amp4)
            do ir = 1, sgmtr%nr
                if (sgmtr%recr(ir)%weight /= 0) then

                    rgx = sgmtr%recr(ir)%hx
                    rgy = sgmtr%recr(ir)%gy
                    rgz = sgmtr%recr(ir)%gz
                    !$omp do collapse(3)
                    do irz = -nkw, nkw
                        do iry = -nkw, nkw
                            do irx = -nkw, nkw

                                if (is_in_block(rgx + irx, rgy + iry, rgz + irz) &
                                        .and. ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then

                                    amp1 = sgmtr%recr(ir)%interp_hx(irx) &
                                        *sgmtr%recr(ir)%interp_iy(iry) &
                                        *sgmtr%recr(ir)%interp_iz(irz) &
                                        *sgmtr%recr(ir)%weight*0.25

                                    vxr_hxiyiz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vxr_hxiyiz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vx(t, ir)*amp1
                                    vyr_hxiyiz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vyr_hxiyiz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vy(t, ir)*amp1
                                    vzr_hxiyiz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vzr_hxiyiz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vz(t, ir)*amp1

                                end if

                            end do
                        end do
                    end do
                    !$omp end do

                    rgx = sgmtr%recr(ir)%gx
                    rgy = sgmtr%recr(ir)%hy
                    rgz = sgmtr%recr(ir)%gz
                    !$omp do collapse(3)
                    do irz = -nkw, nkw
                        do iry = -nkw, nkw
                            do irx = -nkw, nkw

                                if (is_in_block(rgx + irx, rgy + iry, rgz + irz) &
                                        .and. ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then

                                    amp2 = sgmtr%recr(ir)%interp_ix(irx) &
                                        *sgmtr%recr(ir)%interp_hy(iry) &
                                        *sgmtr%recr(ir)%interp_iz(irz) &
                                        *sgmtr%recr(ir)%weight*0.25

                                    vxr_ixhyiz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vxr_ixhyiz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vx(t, ir)*amp2
                                    vyr_ixhyiz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vyr_ixhyiz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vy(t, ir)*amp2
                                    vzr_ixhyiz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vzr_ixhyiz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vz(t, ir)*amp2

                                end if

                            end do
                        end do
                    end do
                    !$omp end do

                    rgx = sgmtr%recr(ir)%gx
                    rgy = sgmtr%recr(ir)%gy
                    rgz = sgmtr%recr(ir)%hz
                    !$omp do collapse(3)
                    do irz = -nkw, nkw
                        do iry = -nkw, nkw
                            do irx = -nkw, nkw

                                if (is_in_block(rgx + irx, rgy + iry, rgz + irz) &
                                        .and. ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then

                                    amp3 = sgmtr%recr(ir)%interp_ix(irx) &
                                        *sgmtr%recr(ir)%interp_iy(iry) &
                                        *sgmtr%recr(ir)%interp_hz(irz) &
                                        *sgmtr%recr(ir)%weight*0.25

                                    vxr_ixiyhz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vxr_ixiyhz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vx(t, ir)*amp3
                                    vyr_ixiyhz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vyr_ixiyhz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vy(t, ir)*amp3
                                    vzr_ixiyhz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vzr_ixiyhz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vz(t, ir)*amp3

                                end if

                            end do
                        end do
                    end do
                    !$omp end do

                    rgx = sgmtr%recr(ir)%hx
                    rgy = sgmtr%recr(ir)%hy
                    rgz = sgmtr%recr(ir)%hz
                    !$omp do collapse(3)
                    do irz = -nkw, nkw
                        do iry = -nkw, nkw
                            do irx = -nkw, nkw

                                if (is_in_block(rgx + irx, rgy + iry, rgz + irz) &
                                        .and. ifelse(yn_free_surface, rgz + irz >= 2, .true.)) then

                                    amp4 = sgmtr%recr(ir)%interp_hx(irx) &
                                        *sgmtr%recr(ir)%interp_hy(iry) &
                                        *sgmtr%recr(ir)%interp_hz(irz) &
                                        *sgmtr%recr(ir)%weight*0.25

                                    vxr_hxhyhz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vxr_hxhyhz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vx(t, ir)*amp4
                                    vyr_hxhyhz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vyr_hxhyhz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vy(t, ir)*amp4
                                    vzr_hxhyhz(rgx + irx, rgy + iry, rgz + irz) = &
                                        vzr_hxhyhz(rgx + irx, rgy + iry, rgz + irz) &
                                        + seis_vz(t, ir)*amp4

                                end if

                            end do
                        end do
                    end do
                    !$omp end do

                end if
            end do
            !$omp end parallel

            ! Compute gradients
            if (mod(t, cc_step_interval) == 0) then
                if (yn_update_medium .and. t >= sgmtr%srcr(1)%hnt) then
                    call compute_gradient
                end if
                if (yn_update_source) then
                    call compute_gradient_source(t)
                end if
            end if

            if (verbose .and. (mod(t, max(nint(nt/10.0), 1)) == 0 .or. t == 1 .or. t == nt)) then

                wnan = group_and(any(isnan(vxr_hxiyiz)) .or. any(isnan(vyr_hxiyiz)) .or. any(isnan(vzr_hxiyiz)))

                wmin1 = group_min(min(minval(vxr_hxiyiz), minval(vxr_ixhyiz), minval(vxr_ixiyhz), minval(vxr_hxhyhz)))
                wmax1 = group_max(max(maxval(vxr_hxiyiz), maxval(vxr_ixhyiz), maxval(vxr_ixiyhz), maxval(vxr_hxhyhz)))

                wmin2 = group_min(min(minval(vyr_hxiyiz), minval(vyr_ixhyiz), minval(vyr_ixiyhz), minval(vyr_hxhyhz)))
                wmax2 = group_max(max(maxval(vyr_hxiyiz), maxval(vyr_ixhyiz), maxval(vyr_ixiyhz), maxval(vyr_hxhyhz)))

                wmin3 = group_min(min(minval(vzr_hxiyiz), minval(vzr_ixhyiz), minval(vzr_ixiyhz), minval(vzr_hxhyhz)))
                wmax3 = group_max(max(maxval(vzr_hxiyiz), maxval(vzr_ixhyiz), maxval(vzr_ixiyhz), maxval(vzr_hxhyhz)))

                if (rankid_group == 0) then
                    call warn(date_time_compact()//' >> Shot '//num2str(sgmtr%id) &
                        //' FWI gradient computation step '//num2str(t)//' of '//num2str(nt))
                    if (wnan) then
                        call warn(date_time_compact()//' >> Vxr, Vyr, Vzr contain NaN!')
                        stop
                    else
                        call warn(date_time_compact()//' >> Vxr, Vyr, Vzr value range = ')
                        call warn(date_time_compact()//'      '//num2str(wmin1, '(es)')//' ~ '//num2str(wmax1, '(es)'))
                        call warn(date_time_compact()//'      '//num2str(wmin2, '(es)')//' ~ '//num2str(wmax2, '(es)'))
                        call warn(date_time_compact()//'      '//num2str(wmin3, '(es)')//' ~ '//num2str(wmax3, '(es)'))
                    end if
                end if

            end if

        end do

        ! Delete temporary files
        call close_boundary_saving(delete=.true.)

        ! Output source parameter gradient
        if (yn_update_source) then
            call allreduce_array_group(grad_mt)
            if (rankid_group == 0) then
                call grd%init(n=[nc_mt, 1, 1], d=[1.0, 1.0, 1.0], o=[0.0, 0.0, 0.0])
                grd%array = -reshape(grad_mt/maxval(abs(grad_mt)), [nc_mt, 1, 1])
                call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_mt.grd')
            end if
        end if
        call mpibarrier_group

        if (.not. yn_update_medium) then
            return
        end if

        ! Output medium parameter gradient
        if (yn_energy_precond) then

            call allreduce_array_group(grad_c11)
            call allreduce_array_group(grad_c12)
            call allreduce_array_group(grad_c13)
            call allreduce_array_group(grad_c14)
            call allreduce_array_group(grad_c15)
            call allreduce_array_group(grad_c16)
            call allreduce_array_group(grad_c22)
            call allreduce_array_group(grad_c23)
            call allreduce_array_group(grad_c24)
            call allreduce_array_group(grad_c25)
            call allreduce_array_group(grad_c26)
            call allreduce_array_group(grad_c33)
            call allreduce_array_group(grad_c34)
            call allreduce_array_group(grad_c35)
            call allreduce_array_group(grad_c36)
            call allreduce_array_group(grad_c44)
            call allreduce_array_group(grad_c45)
            call allreduce_array_group(grad_c46)
            call allreduce_array_group(grad_c55)
            call allreduce_array_group(grad_c56)
            call allreduce_array_group(grad_c66)
            call allreduce_array_group(grad_rho)

            call allreduce_array_group(energy_src_v)
            call allreduce_array_group(energy_rec_v)
            call allreduce_array_group(energy_src_a)
            call allreduce_array_group(energy_rec_a)

            if (kernel_v /= '') then
                energy_src_v = energy_src_v + 1.0e-3*maxval(energy_src_v)
                energy_rec_v = energy_rec_v + 1.0e-3*maxval(energy_rec_v)
                energy_src_v = sqrt(energy_src_v*energy_rec_v)
                grad_c11 = grad_c11/energy_src_v
                grad_c12 = grad_c12/energy_src_v
                grad_c13 = grad_c13/energy_src_v
                grad_c14 = grad_c14/energy_src_v
                grad_c15 = grad_c15/energy_src_v
                grad_c16 = grad_c16/energy_src_v
                grad_c22 = grad_c22/energy_src_v
                grad_c23 = grad_c23/energy_src_v
                grad_c24 = grad_c24/energy_src_v
                grad_c25 = grad_c25/energy_src_v
                grad_c26 = grad_c26/energy_src_v
                grad_c33 = grad_c33/energy_src_v
                grad_c34 = grad_c34/energy_src_v
                grad_c35 = grad_c35/energy_src_v
                grad_c36 = grad_c36/energy_src_v
                grad_c44 = grad_c44/energy_src_v
                grad_c45 = grad_c45/energy_src_v
                grad_c46 = grad_c46/energy_src_v
                grad_c55 = grad_c55/energy_src_v
                grad_c56 = grad_c56/energy_src_v
                grad_c66 = grad_c66/energy_src_v
            end if

            if (kernel_a /= '') then
                energy_src_a = energy_src_a + 1.0e-3*maxval(energy_src_a)
                energy_rec_a = energy_rec_a + 1.0e-3*maxval(energy_rec_a)
                energy_src_a = sqrt(energy_src_a*energy_rec_a)
                grad_rho = grad_rho/energy_src_a
            end if

        end if

        if (rankid_group == 0) then

            nx = this%nx
            ny = this%ny
            nz = this%nz

            ! For free-surface model, map computed gradients to regular mesh
            if (yn_free_surface) then
                if (kernel_v /= '') then
                    call map_irregular_to_regular(grad_c11, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c12, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c13, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c14, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c15, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c16, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c22, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c23, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c24, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c25, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c26, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c33, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c34, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c35, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c36, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c44, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c45, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c46, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c55, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c46, this, [1, nx, 1, ny, 1, nz])
                    call map_irregular_to_regular(grad_c56, this, [1, nx, 1, ny, 1, nz])
                end if
                if (kernel_a /= '') then
                    call map_irregular_to_regular(grad_rho, this, [1, nx, 1, ny, 1, nz])
                end if
            end if

            call grd%init(n=[nz, ny, nx], d=[dz, dy, dx], o=[oz, oy, ox])

            if (yn_free_surface) then
                rho = rho*1.0e-3
                call map_irregular_to_regular(rho, this, [1, nx, 1, ny, 1, nz])
            else
                rho = permute(this%rho, 321)*1.0e-3
            end if

            select case (aniso_param)

                case ('iso')
                    ! Using chain rule to compute gradients of Vp and Vs from
                    ! gradients of Cij. See Vigh et al. (2014)

                    vp = vp(1:nx, 1:ny, 1:nz)*1.0e-3
                    vs = vs(1:nx, 1:ny, 1:nz)*1.0e-3

                    grad_vp = grad_c11 + grad_c12 + grad_c13 + grad_c22 + grad_c23 + grad_c33
                    grad_vp = 2*rho*vp*grad_vp

                    grad_vs = grad_c44 + grad_c55 + grad_c66 - 2*grad_c12 - 2*grad_c13 - 2*grad_c23
                    grad_vs = 2*rho*vs*grad_vs

                    grd%array = permute(return_normal(grad_vp), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vp.grd')

                    grd%array = permute(return_normal(grad_vs), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vs.grd')

                case ('thomsen')

                    vp = vp(1:nx, 1:ny, 1:nz)*1.0e-3
                    vs = vs(1:nx, 1:ny, 1:nz)*1.0e-3
                    tieps = tieps(1:nx, 1:ny, 1:nz)
                    tidel = tidel(1:nx, 1:ny, 1:nz)
                    tigam = tigam(1:nx, 1:ny, 1:nz)
                    tithe = -tithe(1:nx, 1:ny, 1:nz)
                    tiphi = -tiphi(1:nx, 1:ny, 1:nz) ! Make angles consistent with what is used in module_anisotropy for VTI rotation

                    grad = &
                        + thomsen_dc11_dvp*grad_c11 &
                        + thomsen_dc12_dvp*grad_c12 &
                        + thomsen_dc13_dvp*grad_c13 &
                        + thomsen_dc14_dvp*grad_c14 &
                        + thomsen_dc15_dvp*grad_c15 &
                        + thomsen_dc16_dvp*grad_c16 &
                        + thomsen_dc22_dvp*grad_c22 &
                        + thomsen_dc23_dvp*grad_c23 &
                        + thomsen_dc24_dvp*grad_c24 &
                        + thomsen_dc25_dvp*grad_c25 &
                        + thomsen_dc26_dvp*grad_c26 &
                        + thomsen_dc33_dvp*grad_c33 &
                        + thomsen_dc34_dvp*grad_c34 &
                        + thomsen_dc35_dvp*grad_c35 &
                        + thomsen_dc36_dvp*grad_c36 &
                        + thomsen_dc44_dvp*grad_c44 &
                        + thomsen_dc45_dvp*grad_c45 &
                        + thomsen_dc46_dvp*grad_c46 &
                        + thomsen_dc55_dvp*grad_c55 &
                        + thomsen_dc56_dvp*grad_c56 &
                        + thomsen_dc66_dvp*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vp.grd')

                    grad = &
                        + thomsen_dc11_dvs*grad_c11 &
                        + thomsen_dc12_dvs*grad_c12 &
                        + thomsen_dc13_dvs*grad_c13 &
                        + thomsen_dc14_dvs*grad_c14 &
                        + thomsen_dc15_dvs*grad_c15 &
                        + thomsen_dc16_dvs*grad_c16 &
                        + thomsen_dc22_dvs*grad_c22 &
                        + thomsen_dc23_dvs*grad_c23 &
                        + thomsen_dc24_dvs*grad_c24 &
                        + thomsen_dc25_dvs*grad_c25 &
                        + thomsen_dc26_dvs*grad_c26 &
                        + thomsen_dc33_dvs*grad_c33 &
                        + thomsen_dc34_dvs*grad_c34 &
                        + thomsen_dc35_dvs*grad_c35 &
                        + thomsen_dc36_dvs*grad_c36 &
                        + thomsen_dc44_dvs*grad_c44 &
                        + thomsen_dc45_dvs*grad_c45 &
                        + thomsen_dc46_dvs*grad_c46 &
                        + thomsen_dc55_dvs*grad_c55 &
                        + thomsen_dc56_dvs*grad_c56 &
                        + thomsen_dc66_dvs*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vs.grd')

                    grad = &
                        + thomsen_dc11_deps*grad_c11 &
                        + thomsen_dc12_deps*grad_c12 &
                        + thomsen_dc13_deps*grad_c13 &
                        + thomsen_dc14_deps*grad_c14 &
                        + thomsen_dc15_deps*grad_c15 &
                        + thomsen_dc16_deps*grad_c16 &
                        + thomsen_dc22_deps*grad_c22 &
                        + thomsen_dc23_deps*grad_c23 &
                        + thomsen_dc24_deps*grad_c24 &
                        + thomsen_dc25_deps*grad_c25 &
                        + thomsen_dc26_deps*grad_c26 &
                        + thomsen_dc33_deps*grad_c33 &
                        + thomsen_dc34_deps*grad_c34 &
                        + thomsen_dc35_deps*grad_c35 &
                        + thomsen_dc36_deps*grad_c36 &
                        + thomsen_dc44_deps*grad_c44 &
                        + thomsen_dc45_deps*grad_c45 &
                        + thomsen_dc46_deps*grad_c46 &
                        + thomsen_dc55_deps*grad_c55 &
                        + thomsen_dc56_deps*grad_c56 &
                        + thomsen_dc66_deps*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_epsilon.grd')

                    grad = &
                        + thomsen_dc11_ddel*grad_c11 &
                        + thomsen_dc12_ddel*grad_c12 &
                        + thomsen_dc13_ddel*grad_c13 &
                        + thomsen_dc14_ddel*grad_c14 &
                        + thomsen_dc15_ddel*grad_c15 &
                        + thomsen_dc16_ddel*grad_c16 &
                        + thomsen_dc22_ddel*grad_c22 &
                        + thomsen_dc23_ddel*grad_c23 &
                        + thomsen_dc24_ddel*grad_c24 &
                        + thomsen_dc25_ddel*grad_c25 &
                        + thomsen_dc26_ddel*grad_c26 &
                        + thomsen_dc33_ddel*grad_c33 &
                        + thomsen_dc34_ddel*grad_c34 &
                        + thomsen_dc35_ddel*grad_c35 &
                        + thomsen_dc36_ddel*grad_c36 &
                        + thomsen_dc44_ddel*grad_c44 &
                        + thomsen_dc45_ddel*grad_c45 &
                        + thomsen_dc46_ddel*grad_c46 &
                        + thomsen_dc55_ddel*grad_c55 &
                        + thomsen_dc56_ddel*grad_c56 &
                        + thomsen_dc66_ddel*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_delta.grd')

                    grad = &
                        + thomsen_dc11_dgam*grad_c11 &
                        + thomsen_dc12_dgam*grad_c12 &
                        + thomsen_dc13_dgam*grad_c13 &
                        + thomsen_dc14_dgam*grad_c14 &
                        + thomsen_dc15_dgam*grad_c15 &
                        + thomsen_dc16_dgam*grad_c16 &
                        + thomsen_dc22_dgam*grad_c22 &
                        + thomsen_dc23_dgam*grad_c23 &
                        + thomsen_dc24_dgam*grad_c24 &
                        + thomsen_dc25_dgam*grad_c25 &
                        + thomsen_dc26_dgam*grad_c26 &
                        + thomsen_dc33_dgam*grad_c33 &
                        + thomsen_dc34_dgam*grad_c34 &
                        + thomsen_dc35_dgam*grad_c35 &
                        + thomsen_dc36_dgam*grad_c36 &
                        + thomsen_dc44_dgam*grad_c44 &
                        + thomsen_dc45_dgam*grad_c45 &
                        + thomsen_dc46_dgam*grad_c46 &
                        + thomsen_dc55_dgam*grad_c55 &
                        + thomsen_dc56_dgam*grad_c56 &
                        + thomsen_dc66_dgam*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_gamma.grd')


                case ('a-t')

                    vp = vp(1:nx, 1:ny, 1:nz)*1.0e-3
                    vs = vs(1:nx, 1:ny, 1:nz)*1.0e-3
                    tieps = tieps(1:nx, 1:ny, 1:nz)
                    tieta = tieta(1:nx, 1:ny, 1:nz)
                    tigam = tigam(1:nx, 1:ny, 1:nz)
                    tithe = -tithe(1:nx, 1:ny, 1:nz)
                    tiphi = -tiphi(1:nx, 1:ny, 1:nz) ! Make angles consistent with what is used in module_anisotropy for VTI rotation

                    grad = &
                        + alkhalifah_tsvankin_dc11_dvp*grad_c11 &
                        + alkhalifah_tsvankin_dc12_dvp*grad_c12 &
                        + alkhalifah_tsvankin_dc13_dvp*grad_c13 &
                        + alkhalifah_tsvankin_dc14_dvp*grad_c14 &
                        + alkhalifah_tsvankin_dc15_dvp*grad_c15 &
                        + alkhalifah_tsvankin_dc16_dvp*grad_c16 &
                        + alkhalifah_tsvankin_dc22_dvp*grad_c22 &
                        + alkhalifah_tsvankin_dc23_dvp*grad_c23 &
                        + alkhalifah_tsvankin_dc24_dvp*grad_c24 &
                        + alkhalifah_tsvankin_dc25_dvp*grad_c25 &
                        + alkhalifah_tsvankin_dc26_dvp*grad_c26 &
                        + alkhalifah_tsvankin_dc33_dvp*grad_c33 &
                        + alkhalifah_tsvankin_dc34_dvp*grad_c34 &
                        + alkhalifah_tsvankin_dc35_dvp*grad_c35 &
                        + alkhalifah_tsvankin_dc36_dvp*grad_c36 &
                        + alkhalifah_tsvankin_dc44_dvp*grad_c44 &
                        + alkhalifah_tsvankin_dc45_dvp*grad_c45 &
                        + alkhalifah_tsvankin_dc46_dvp*grad_c46 &
                        + alkhalifah_tsvankin_dc55_dvp*grad_c55 &
                        + alkhalifah_tsvankin_dc56_dvp*grad_c56 &
                        + alkhalifah_tsvankin_dc66_dvp*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vp.grd')

                    grad = &
                        + alkhalifah_tsvankin_dc11_dvs*grad_c11 &
                        + alkhalifah_tsvankin_dc12_dvs*grad_c12 &
                        + alkhalifah_tsvankin_dc13_dvs*grad_c13 &
                        + alkhalifah_tsvankin_dc14_dvs*grad_c14 &
                        + alkhalifah_tsvankin_dc15_dvs*grad_c15 &
                        + alkhalifah_tsvankin_dc16_dvs*grad_c16 &
                        + alkhalifah_tsvankin_dc22_dvs*grad_c22 &
                        + alkhalifah_tsvankin_dc23_dvs*grad_c23 &
                        + alkhalifah_tsvankin_dc24_dvs*grad_c24 &
                        + alkhalifah_tsvankin_dc25_dvs*grad_c25 &
                        + alkhalifah_tsvankin_dc26_dvs*grad_c26 &
                        + alkhalifah_tsvankin_dc33_dvs*grad_c33 &
                        + alkhalifah_tsvankin_dc34_dvs*grad_c34 &
                        + alkhalifah_tsvankin_dc35_dvs*grad_c35 &
                        + alkhalifah_tsvankin_dc36_dvs*grad_c36 &
                        + alkhalifah_tsvankin_dc44_dvs*grad_c44 &
                        + alkhalifah_tsvankin_dc45_dvs*grad_c45 &
                        + alkhalifah_tsvankin_dc46_dvs*grad_c46 &
                        + alkhalifah_tsvankin_dc55_dvs*grad_c55 &
                        + alkhalifah_tsvankin_dc56_dvs*grad_c56 &
                        + alkhalifah_tsvankin_dc66_dvs*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_vs.grd')

                    grad = &
                        + alkhalifah_tsvankin_dc11_deps*grad_c11 &
                        + alkhalifah_tsvankin_dc12_deps*grad_c12 &
                        + alkhalifah_tsvankin_dc13_deps*grad_c13 &
                        + alkhalifah_tsvankin_dc14_deps*grad_c14 &
                        + alkhalifah_tsvankin_dc15_deps*grad_c15 &
                        + alkhalifah_tsvankin_dc16_deps*grad_c16 &
                        + alkhalifah_tsvankin_dc22_deps*grad_c22 &
                        + alkhalifah_tsvankin_dc23_deps*grad_c23 &
                        + alkhalifah_tsvankin_dc24_deps*grad_c24 &
                        + alkhalifah_tsvankin_dc25_deps*grad_c25 &
                        + alkhalifah_tsvankin_dc26_deps*grad_c26 &
                        + alkhalifah_tsvankin_dc33_deps*grad_c33 &
                        + alkhalifah_tsvankin_dc34_deps*grad_c34 &
                        + alkhalifah_tsvankin_dc35_deps*grad_c35 &
                        + alkhalifah_tsvankin_dc36_deps*grad_c36 &
                        + alkhalifah_tsvankin_dc44_deps*grad_c44 &
                        + alkhalifah_tsvankin_dc45_deps*grad_c45 &
                        + alkhalifah_tsvankin_dc46_deps*grad_c46 &
                        + alkhalifah_tsvankin_dc55_deps*grad_c55 &
                        + alkhalifah_tsvankin_dc56_deps*grad_c56 &
                        + alkhalifah_tsvankin_dc66_deps*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_epsilon.grd')

                    grad = &
                        + alkhalifah_tsvankin_dc11_deta*grad_c11 &
                        + alkhalifah_tsvankin_dc12_deta*grad_c12 &
                        + alkhalifah_tsvankin_dc13_deta*grad_c13 &
                        + alkhalifah_tsvankin_dc14_deta*grad_c14 &
                        + alkhalifah_tsvankin_dc15_deta*grad_c15 &
                        + alkhalifah_tsvankin_dc16_deta*grad_c16 &
                        + alkhalifah_tsvankin_dc22_deta*grad_c22 &
                        + alkhalifah_tsvankin_dc23_deta*grad_c23 &
                        + alkhalifah_tsvankin_dc24_deta*grad_c24 &
                        + alkhalifah_tsvankin_dc25_deta*grad_c25 &
                        + alkhalifah_tsvankin_dc26_deta*grad_c26 &
                        + alkhalifah_tsvankin_dc33_deta*grad_c33 &
                        + alkhalifah_tsvankin_dc34_deta*grad_c34 &
                        + alkhalifah_tsvankin_dc35_deta*grad_c35 &
                        + alkhalifah_tsvankin_dc36_deta*grad_c36 &
                        + alkhalifah_tsvankin_dc44_deta*grad_c44 &
                        + alkhalifah_tsvankin_dc45_deta*grad_c45 &
                        + alkhalifah_tsvankin_dc46_deta*grad_c46 &
                        + alkhalifah_tsvankin_dc55_deta*grad_c55 &
                        + alkhalifah_tsvankin_dc56_deta*grad_c56 &
                        + alkhalifah_tsvankin_dc66_deta*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_eta.grd')

                    grad = &
                        + alkhalifah_tsvankin_dc11_dgam*grad_c11 &
                        + alkhalifah_tsvankin_dc12_dgam*grad_c12 &
                        + alkhalifah_tsvankin_dc13_dgam*grad_c13 &
                        + alkhalifah_tsvankin_dc14_dgam*grad_c14 &
                        + alkhalifah_tsvankin_dc15_dgam*grad_c15 &
                        + alkhalifah_tsvankin_dc16_dgam*grad_c16 &
                        + alkhalifah_tsvankin_dc22_dgam*grad_c22 &
                        + alkhalifah_tsvankin_dc23_dgam*grad_c23 &
                        + alkhalifah_tsvankin_dc24_dgam*grad_c24 &
                        + alkhalifah_tsvankin_dc25_dgam*grad_c25 &
                        + alkhalifah_tsvankin_dc26_dgam*grad_c26 &
                        + alkhalifah_tsvankin_dc33_dgam*grad_c33 &
                        + alkhalifah_tsvankin_dc34_dgam*grad_c34 &
                        + alkhalifah_tsvankin_dc35_dgam*grad_c35 &
                        + alkhalifah_tsvankin_dc36_dgam*grad_c36 &
                        + alkhalifah_tsvankin_dc44_dgam*grad_c44 &
                        + alkhalifah_tsvankin_dc45_dgam*grad_c45 &
                        + alkhalifah_tsvankin_dc46_dgam*grad_c46 &
                        + alkhalifah_tsvankin_dc55_dgam*grad_c55 &
                        + alkhalifah_tsvankin_dc56_dgam*grad_c56 &
                        + alkhalifah_tsvankin_dc66_dgam*grad_c66
                    grd%array = permute(return_normal(grad), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_gamma.grd')

                case ('cij')
                    grd%array = permute(return_normal(grad_c11), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c11.grd')
                    grd%array = permute(return_normal(grad_c12), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c12.grd')
                    grd%array = permute(return_normal(grad_c13), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c13.grd')
                    grd%array = permute(return_normal(grad_c14), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c14.grd')
                    grd%array = permute(return_normal(grad_c15), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c15.grd')
                    grd%array = permute(return_normal(grad_c16), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c16.grd')
                    grd%array = permute(return_normal(grad_c22), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c22.grd')
                    grd%array = permute(return_normal(grad_c23), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c23.grd')
                    grd%array = permute(return_normal(grad_c24), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c24.grd')
                    grd%array = permute(return_normal(grad_c25), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c25.grd')
                    grd%array = permute(return_normal(grad_c26), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c26.grd')
                    grd%array = permute(return_normal(grad_c33), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c33.grd')
                    grd%array = permute(return_normal(grad_c34), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c34.grd')
                    grd%array = permute(return_normal(grad_c35), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c35.grd')
                    grd%array = permute(return_normal(grad_c36), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c36.grd')
                    grd%array = permute(return_normal(grad_c44), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c44.grd')
                    grd%array = permute(return_normal(grad_c45), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c45.grd')
                    grd%array = permute(return_normal(grad_c46), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c46.grd')
                    grd%array = permute(return_normal(grad_c55), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c55.grd')
                    grd%array = permute(return_normal(grad_c56), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c56.grd')
                    grd%array = permute(return_normal(grad_c66), 321)
                    call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_c66.grd')

            end select

            grd%array = permute(return_normal(grad_rho), 321)
            call grd%output(tidy(dir_working)//'/shot_'//num2str(sgmtr%id)//'_grad_rho.grd')

        end if

        call mpibarrier_group

    end subroutine

    !
    !> Compute medium parameter gradients
    !
    subroutine compute_gradient

        integer :: i, j, k
        integer :: sgnh
        real :: pxx, pxy, pxz, pyy, pyz, pzz
        real :: pxxr, pxyr, pxzr, pyyr, pyzr, pzzr

        !$omp parallel do private(i, j, k, pxx, pxy, pxz, pyy, pyz, pzz, &
            !$omp   pxxr, pxyr, pxzr, pyyr, pyzr, pzzr) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1, nx2

                    ! Forward wavefields
                    pxx = &
                        stressxx_ixiyiz(i, j, k) &
                        + 0.25*sum(stressxx_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stressxx_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stressxx_ixhyhz(i, j:j + 1, k:k + 1)) &
                        - prev_stressxx(i, j, k)
                    pyy = &
                        stressyy_ixiyiz(i, j, k) &
                        + 0.25*sum(stressyy_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stressyy_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stressyy_ixhyhz(i, j:j + 1, k:k + 1)) &
                        - prev_stressyy(i, j, k)
                    pzz = &
                        stresszz_ixiyiz(i, j, k) &
                        + 0.25*sum(stresszz_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stresszz_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stresszz_ixhyhz(i, j:j + 1, k:k + 1)) &
                        - prev_stresszz(i, j, k)
                    pyz = &
                        stressyz_ixiyiz(i, j, k) &
                        + 0.25*sum(stressyz_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stressyz_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stressyz_ixhyhz(i, j:j + 1, k:k + 1)) &
                        - prev_stressyz(i, j, k)
                    pxz = &
                        stressxz_ixiyiz(i, j, k) &
                        + 0.25*sum(stressxz_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stressxz_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stressxz_ixhyhz(i, j:j + 1, k:k + 1)) &
                        - prev_stressxz(i, j, k)
                    pxy = &
                        stressxy_ixiyiz(i, j, k) &
                        + 0.25*sum(stressxy_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stressxy_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stressxy_ixhyhz(i, j:j + 1, k:k + 1)) &
                        - prev_stressxy(i, j, k)

                    strainxx(i, j, k) = &
                        s11(i, j, k)*pxx &
                        + s12(i, j, k)*pyy &
                        + s13(i, j, k)*pzz &
                        + s14(i, j, k)*pyz &
                        + s15(i, j, k)*pxz &
                        + s16(i, j, k)*pxy
                    strainyy(i, j, k) = &
                        s12(i, j, k)*pxx &
                        + s22(i, j, k)*pyy &
                        + s23(i, j, k)*pzz &
                        + s24(i, j, k)*pyz &
                        + s25(i, j, k)*pxz &
                        + s26(i, j, k)*pxy
                    strainzz(i, j, k) = &
                        s13(i, j, k)*pxx &
                        + s23(i, j, k)*pyy &
                        + s33(i, j, k)*pzz &
                        + s34(i, j, k)*pyz &
                        + s35(i, j, k)*pxz &
                        + s36(i, j, k)*pxy
                    strainyz(i, j, k) = &
                        s14(i, j, k)*pxx &
                        + s24(i, j, k)*pyy &
                        + s34(i, j, k)*pzz &
                        + s44(i, j, k)*pyz &
                        + s45(i, j, k)*pxz &
                        + s46(i, j, k)*pxy
                    strainxz(i, j, k) = &
                        s15(i, j, k)*pxx &
                        + s25(i, j, k)*pyy &
                        + s35(i, j, k)*pzz &
                        + s45(i, j, k)*pyz &
                        + s55(i, j, k)*pxz &
                        + s56(i, j, k)*pxy
                    strainxy(i, j, k) = &
                        s16(i, j, k)*pxx &
                        + s26(i, j, k)*pyy &
                        + s36(i, j, k)*pzz &
                        + s46(i, j, k)*pyz &
                        + s56(i, j, k)*pxz &
                        + s66(i, j, k)*pxy

                    ! Adjoint wavefields
                    pxxr = &
                        stressxxr_ixiyiz(i, j, k) &
                        + 0.25*sum(stressxxr_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stressxxr_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stressxxr_ixhyhz(i, j:j + 1, k:k + 1))
                    pyyr = &
                        stressyyr_ixiyiz(i, j, k) &
                        + 0.25*sum(stressyyr_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stressyyr_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stressyyr_ixhyhz(i, j:j + 1, k:k + 1))
                    pzzr = &
                        stresszzr_ixiyiz(i, j, k) &
                        + 0.25*sum(stresszzr_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stresszzr_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stresszzr_ixhyhz(i, j:j + 1, k:k + 1))
                    pyzr = &
                        stressyzr_ixiyiz(i, j, k) &
                        + 0.25*sum(stressyzr_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stressyzr_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stressyzr_ixhyhz(i, j:j + 1, k:k + 1))
                    pxzr = &
                        stressxzr_ixiyiz(i, j, k) &
                        + 0.25*sum(stressxzr_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stressxzr_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stressxzr_ixhyhz(i, j:j + 1, k:k + 1))
                    pxyr = &
                        stressxyr_ixiyiz(i, j, k) &
                        + 0.25*sum(stressxyr_hxhyiz(i:i + 1, j:j + 1, k)) &
                        + 0.25*sum(stressxyr_hxiyhz(i:i + 1, j, k:k + 1)) &
                        + 0.25*sum(stressxyr_ixhyhz(i, j:j + 1, k:k + 1))

                    strainxxr(i, j, k) = &
                        s11(i, j, k)*pxxr &
                        + s12(i, j, k)*pyyr &
                        + s13(i, j, k)*pzzr &
                        + s14(i, j, k)*pyzr &
                        + s15(i, j, k)*pxzr &
                        + s16(i, j, k)*pxyr
                    strainyyr(i, j, k) = &
                        s12(i, j, k)*pxxr &
                        + s22(i, j, k)*pyyr &
                        + s23(i, j, k)*pzzr &
                        + s24(i, j, k)*pyzr &
                        + s25(i, j, k)*pxzr &
                        + s26(i, j, k)*pxyr
                    strainzzr(i, j, k) = &
                        s13(i, j, k)*pxxr &
                        + s23(i, j, k)*pyyr &
                        + s33(i, j, k)*pzzr &
                        + s34(i, j, k)*pyzr &
                        + s35(i, j, k)*pxzr &
                        + s36(i, j, k)*pxyr
                    strainyzr(i, j, k) = &
                        s14(i, j, k)*pxxr &
                        + s24(i, j, k)*pyyr &
                        + s34(i, j, k)*pzzr &
                        + s44(i, j, k)*pyzr &
                        + s45(i, j, k)*pxzr &
                        + s46(i, j, k)*pxyr
                    strainxzr(i, j, k) = &
                        s15(i, j, k)*pxxr &
                        + s25(i, j, k)*pyyr &
                        + s35(i, j, k)*pzzr &
                        + s45(i, j, k)*pyzr &
                        + s55(i, j, k)*pxzr &
                        + s56(i, j, k)*pxyr
                    strainxyr(i, j, k) = &
                        s16(i, j, k)*pxxr &
                        + s26(i, j, k)*pyyr &
                        + s36(i, j, k)*pzzr &
                        + s46(i, j, k)*pyzr &
                        + s56(i, j, k)*pxzr &
                        + s66(i, j, k)*pxyr

                end do
            end do
        end do
        !$omp end parallel do

        if (kernel_v /= '') then

            if (kernel_v == 'full') then

                !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
                do k = nz1_interior, nz2_interior
                    do j = ny1_interior, ny2_interior
                        do i = nx1_interior, nx2_interior

                            grad_c11(i, j, k) = grad_c11(i, j, k) - strainxx(i, j, k)*strainxxr(i, j, k)
                            grad_c12(i, j, k) = grad_c12(i, j, k) - strainxx(i, j, k)*strainyyr(i, j, k) - strainyy(i, j, k)*strainxxr(i, j, k)
                            grad_c13(i, j, k) = grad_c13(i, j, k) - strainxx(i, j, k)*strainzzr(i, j, k) - strainzz(i, j, k)*strainxxr(i, j, k)
                            grad_c14(i, j, k) = grad_c14(i, j, k) - strainxx(i, j, k)*strainyzr(i, j, k) - strainyz(i, j, k)*strainxxr(i, j, k)
                            grad_c15(i, j, k) = grad_c15(i, j, k) - strainxx(i, j, k)*strainxzr(i, j, k) - strainxz(i, j, k)*strainxxr(i, j, k)
                            grad_c16(i, j, k) = grad_c16(i, j, k) - strainxx(i, j, k)*strainxyr(i, j, k) - strainxy(i, j, k)*strainxxr(i, j, k)
                            grad_c22(i, j, k) = grad_c22(i, j, k) - strainyy(i, j, k)*strainyyr(i, j, k)
                            grad_c23(i, j, k) = grad_c23(i, j, k) - strainyy(i, j, k)*strainzzr(i, j, k) - strainzz(i, j, k)*strainyyr(i, j, k)
                            grad_c24(i, j, k) = grad_c24(i, j, k) - strainyy(i, j, k)*strainyzr(i, j, k) - strainyz(i, j, k)*strainyyr(i, j, k)
                            grad_c25(i, j, k) = grad_c25(i, j, k) - strainyy(i, j, k)*strainxzr(i, j, k) - strainxz(i, j, k)*strainyyr(i, j, k)
                            grad_c26(i, j, k) = grad_c26(i, j, k) - strainyy(i, j, k)*strainxyr(i, j, k) - strainxy(i, j, k)*strainyyr(i, j, k)
                            grad_c33(i, j, k) = grad_c33(i, j, k) - strainzz(i, j, k)*strainzzr(i, j, k)
                            grad_c34(i, j, k) = grad_c34(i, j, k) - strainzz(i, j, k)*strainyzr(i, j, k) - strainyz(i, j, k)*strainzzr(i, j, k)
                            grad_c35(i, j, k) = grad_c35(i, j, k) - strainzz(i, j, k)*strainxzr(i, j, k) - strainxz(i, j, k)*strainzzr(i, j, k)
                            grad_c36(i, j, k) = grad_c36(i, j, k) - strainzz(i, j, k)*strainxyr(i, j, k) - strainxy(i, j, k)*strainzzr(i, j, k)
                            grad_c44(i, j, k) = grad_c44(i, j, k) - strainyz(i, j, k)*strainyzr(i, j, k)
                            grad_c45(i, j, k) = grad_c45(i, j, k) - strainyz(i, j, k)*strainxzr(i, j, k) - strainxz(i, j, k)*strainyzr(i, j, k)
                            grad_c46(i, j, k) = grad_c46(i, j, k) - strainyz(i, j, k)*strainxyr(i, j, k) - strainxy(i, j, k)*strainyzr(i, j, k)
                            grad_c55(i, j, k) = grad_c55(i, j, k) - strainxz(i, j, k)*strainxzr(i, j, k)
                            grad_c56(i, j, k) = grad_c56(i, j, k) - strainxz(i, j, k)*strainxyr(i, j, k) - strainxy(i, j, k)*strainxzr(i, j, k)
                            grad_c66(i, j, k) = grad_c66(i, j, k) - strainxy(i, j, k)*strainxyr(i, j, k)

                        end do
                    end do
                end do
                !$omp end parallel do

            else

                ! Along x
                if (index(kernel_v, 'lowx') /= 0) then
                    sgnh = 1
                else if (index(kernel_v, 'highx') /= 0) then
                    sgnh = -1
                else
                    sgnh = 0
                end if

                if (sgnh /= 0) then

                    call commute_array(strainxx, htlen, dim=1)
                    call commute_array(strainyy, htlen, dim=1)
                    call commute_array(strainzz, htlen, dim=1)
                    call commute_array(strainyz, htlen, dim=1)
                    call commute_array(strainxz, htlen, dim=1)
                    call commute_array(strainxy, htlen, dim=1)

                    call commute_array(strainxxr, htlen, dim=1)
                    call commute_array(strainyyr, htlen, dim=1)
                    call commute_array(strainzzr, htlen, dim=1)
                    call commute_array(strainyzr, htlen, dim=1)
                    call commute_array(strainxzr, htlen, dim=1)
                    call commute_array(strainxyr, htlen, dim=1)

                    strainxx_hilbert = compute_hilbert_transform(strainxx, dim=1)
                    strainyy_hilbert = compute_hilbert_transform(strainyy, dim=1)
                    strainzz_hilbert = compute_hilbert_transform(strainzz, dim=1)
                    strainyz_hilbert = compute_hilbert_transform(strainyz, dim=1)
                    strainxz_hilbert = compute_hilbert_transform(strainxz, dim=1)
                    strainxy_hilbert = compute_hilbert_transform(strainxy, dim=1)

                    strainxxr_hilbert = compute_hilbert_transform(strainxxr, dim=1)
                    strainyyr_hilbert = compute_hilbert_transform(strainyyr, dim=1)
                    strainzzr_hilbert = compute_hilbert_transform(strainzzr, dim=1)
                    strainyzr_hilbert = compute_hilbert_transform(strainyzr, dim=1)
                    strainxzr_hilbert = compute_hilbert_transform(strainxzr, dim=1)
                    strainxyr_hilbert = compute_hilbert_transform(strainxyr, dim=1)

                    grad_c11(interior_region) = grad_c11(interior_region) &
                        - compute_directional_gradient(strainxx, strainxxr, strainxx_hilbert, strainxxr_hilbert, sgnh, dim=1)
                    grad_c12(interior_region) = grad_c12(interior_region) &
                        - compute_directional_gradient(strainxx, strainyyr, strainxx_hilbert, strainyyr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainyy, strainxxr, strainyy_hilbert, strainxxr_hilbert, sgnh, dim=1)
                    grad_c13(interior_region) = grad_c13(interior_region) &
                        - compute_directional_gradient(strainxx, strainzzr, strainxx_hilbert, strainzzr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainzz, strainxxr, strainzz_hilbert, strainxxr_hilbert, sgnh, dim=1)
                    grad_c14(interior_region) = grad_c14(interior_region) &
                        - compute_directional_gradient(strainxx, strainyzr, strainxx_hilbert, strainyzr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainyz, strainxxr, strainyz_hilbert, strainxxr_hilbert, sgnh, dim=1)
                    grad_c15(interior_region) = grad_c15(interior_region) &
                        - compute_directional_gradient(strainxx, strainxzr, strainxx_hilbert, strainxzr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainxz, strainxxr, strainxz_hilbert, strainxxr_hilbert, sgnh, dim=1)
                    grad_c16(interior_region) = grad_c16(interior_region) &
                        - compute_directional_gradient(strainxx, strainxyr, strainxx_hilbert, strainxyr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainxy, strainxxr, strainxy_hilbert, strainxxr_hilbert, sgnh, dim=1)

                    grad_c22(interior_region) = grad_c22(interior_region) &
                        - compute_directional_gradient(strainyy, strainyyr, strainyy_hilbert, strainyyr_hilbert, sgnh, dim=1)
                    grad_c23(interior_region) = grad_c23(interior_region) &
                        - compute_directional_gradient(strainyy, strainzzr, strainyy_hilbert, strainzzr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainzz, strainyyr, strainzz_hilbert, strainyyr_hilbert, sgnh, dim=1)
                    grad_c24(interior_region) = grad_c24(interior_region) &
                        - compute_directional_gradient(strainyy, strainyzr, strainyy_hilbert, strainyzr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainyz, strainyyr, strainyz_hilbert, strainyyr_hilbert, sgnh, dim=1)
                    grad_c25(interior_region) = grad_c25(interior_region) &
                        - compute_directional_gradient(strainyy, strainxzr, strainyy_hilbert, strainxzr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainxz, strainyyr, strainxz_hilbert, strainyyr_hilbert, sgnh, dim=1)
                    grad_c26(interior_region) = grad_c26(interior_region) &
                        - compute_directional_gradient(strainyy, strainxyr, strainyy_hilbert, strainxyr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainxy, strainyyr, strainxy_hilbert, strainyyr_hilbert, sgnh, dim=1)

                    grad_c33(interior_region) = grad_c33(interior_region) &
                        - compute_directional_gradient(strainzz, strainzzr, strainzz_hilbert, strainzzr_hilbert, sgnh, dim=1)
                    grad_c34(interior_region) = grad_c34(interior_region) &
                        - compute_directional_gradient(strainzz, strainyzr, strainzz_hilbert, strainyzr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainyz, strainzzr, strainyz_hilbert, strainzzr_hilbert, sgnh, dim=1)
                    grad_c35(interior_region) = grad_c35(interior_region) &
                        - compute_directional_gradient(strainzz, strainxzr, strainzz_hilbert, strainxzr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainxz, strainzzr, strainxz_hilbert, strainzzr_hilbert, sgnh, dim=1)
                    grad_c36(interior_region) = grad_c36(interior_region) &
                        - compute_directional_gradient(strainzz, strainxyr, strainzz_hilbert, strainxyr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainxy, strainzzr, strainxy_hilbert, strainzzr_hilbert, sgnh, dim=1)

                    grad_c44(interior_region) = grad_c44(interior_region) &
                        - compute_directional_gradient(strainyz, strainyzr, strainyz_hilbert, strainyzr_hilbert, sgnh, dim=1)
                    grad_c45(interior_region) = grad_c45(interior_region) &
                        - compute_directional_gradient(strainyz, strainxzr, strainyz_hilbert, strainxzr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainxz, strainyzr, strainxz_hilbert, strainyzr_hilbert, sgnh, dim=1)
                    grad_c46(interior_region) = grad_c46(interior_region) &
                        - compute_directional_gradient(strainyz, strainxyr, strainyz_hilbert, strainxyr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainxy, strainyzr, strainxy_hilbert, strainyzr_hilbert, sgnh, dim=1)

                    grad_c55(interior_region) = grad_c55(interior_region) &
                        - compute_directional_gradient(strainxz, strainxzr, strainxz_hilbert, strainxzr_hilbert, sgnh, dim=1)
                    grad_c56(interior_region) = grad_c56(interior_region) &
                        - compute_directional_gradient(strainxz, strainxyr, strainxz_hilbert, strainxyr_hilbert, sgnh, dim=1) &
                        - compute_directional_gradient(strainxy, strainxzr, strainxy_hilbert, strainxzr_hilbert, sgnh, dim=1)

                    grad_c66(interior_region) = grad_c66(interior_region) &
                        - compute_directional_gradient(strainxy, strainxyr, strainxy_hilbert, strainxyr_hilbert, sgnh, dim=1)

                end if

                ! Along y
                if (index(kernel_v, 'lowy') /= 0) then
                    sgnh = 1
                else if (index(kernel_v, 'highy') /= 0) then
                    sgnh = -1
                else
                    sgnh = 0
                end if

                if (sgnh /= 0) then

                    call commute_array(strainxx, htlen, dim=2)
                    call commute_array(strainyy, htlen, dim=2)
                    call commute_array(strainzz, htlen, dim=2)
                    call commute_array(strainyz, htlen, dim=2)
                    call commute_array(strainxz, htlen, dim=2)
                    call commute_array(strainxy, htlen, dim=2)

                    call commute_array(strainxxr, htlen, dim=2)
                    call commute_array(strainyyr, htlen, dim=2)
                    call commute_array(strainzzr, htlen, dim=2)
                    call commute_array(strainyzr, htlen, dim=2)
                    call commute_array(strainxzr, htlen, dim=2)
                    call commute_array(strainxyr, htlen, dim=2)

                    strainxx_hilbert = compute_hilbert_transform(strainxx, dim=2)
                    strainyy_hilbert = compute_hilbert_transform(strainyy, dim=2)
                    strainzz_hilbert = compute_hilbert_transform(strainzz, dim=2)
                    strainyz_hilbert = compute_hilbert_transform(strainyz, dim=2)
                    strainxz_hilbert = compute_hilbert_transform(strainxz, dim=2)
                    strainxy_hilbert = compute_hilbert_transform(strainxy, dim=2)

                    strainxxr_hilbert = compute_hilbert_transform(strainxxr, dim=2)
                    strainyyr_hilbert = compute_hilbert_transform(strainyyr, dim=2)
                    strainzzr_hilbert = compute_hilbert_transform(strainzzr, dim=2)
                    strainyzr_hilbert = compute_hilbert_transform(strainyzr, dim=2)
                    strainxzr_hilbert = compute_hilbert_transform(strainxzr, dim=2)
                    strainxyr_hilbert = compute_hilbert_transform(strainxyr, dim=2)

                    grad_c11(interior_region) = grad_c11(interior_region) &
                        - compute_directional_gradient(strainxx, strainxxr, strainxx_hilbert, strainxxr_hilbert, sgnh, dim=2)
                    grad_c12(interior_region) = grad_c12(interior_region) &
                        - compute_directional_gradient(strainxx, strainyyr, strainxx_hilbert, strainyyr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainyy, strainxxr, strainyy_hilbert, strainxxr_hilbert, sgnh, dim=2)
                    grad_c13(interior_region) = grad_c13(interior_region) &
                        - compute_directional_gradient(strainxx, strainzzr, strainxx_hilbert, strainzzr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainzz, strainxxr, strainzz_hilbert, strainxxr_hilbert, sgnh, dim=2)
                    grad_c14(interior_region) = grad_c14(interior_region) &
                        - compute_directional_gradient(strainxx, strainyzr, strainxx_hilbert, strainyzr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainyz, strainxxr, strainyz_hilbert, strainxxr_hilbert, sgnh, dim=2)
                    grad_c15(interior_region) = grad_c15(interior_region) &
                        - compute_directional_gradient(strainxx, strainxzr, strainxx_hilbert, strainxzr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainxz, strainxxr, strainxz_hilbert, strainxxr_hilbert, sgnh, dim=2)
                    grad_c16(interior_region) = grad_c16(interior_region) &
                        - compute_directional_gradient(strainxx, strainxyr, strainxx_hilbert, strainxyr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainxy, strainxxr, strainxy_hilbert, strainxxr_hilbert, sgnh, dim=2)

                    grad_c22(interior_region) = grad_c22(interior_region) &
                        - compute_directional_gradient(strainyy, strainyyr, strainyy_hilbert, strainyyr_hilbert, sgnh, dim=2)
                    grad_c23(interior_region) = grad_c23(interior_region) &
                        - compute_directional_gradient(strainyy, strainzzr, strainyy_hilbert, strainzzr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainzz, strainyyr, strainzz_hilbert, strainyyr_hilbert, sgnh, dim=2)
                    grad_c24(interior_region) = grad_c24(interior_region) &
                        - compute_directional_gradient(strainyy, strainyzr, strainyy_hilbert, strainyzr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainyz, strainyyr, strainyz_hilbert, strainyyr_hilbert, sgnh, dim=2)
                    grad_c25(interior_region) = grad_c25(interior_region) &
                        - compute_directional_gradient(strainyy, strainxzr, strainyy_hilbert, strainxzr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainxz, strainyyr, strainxz_hilbert, strainyyr_hilbert, sgnh, dim=2)
                    grad_c26(interior_region) = grad_c26(interior_region) &
                        - compute_directional_gradient(strainyy, strainxyr, strainyy_hilbert, strainxyr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainxy, strainyyr, strainxy_hilbert, strainyyr_hilbert, sgnh, dim=2)

                    grad_c33(interior_region) = grad_c33(interior_region) &
                        - compute_directional_gradient(strainzz, strainzzr, strainzz_hilbert, strainzzr_hilbert, sgnh, dim=2)
                    grad_c34(interior_region) = grad_c34(interior_region) &
                        - compute_directional_gradient(strainzz, strainyzr, strainzz_hilbert, strainyzr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainyz, strainzzr, strainyz_hilbert, strainzzr_hilbert, sgnh, dim=2)
                    grad_c35(interior_region) = grad_c35(interior_region) &
                        - compute_directional_gradient(strainzz, strainxzr, strainzz_hilbert, strainxzr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainxz, strainzzr, strainxz_hilbert, strainzzr_hilbert, sgnh, dim=2)
                    grad_c36(interior_region) = grad_c36(interior_region) &
                        - compute_directional_gradient(strainzz, strainxyr, strainzz_hilbert, strainxyr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainxy, strainzzr, strainxy_hilbert, strainzzr_hilbert, sgnh, dim=2)

                    grad_c44(interior_region) = grad_c44(interior_region) &
                        - compute_directional_gradient(strainyz, strainyzr, strainyz_hilbert, strainyzr_hilbert, sgnh, dim=2)
                    grad_c45(interior_region) = grad_c45(interior_region) &
                        - compute_directional_gradient(strainyz, strainxzr, strainyz_hilbert, strainxzr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainxz, strainyzr, strainxz_hilbert, strainyzr_hilbert, sgnh, dim=2)
                    grad_c46(interior_region) = grad_c46(interior_region) &
                        - compute_directional_gradient(strainyz, strainxyr, strainyz_hilbert, strainxyr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainxy, strainyzr, strainxy_hilbert, strainyzr_hilbert, sgnh, dim=2)

                    grad_c55(interior_region) = grad_c55(interior_region) &
                        - compute_directional_gradient(strainxz, strainxzr, strainxz_hilbert, strainxzr_hilbert, sgnh, dim=2)
                    grad_c56(interior_region) = grad_c56(interior_region) &
                        - compute_directional_gradient(strainxz, strainxyr, strainxz_hilbert, strainxyr_hilbert, sgnh, dim=2) &
                        - compute_directional_gradient(strainxy, strainxzr, strainxy_hilbert, strainxzr_hilbert, sgnh, dim=2)

                    grad_c66(interior_region) = grad_c66(interior_region) &
                        - compute_directional_gradient(strainxy, strainxyr, strainxy_hilbert, strainxyr_hilbert, sgnh, dim=2)

                end if

                ! Along z
                if (index(kernel_v, 'lowz') /= 0) then
                    sgnh = 1
                else if (index(kernel_v, 'highz') /= 0) then
                    sgnh = -1
                else
                    sgnh = 0
                end if

                if (sgnh /= 0) then

                    call commute_array(strainxx, htlen, dim=3)
                    call commute_array(strainyy, htlen, dim=3)
                    call commute_array(strainzz, htlen, dim=3)
                    call commute_array(strainyz, htlen, dim=3)
                    call commute_array(strainxz, htlen, dim=3)
                    call commute_array(strainxy, htlen, dim=3)

                    call commute_array(strainxxr, htlen, dim=3)
                    call commute_array(strainyyr, htlen, dim=3)
                    call commute_array(strainzzr, htlen, dim=3)
                    call commute_array(strainyzr, htlen, dim=3)
                    call commute_array(strainxzr, htlen, dim=3)
                    call commute_array(strainxyr, htlen, dim=3)

                    strainxx_hilbert = compute_hilbert_transform(strainxx, dim=3)
                    strainyy_hilbert = compute_hilbert_transform(strainyy, dim=3)
                    strainzz_hilbert = compute_hilbert_transform(strainzz, dim=3)
                    strainyz_hilbert = compute_hilbert_transform(strainyz, dim=3)
                    strainxz_hilbert = compute_hilbert_transform(strainxz, dim=3)
                    strainxy_hilbert = compute_hilbert_transform(strainxy, dim=3)

                    strainxxr_hilbert = compute_hilbert_transform(strainxxr, dim=3)
                    strainyyr_hilbert = compute_hilbert_transform(strainyyr, dim=3)
                    strainzzr_hilbert = compute_hilbert_transform(strainzzr, dim=3)
                    strainyzr_hilbert = compute_hilbert_transform(strainyzr, dim=3)
                    strainxzr_hilbert = compute_hilbert_transform(strainxzr, dim=3)
                    strainxyr_hilbert = compute_hilbert_transform(strainxyr, dim=3)

                    grad_c11(interior_region) = grad_c11(interior_region) &
                        - compute_directional_gradient(strainxx, strainxxr, strainxx_hilbert, strainxxr_hilbert, sgnh, dim=3)
                    grad_c12(interior_region) = grad_c12(interior_region) &
                        - compute_directional_gradient(strainxx, strainyyr, strainxx_hilbert, strainyyr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainyy, strainxxr, strainyy_hilbert, strainxxr_hilbert, sgnh, dim=3)
                    grad_c13(interior_region) = grad_c13(interior_region) &
                        - compute_directional_gradient(strainxx, strainzzr, strainxx_hilbert, strainzzr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainzz, strainxxr, strainzz_hilbert, strainxxr_hilbert, sgnh, dim=3)
                    grad_c14(interior_region) = grad_c14(interior_region) &
                        - compute_directional_gradient(strainxx, strainyzr, strainxx_hilbert, strainyzr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainyz, strainxxr, strainyz_hilbert, strainxxr_hilbert, sgnh, dim=3)
                    grad_c15(interior_region) = grad_c15(interior_region) &
                        - compute_directional_gradient(strainxx, strainxzr, strainxx_hilbert, strainxzr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainxz, strainxxr, strainxz_hilbert, strainxxr_hilbert, sgnh, dim=3)
                    grad_c16(interior_region) = grad_c16(interior_region) &
                        - compute_directional_gradient(strainxx, strainxyr, strainxx_hilbert, strainxyr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainxy, strainxxr, strainxy_hilbert, strainxxr_hilbert, sgnh, dim=3)

                    grad_c22(interior_region) = grad_c22(interior_region) &
                        - compute_directional_gradient(strainyy, strainyyr, strainyy_hilbert, strainyyr_hilbert, sgnh, dim=3)
                    grad_c23(interior_region) = grad_c23(interior_region) &
                        - compute_directional_gradient(strainyy, strainzzr, strainyy_hilbert, strainzzr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainzz, strainyyr, strainzz_hilbert, strainyyr_hilbert, sgnh, dim=3)
                    grad_c24(interior_region) = grad_c24(interior_region) &
                        - compute_directional_gradient(strainyy, strainyzr, strainyy_hilbert, strainyzr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainyz, strainyyr, strainyz_hilbert, strainyyr_hilbert, sgnh, dim=3)
                    grad_c25(interior_region) = grad_c25(interior_region) &
                        - compute_directional_gradient(strainyy, strainxzr, strainyy_hilbert, strainxzr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainxz, strainyyr, strainxz_hilbert, strainyyr_hilbert, sgnh, dim=3)
                    grad_c26(interior_region) = grad_c26(interior_region) &
                        - compute_directional_gradient(strainyy, strainxyr, strainyy_hilbert, strainxyr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainxy, strainyyr, strainxy_hilbert, strainyyr_hilbert, sgnh, dim=3)

                    grad_c33(interior_region) = grad_c33(interior_region) &
                        - compute_directional_gradient(strainzz, strainzzr, strainzz_hilbert, strainzzr_hilbert, sgnh, dim=3)
                    grad_c34(interior_region) = grad_c34(interior_region) &
                        - compute_directional_gradient(strainzz, strainyzr, strainzz_hilbert, strainyzr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainyz, strainzzr, strainyz_hilbert, strainzzr_hilbert, sgnh, dim=3)
                    grad_c35(interior_region) = grad_c35(interior_region) &
                        - compute_directional_gradient(strainzz, strainxzr, strainzz_hilbert, strainxzr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainxz, strainzzr, strainxz_hilbert, strainzzr_hilbert, sgnh, dim=3)
                    grad_c36(interior_region) = grad_c36(interior_region) &
                        - compute_directional_gradient(strainzz, strainxyr, strainzz_hilbert, strainxyr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainxy, strainzzr, strainxy_hilbert, strainzzr_hilbert, sgnh, dim=3)

                    grad_c44(interior_region) = grad_c44(interior_region) &
                        - compute_directional_gradient(strainyz, strainyzr, strainyz_hilbert, strainyzr_hilbert, sgnh, dim=3)
                    grad_c45(interior_region) = grad_c45(interior_region) &
                        - compute_directional_gradient(strainyz, strainxzr, strainyz_hilbert, strainxzr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainxz, strainyzr, strainxz_hilbert, strainyzr_hilbert, sgnh, dim=3)
                    grad_c46(interior_region) = grad_c46(interior_region) &
                        - compute_directional_gradient(strainyz, strainxyr, strainyz_hilbert, strainxyr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainxy, strainyzr, strainxy_hilbert, strainyzr_hilbert, sgnh, dim=3)

                    grad_c55(interior_region) = grad_c55(interior_region) &
                        - compute_directional_gradient(strainxz, strainxzr, strainxz_hilbert, strainxzr_hilbert, sgnh, dim=3)
                    grad_c56(interior_region) = grad_c56(interior_region) &
                        - compute_directional_gradient(strainxz, strainxyr, strainxz_hilbert, strainxyr_hilbert, sgnh, dim=3) &
                        - compute_directional_gradient(strainxy, strainxzr, strainxy_hilbert, strainxzr_hilbert, sgnh, dim=3)

                    grad_c66(interior_region) = grad_c66(interior_region) &
                        - compute_directional_gradient(strainxy, strainxyr, strainxy_hilbert, strainxyr_hilbert, sgnh, dim=3)

                end if

            end if

            if (yn_energy_precond) then

                !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
                do k = nz1_interior, nz2_interior
                    do j = ny1_interior, ny2_interior
                        do i = nx1_interior, nx2_interior
                            energy_src_v(i, j, k) = energy_src_v(i, j, k) &
                                + strainxx(i, j, k)**2 + strainyy(i, j, k)**2 + strainzz(i, j, k)**2 &
                                + 2*strainyz(i, j, k)**2 + 2*strainxz(i, j, k)**2 + 2*strainxy(i, j, k)**2
                            energy_rec_v(i, j, k) = energy_rec_v(i, j, k) &
                                + strainxxr(i, j, k)**2 + strainyyr(i, j, k)**2 + strainzzr(i, j, k)**2 &
                                + 2*strainyzr(i, j, k)**2 + 2*strainxzr(i, j, k)**2 + 2*strainxyr(i, j, k)**2
                        end do
                    end do
                end do
                !$omp end parallel do

            end if

        end if

        if (kernel_a /= '') then

            !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
            do k = nz1, nz2
                do j = ny1, ny2
                    do i = nx1, nx2

                        src_vx(i, j, k) = &
                            0.5*sum(vx_hxiyiz(i:i + 1, j, k)) &
                            + 0.5*sum(vx_ixhyiz(i, j:j + 1, k)) &
                            + 0.5*sum(vx_ixiyhz(i, j, k:k + 1)) &
                            + 0.125*sum(vx_hxhyhz(i:i + 1, j:j + 1, k:k + 1)) &
                            - prev_vx(i, j, k)
                        src_vy(i, j, k) = &
                            0.5*sum(vy_hxiyiz(i:i + 1, j, k)) &
                            + 0.5*sum(vy_ixhyiz(i, j:j + 1, k)) &
                            + 0.5*sum(vy_ixiyhz(i, j, k:k + 1)) &
                            + 0.125*sum(vy_hxhyhz(i:i + 1, j:j + 1, k:k + 1)) &
                            - prev_vy(i, j, k)
                        src_vz(i, j, k) = &
                            0.5*sum(vz_hxiyiz(i:i + 1, j, k)) &
                            + 0.5*sum(vz_ixhyiz(i, j:j + 1, k)) &
                            + 0.5*sum(vz_ixiyhz(i, j, k:k + 1)) &
                            + 0.125*sum(vz_hxhyhz(i:i + 1, j:j + 1, k:k + 1)) &
                            - prev_vz(i, j, k)

                        rec_vx(i, j, k) = &
                            0.5*sum(vxr_hxiyiz(i:i + 1, j, k)) &
                            + 0.5*sum(vxr_ixhyiz(i, j:j + 1, k)) &
                            + 0.5*sum(vxr_ixiyhz(i, j, k:k + 1)) &
                            + 0.125*sum(vxr_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                        rec_vy(i, j, k) = &
                            0.5*sum(vyr_hxiyiz(i:i + 1, j, k)) &
                            + 0.5*sum(vyr_ixhyiz(i, j:j + 1, k)) &
                            + 0.5*sum(vyr_ixiyhz(i, j, k:k + 1)) &
                            + 0.125*sum(vyr_hxhyhz(i:i + 1, j:j + 1, k:k + 1))
                        rec_vz(i, j, k) = &
                            0.5*sum(vzr_hxiyiz(i:i + 1, j, k)) &
                            + 0.5*sum(vzr_ixhyiz(i, j:j + 1, k)) &
                            + 0.5*sum(vzr_ixiyhz(i, j, k:k + 1)) &
                            + 0.125*sum(vzr_hxhyhz(i:i + 1, j:j + 1, k:k + 1))

                    end do
                end do
            end do
            !$omp end parallel do

            if (kernel_a == 'full') then

                !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
                do k = nz1_interior, nz2_interior
                    do j = ny1_interior, ny2_interior
                        do i = nx1_interior, nx2_interior
                            grad_rho(i, j, k) = grad_rho(i, j, k) + src_vx(i, j, k)*rec_vx(i, j, k) &
                                + src_vy(i, j, k)*rec_vy(i, j, k) + src_vz(i, j, k)*rec_vz(i, j, k)
                        end do
                    end do
                end do
                !$omp end parallel do

            else

                ! Along x
                if (index(kernel_a, 'lowx') /= 0) then
                    sgnh = 1
                else if (index(kernel_a, 'highx') /= 0) then
                    sgnh = -1
                else
                    sgnh = 0
                end if

                if (sgnh /= 0) then

                    call commute_array(src_vx, htlen, dim=1)
                    call commute_array(rec_vx, htlen, dim=1)
                    src_hilbert = compute_hilbert_transform(src_vx, dim=1)
                    rec_hilbert = compute_hilbert_transform(rec_vx, dim=1)
                    grad_rho(interior_region) = grad_rho(interior_region) &
                        + compute_directional_gradient(src_vx, rec_vx, src_hilbert, rec_hilbert, sgnh, dim=1)

                    call commute_array(src_vy, htlen, dim=1)
                    call commute_array(rec_vy, htlen, dim=1)
                    src_hilbert = compute_hilbert_transform(src_vy, dim=1)
                    rec_hilbert = compute_hilbert_transform(rec_vy, dim=1)
                    grad_rho(interior_region) = grad_rho(interior_region) &
                        + compute_directional_gradient(src_vy, rec_vy, src_hilbert, rec_hilbert, sgnh, dim=1)

                    call commute_array(src_vz, htlen, dim=1)
                    call commute_array(rec_vz, htlen, dim=1)
                    src_hilbert = compute_hilbert_transform(src_vz, dim=1)
                    rec_hilbert = compute_hilbert_transform(rec_vz, dim=1)
                    grad_rho(interior_region) = grad_rho(interior_region) &
                        + compute_directional_gradient(src_vz, rec_vz, src_hilbert, rec_hilbert, sgnh, dim=1)

                end if

                ! Along x
                if (index(kernel_a, 'lowy') /= 0) then
                    sgnh = 1
                else if (index(kernel_a, 'highy') /= 0) then
                    sgnh = -1
                else
                    sgnh = 0
                end if

                if (sgnh /= 0) then

                    call commute_array(src_vx, htlen, dim=2)
                    call commute_array(rec_vx, htlen, dim=2)
                    src_hilbert = compute_hilbert_transform(src_vx, dim=2)
                    rec_hilbert = compute_hilbert_transform(rec_vx, dim=2)
                    grad_rho(interior_region) = grad_rho(interior_region) &
                        + compute_directional_gradient(src_vx, rec_vx, src_hilbert, rec_hilbert, sgnh, dim=2)

                    call commute_array(src_vy, htlen, dim=2)
                    call commute_array(rec_vy, htlen, dim=2)
                    src_hilbert = compute_hilbert_transform(src_vy, dim=2)
                    rec_hilbert = compute_hilbert_transform(rec_vy, dim=2)
                    grad_rho(interior_region) = grad_rho(interior_region) &
                        + compute_directional_gradient(src_vy, rec_vy, src_hilbert, rec_hilbert, sgnh, dim=2)

                    call commute_array(src_vz, htlen, dim=2)
                    call commute_array(rec_vz, htlen, dim=2)
                    src_hilbert = compute_hilbert_transform(src_vz, dim=2)
                    rec_hilbert = compute_hilbert_transform(rec_vz, dim=2)
                    grad_rho(interior_region) = grad_rho(interior_region) &
                        + compute_directional_gradient(src_vz, rec_vz, src_hilbert, rec_hilbert, sgnh, dim=2)

                end if

                ! Along x
                if (index(kernel_a, 'lowz') /= 0) then
                    sgnh = 1
                else if (index(kernel_a, 'highz') /= 0) then
                    sgnh = -1
                else
                    sgnh = 0
                end if

                if (sgnh /= 0) then

                    call commute_array(src_vx, htlen, dim=3)
                    call commute_array(rec_vx, htlen, dim=3)
                    src_hilbert = compute_hilbert_transform(src_vx, dim=3)
                    rec_hilbert = compute_hilbert_transform(rec_vx, dim=3)
                    grad_rho(interior_region) = grad_rho(interior_region) &
                        + compute_directional_gradient(src_vx, rec_vx, src_hilbert, rec_hilbert, sgnh, dim=3)

                    call commute_array(src_vy, htlen, dim=3)
                    call commute_array(rec_vy, htlen, dim=3)
                    src_hilbert = compute_hilbert_transform(src_vy, dim=3)
                    rec_hilbert = compute_hilbert_transform(rec_vy, dim=3)
                    grad_rho(interior_region) = grad_rho(interior_region) &
                        + compute_directional_gradient(src_vy, rec_vy, src_hilbert, rec_hilbert, sgnh, dim=3)

                    call commute_array(src_vz, htlen, dim=3)
                    call commute_array(rec_vz, htlen, dim=3)
                    src_hilbert = compute_hilbert_transform(src_vz, dim=3)
                    rec_hilbert = compute_hilbert_transform(rec_vz, dim=3)
                    grad_rho(interior_region) = grad_rho(interior_region) &
                        + compute_directional_gradient(src_vz, rec_vz, src_hilbert, rec_hilbert, sgnh, dim=3)

                end if

            end if

            if (yn_energy_precond) then

                !$omp parallel do private(i, j, k) collapse(3) schedule(auto)
                do k = nz1_interior, nz2_interior
                    do j = ny1_interior, ny2_interior
                        do i = nx1_interior, nx2_interior
                            energy_src_a(i, j, k) = energy_src_a(i, j, k) + src_vx(i, j, k)**2 + src_vy(i, j, k)**2 + src_vz(i, j, k)**2
                            energy_rec_a(i, j, k) = energy_rec_a(i, j, k) + rec_vx(i, j, k)**2 + rec_vy(i, j, k)**2 + rec_vz(i, j, k)**2
                        end do
                    end do
                end do
                !$omp end parallel do

            end if

        end if

    end subroutine

#include '../../lib/inc_directional_gradient.f90'

    !
    !> Compute source parameter gradients
    !
    subroutine compute_gradient_source(t)

        integer, intent(in) :: t

        integer :: sgx, sgy, sgz, irx, iry, irz
        integer :: i

        do i = 1, sgmtr%ns

            sgx = sgmtr%srcr(i)%gx
            sgy = sgmtr%srcr(i)%gy
            sgz = sgmtr%srcr(i)%gz

            do irz = -nkw, nkw
                do iry = -nkw, nkw
                    do irx = -nkw, nkw
                        if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                            ! ix-iy-iz
                            grad_mt(1) = grad_mt(1) - &
                                stressxxr_ixiyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(2) = grad_mt(2) - &
                                stressyyr_ixiyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(3) = grad_mt(3) - &
                                stresszzr_ixiyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(4) = grad_mt(4) - &
                                stressxyr_ixiyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(5) = grad_mt(5) - &
                                stressxzr_ixiyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(6) = grad_mt(6) - &
                                stressyzr_ixiyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                        end if
                    end do
                end do
            end do

            sgx = sgmtr%srcr(i)%hx
            sgy = sgmtr%srcr(i)%hy
            sgz = sgmtr%srcr(i)%gz
            do irz = -nkw, nkw
                do iry = -nkw, nkw
                    do irx = -nkw, nkw
                        if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                            ! hx-hy-iz
                            grad_mt(1) = grad_mt(1) - &
                                stressxxr_hxhyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(2) = grad_mt(2) - &
                                stressyyr_hxhyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(3) = grad_mt(3) - &
                                stresszzr_hxhyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(4) = grad_mt(4) - &
                                stressxyr_hxhyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(5) = grad_mt(5) - &
                                stressxzr_hxhyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                            grad_mt(6) = grad_mt(6) - &
                                stressyzr_hxhyiz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_iz(irz)*dstf_dt(t, i)
                        end if
                    end do
                end do
            end do

            sgx = sgmtr%srcr(i)%hx
            sgy = sgmtr%srcr(i)%gy
            sgz = sgmtr%srcr(i)%hz
            do irz = -nkw, nkw
                do iry = -nkw, nkw
                    do irx = -nkw, nkw
                        if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                            ! hx-iy-hz
                            grad_mt(1) = grad_mt(1) - &
                                stressxxr_hxiyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                            grad_mt(2) = grad_mt(2) - &
                                stressyyr_hxiyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                            grad_mt(3) = grad_mt(3) - &
                                stresszzr_hxiyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                            grad_mt(4) = grad_mt(4) - &
                                stressxyr_hxiyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                            grad_mt(5) = grad_mt(5) - &
                                stressxzr_hxiyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                            grad_mt(6) = grad_mt(6) - &
                                stressyzr_hxiyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_hx(irx) &
                                *sgmtr%srcr(i)%interp_iy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                        end if
                    end do
                end do
            end do

            sgx = sgmtr%srcr(i)%gx
            sgy = sgmtr%srcr(i)%hy
            sgz = sgmtr%srcr(i)%hz
            do irz = -nkw, nkw
                do iry = -nkw, nkw
                    do irx = -nkw, nkw
                        if (is_in_block(sgx + irx, sgy + iry, sgz + irz) .and. ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                            ! ix-hy-hz
                            grad_mt(1) = grad_mt(1) - &
                                stressxxr_ixhyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                            grad_mt(2) = grad_mt(2) - &
                                stressyyr_ixhyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                            grad_mt(3) = grad_mt(3) - &
                                stresszzr_ixhyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                            grad_mt(4) = grad_mt(4) - &
                                stressxyr_ixhyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                            grad_mt(5) = grad_mt(5) - &
                                stressxzr_ixhyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                            grad_mt(6) = grad_mt(6) - &
                                stressyzr_ixhyhz(sgx + irx, sgy + iry, sgz + irz) &
                                *sgmtr%srcr(i)%interp_ix(irx) &
                                *sgmtr%srcr(i)%interp_hy(iry) &
                                *sgmtr%srcr(i)%interp_hz(irz)*dstf_dt(t, i)
                        end if
                    end do
                end do
            end do

        end do

    end subroutine

end submodule
