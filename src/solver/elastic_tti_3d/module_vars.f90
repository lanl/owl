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


module elastic_tti_3d_vars

    use libflit
    use mod_fdcoef
    use mod_su
    use mod_source_receiver, only: source_receiver_geometry
    use mod_utility, only: check_dt_f0
    use mod_source_receiver, only : nkw
    use mod_anisotropy, only: thomsen_to_cij, alkhalifah_tsvankin_to_cij, min_max_phase_velocity_3d
    use mpi_f08

    use elastic_tti_3d

    implicit none

#include 'libflit.macro'

    ! Forward wavefield
    ! stress -- set a
    real, allocatable, dimension(:, :, :) :: stressxx_ixiyiz, stressyy_ixiyiz, stresszz_ixiyiz
    real, allocatable, dimension(:, :, :) :: stressxy_ixiyiz, stressxz_ixiyiz, stressyz_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvx_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvy_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvz_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvx_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvy_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvz_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvx_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvy_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvz_ixiyiz

    ! stress -- set b
    real, allocatable, dimension(:, :, :) :: stressxx_hxhyiz, stressyy_hxhyiz, stresszz_hxhyiz
    real, allocatable, dimension(:, :, :) :: stressxy_hxhyiz, stressxz_hxhyiz, stressyz_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvx_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvy_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvz_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvx_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvy_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvz_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvx_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvy_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvz_hxhyiz

    ! stress -- set c
    real, allocatable, dimension(:, :, :) :: stressxx_hxiyhz, stressyy_hxiyhz, stresszz_hxiyhz
    real, allocatable, dimension(:, :, :) :: stressxy_hxiyhz, stressxz_hxiyhz, stressyz_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvx_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvy_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvz_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvx_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvy_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvz_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvx_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvy_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvz_hxiyhz

    ! stress -- set d
    real, allocatable, dimension(:, :, :) :: stressxx_ixhyhz, stressyy_ixhyhz, stresszz_ixhyhz
    real, allocatable, dimension(:, :, :) :: stressxy_ixhyhz, stressxz_ixhyhz, stressyz_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvx_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvy_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvz_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvx_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvy_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvz_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvx_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvy_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvz_ixhyhz

    ! velocity -- set a
    real, allocatable, dimension(:, :, :) :: vx_hxiyiz, vy_hxiyiz, vz_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxx_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxy_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxz_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyxy_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyyy_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyyz_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzxx_hxiyiz, memory_pdzxy_hxiyiz, memory_pdzxz_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzyy_hxiyiz, memory_pdzyz_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzzz_hxiyiz

    ! velocity -- set b
    real, allocatable, dimension(:, :, :) :: vx_ixhyiz, vy_ixhyiz, vz_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxx_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxy_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxz_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyxy_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyyy_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyyz_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzxx_ixhyiz, memory_pdzxy_ixhyiz, memory_pdzxz_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzyy_ixhyiz, memory_pdzyz_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzzz_ixhyiz

    ! velocity -- set c
    real, allocatable, dimension(:, :, :) :: vx_ixiyhz, vy_ixiyhz, vz_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxx_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxy_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxz_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyxy_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyyy_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyyz_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzxx_ixiyhz, memory_pdzxy_ixiyhz, memory_pdzxz_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzyy_ixiyhz, memory_pdzyz_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzzz_ixiyhz

    ! velocity -- set d
    real, allocatable, dimension(:, :, :) :: vx_hxhyhz, vy_hxhyhz, vz_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxx_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxy_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxz_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyxy_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyyy_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyyz_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzxx_hxhyhz, memory_pdzxy_hxhyhz, memory_pdzxz_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzyy_hxhyhz, memory_pdzyz_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzzz_hxhyhz

    ! Adjoint wavefields
    ! stress -- set a
    real, allocatable, dimension(:, :, :) :: stressxxr_ixiyiz, stressyyr_ixiyiz, stresszzr_ixiyiz
    real, allocatable, dimension(:, :, :) :: stressxyr_ixiyiz, stressxzr_ixiyiz, stressyzr_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvxr_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvyr_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvzr_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvxr_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvyr_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvzr_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvxr_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvyr_ixiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvzr_ixiyiz

    ! stress -- set b
    real, allocatable, dimension(:, :, :) :: stressxxr_hxhyiz, stressyyr_hxhyiz, stresszzr_hxhyiz
    real, allocatable, dimension(:, :, :) :: stressxyr_hxhyiz, stressxzr_hxhyiz, stressyzr_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvxr_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvyr_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxvzr_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvxr_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvyr_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyvzr_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvxr_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvyr_hxhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzvzr_hxhyiz

    ! stress -- set c
    real, allocatable, dimension(:, :, :) :: stressxxr_hxiyhz, stressyyr_hxiyhz, stresszzr_hxiyhz
    real, allocatable, dimension(:, :, :) :: stressxyr_hxiyhz, stressxzr_hxiyhz, stressyzr_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvxr_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvyr_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvzr_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvxr_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvyr_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvzr_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvxr_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvyr_hxiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvzr_hxiyhz

    ! stress -- set d
    real, allocatable, dimension(:, :, :) :: stressxxr_ixhyhz, stressyyr_ixhyhz, stresszzr_ixhyhz
    real, allocatable, dimension(:, :, :) :: stressxyr_ixhyhz, stressxzr_ixhyhz, stressyzr_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvxr_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvyr_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxvzr_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvxr_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvyr_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyvzr_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvxr_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvyr_ixhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzvzr_ixhyhz

    ! velocity -- set a
    real, allocatable, dimension(:, :, :) :: vxr_hxiyiz, vyr_hxiyiz, vzr_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxxr_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxyr_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxzr_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyxyr_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyyyr_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyyzr_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzxxr_hxiyiz, memory_pdzxyr_hxiyiz, memory_pdzxzr_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzyyr_hxiyiz, memory_pdzyzr_hxiyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzzzr_hxiyiz

    ! velocity -- set b
    real, allocatable, dimension(:, :, :) :: vxr_ixhyiz, vyr_ixhyiz, vzr_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxxr_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxyr_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdxxzr_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyxyr_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyyyr_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdyyzr_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzxxr_ixhyiz, memory_pdzxyr_ixhyiz, memory_pdzxzr_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzyyr_ixhyiz, memory_pdzyzr_ixhyiz
    real, allocatable, dimension(:, :, :) :: memory_pdzzzr_ixhyiz

    ! velocity -- set c
    real, allocatable, dimension(:, :, :) :: vxr_ixiyhz, vyr_ixiyhz, vzr_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxxr_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxyr_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxzr_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyxyr_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyyyr_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyyzr_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzxxr_ixiyhz, memory_pdzxyr_ixiyhz, memory_pdzxzr_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzyyr_ixiyhz, memory_pdzyzr_ixiyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzzzr_ixiyhz

    ! velocity -- set d
    real, allocatable, dimension(:, :, :) :: vxr_hxhyhz, vyr_hxhyhz, vzr_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxxr_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxyr_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdxxzr_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyxyr_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyyyr_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdyyzr_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzxxr_hxhyhz, memory_pdzxyr_hxhyhz, memory_pdzxzr_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzyyr_hxhyhz, memory_pdzyzr_hxhyhz
    real, allocatable, dimension(:, :, :) :: memory_pdzzzr_hxhyhz

    real, allocatable, dimension(:, :, :) :: energy_src_v, energy_rec_v
    real, allocatable, dimension(:, :, :) :: energy_src_a, energy_rec_a

    real, allocatable, dimension(:, :, :) :: snapvx, snapvy, snapvz
    real, allocatable, dimension(:, :, :) :: snapvxr, snapvyr, snapvzr

    real, allocatable, dimension(:, :, :) :: prev_vx, prev_vy, prev_vz
    real, allocatable, dimension(:, :, :) :: prev_stressxx, prev_stressyy, prev_stresszz
    real, allocatable, dimension(:, :, :) :: prev_stressyz, prev_stressxz, prev_stressxy

    real, allocatable, dimension(:, :, :) :: prev_pr, prev_stressxxr, prev_stressyyr, prev_stresszzr
    real, allocatable, dimension(:, :, :) :: prev_stressyzr, prev_stressxzr, prev_stressxyr

    real, allocatable, dimension(:, :, :) :: strainxx, strainyy, strainzz, strainxz, strainyz, strainxy
    real, allocatable, dimension(:, :, :) :: strainxxr, strainyyr, strainzzr, strainxzr, strainyzr, strainxyr
    real, allocatable, dimension(:, :, :) :: src_vx, src_vy, src_vz, rec_vx, rec_vy, rec_vz

    real, allocatable, dimension(:, :, :) :: strainxx_hilbert, strainyy_hilbert, strainzz_hilbert, strainxz_hilbert, strainyz_hilbert, strainxy_hilbert
    real, allocatable, dimension(:, :, :) :: strainxxr_hilbert, strainyyr_hilbert, strainzzr_hilbert, strainxzr_hilbert, strainyzr_hilbert, strainxyr_hilbert
    real, allocatable, dimension(:, :, :) :: src_hilbert, rec_hilbert

    real, allocatable, dimension(:, :, :) :: vp, vs, rho, tieps, tidel, tigam, tithe, tiphi, tieta

    real, allocatable, dimension(:, :, :) :: c11, c12, c13, c14, c15, c16, &
        c22, c23, c24, c25, c26, &
        c33, c34, c35, c36, &
        c44, c45, c46, &
        c55, c56, &
        c66

    real, allocatable, dimension(:, :, :) :: s11, s12, s13, s14, s15, s16, &
        s22, s23, s24, s25, s26, &
        s33, s34, s35, s36, &
        s44, s45, s46, &
        s55, s56, &
        s66

    real, allocatable, dimension(:, :, :) :: grad_c11, grad_c12, grad_c13, grad_c14, grad_c15, grad_c16, &
        grad_c22, grad_c23, grad_c24, grad_c25, grad_c26, &
        grad_c33, grad_c34, grad_c35, grad_c36, &
        grad_c44, grad_c45, grad_c46, &
        grad_c55, grad_c56, &
        grad_c66

    real, allocatable, dimension(:, :, :) :: grad_vp, grad_vs, grad_eps, grad_del, grad_gam, grad_rho, grad_eta

    real, allocatable, dimension(:, :, :) :: energy_src, energy_rec

    real, allocatable, dimension(:) :: snaps

    integer :: np

    integer :: nx, ny, nz
    real :: dx, dy, dz
    real :: ox, oy, oz

    integer :: nt
    real :: dt
    real :: tmax

    integer :: data_nt
    real :: data_dt
    real :: data_tmax

    integer :: pml

    logical :: yn_reconstruct = .false.
    logical :: yn_free_surface = .false.
    type(source_receiver_geometry) :: sgmtr

    integer :: cc_step_interval
    logical :: verbose = .false.

    logical :: yn_energy_precond = .true.

    character(len=1024) :: dir_synthetic, dir_snapshot, dir_working, dir_adjoint

    logical :: yn_compx = .true.
    logical :: yn_compy = .true.
    logical :: yn_compz = .true.

    character(len=12) :: aniso_param = 'iso'

    integer :: nx1, nx2, ny1, ny2, nz1, nz2
    integer :: nx1_interior, nx2_interior, ny1_interior, ny2_interior, nz1_interior, nz2_interior

    real :: pmlvp

    real, allocatable, dimension(:, :, :) :: zz_i, zz_h, dz_i, dz_h
    real, allocatable, dimension(:, :, :) :: dz_scaling_i, dz_scaling_h

    real, allocatable, dimension(:, :) :: topo_ixiy, topo_ixhy, topo_hxiy, topo_hxhy
    real, allocatable, dimension(:, :) :: slopex_ixiy, slopex_ixhy, slopex_hxiy, slopex_hxhy
    real, allocatable, dimension(:, :) :: slopey_ixiy, slopey_ixhy, slopey_hxiy, slopey_hxhy
    type(meta_array2_real), dimension(:, :), allocatable :: coefiix, coefiiy, coefiiz, coefiizr
    type(meta_array2_real), dimension(:, :), allocatable :: coefhhx, coefhhy, coefhhz, coefhhzr

    real :: depth_max, eta_max

    real, allocatable, dimension(:) :: eta_zz_i, eta_zz_h
    real, allocatable, dimension(:) :: eta_dz_i, eta_dz_h
    real, allocatable, dimension(:) :: eta_dz_scaling_i, eta_dz_scaling_h

    real, allocatable, dimension(:) :: topox, topoy
    real, allocatable, dimension(:, :) :: topoz

    integer :: htlen = 30
    real :: topo_max
    integer, allocatable, dimension(:, :) :: trace_range

    integer :: nc_mt
    real, allocatable, dimension(:) :: grad_mt
    logical :: yn_update_medium = .true.
    logical :: yn_update_source = .false.
    real, allocatable, dimension(:, :) :: dstf_dt

contains

    !
    !> Map the model on the original regular mesh to the depth-varying mesh
    !
    subroutine map_regular_to_irregular(v, this, rr)

        real, allocatable, dimension(:, :, :), intent(inout) :: v
        type(wave_solver_elastic_tti_3d), intent(in) :: this
        integer, dimension(1:6), intent(in) :: rr

        real, allocatable, dimension(:, :, :) :: m
        real, allocatable, dimension(:) :: z
        integer, allocatable, dimension(:) :: r
        integer :: i, j

        r = zeros(6)
        r(1) = lbound(v, dim=1)
        r(2) = ubound(v, dim=1)
        r(3) = lbound(v, dim=2)
        r(4) = ubound(v, dim=2)
        r(5) = lbound(v, dim=3)
        r(6) = ubound(v, dim=3)

        ! Depth coordinates of regular-mesh array
        z = regspace(r(5) - 1.0, 1.0, r(6) - 1.0)*this%dz - maxval(topo_ixiy)

        ! Interpolate to irregular-mesh array
        call alloc_array(m, rr)

        !$omp parallel do private(i, j)
        do j = max(r(3), rr(3)), min(r(4), rr(4))
            do i = max(r(1), rr(1)), min(r(2), rr(2))

                m(i, j, :) = ginterp(z, v(i, j, r(5):r(6)), zz_i(i, j, :), 'linear')

            end do
        end do
        !$omp end parallel do

        ! Overwrite original
        v = m

    end subroutine

    !
    !> Map the model on the depth-varying mesh to the original regular mesh
    !
    subroutine map_irregular_to_regular(v, this, rr)

        real, allocatable, dimension(:, :, :), intent(inout) :: v
        type(wave_solver_elastic_tti_3d), intent(in) :: this
        integer, dimension(1:6), intent(in) :: rr

        real, allocatable, dimension(:, :, :) :: m
        real, allocatable, dimension(:) :: zz
        integer, allocatable, dimension(:) :: r
        integer :: i, j, k
        real :: max_topo

        r = zeros(6)
        r(1) = lbound(v, dim=1)
        r(2) = ubound(v, dim=1)
        r(3) = lbound(v, dim=2)
        r(4) = ubound(v, dim=2)
        r(5) = lbound(v, dim=3)
        r(6) = ubound(v, dim=3)

        zz = regspace(rr(5) - 1.0, 1.0, rr(6) - 1.0)*this%dz

        call alloc_array(m, rr)

        max_topo = maxval(topo_ixiy)

        !$omp parallel do private(i, j, k)
        do j = max(r(3), rr(3)), min(r(4), rr(4))
            do i = max(r(1), rr(1)), min(r(2), rr(2))

                m(i, j, :) = ginterp(zz_i(i, j, r(5):r(6)), v(i, j, r(5):r(6)), zz - max_topo, 'cubic')

                ! Mask out points above the free surface
                do k = 1, size(zz)
                    if (zz(k) < max_topo - topo_ixiy(i, j)) then
                        m(i, j, k) = 0
                    end if
                end do

            end do
        end do
        !$omp end parallel do

        v = m

    end subroutine

    !
    !> Prepare for modeling
    !
    subroutine prepare_modeling(this)

        type(wave_solver_elastic_tti_3d), intent(in) :: this

        real :: minv, maxv, dtstable, f0clean, temp
        integer :: i, j, k, l, refine_nz, ntp
        real :: rr, dz_max, rayleigh_wavelength
        real :: tmpmin, tmpmax
        real, allocatable, dimension(:, :, :) :: mdispersion, mstability
        real, allocatable, dimension(:, :) :: topo
        real, allocatable, dimension(:) :: stp, rtp
        real :: dz0
        real :: wmin, wmax
        real :: stable_eta_max, rayleigh_vs
        integer :: nbeg, nend

        nx = this%nx
        ny = this%ny
        nz = this%nz
        dx = this%dx
        dy = this%dy
        dz = this%dz
        ox = this%ox
        oy = this%oy
        oz = this%oz

        pml = this%pml

        ! Allocate memory for Cij models
        c11 = zeros(nx, ny, nz)
        c12 = zeros(nx, ny, nz)
        c13 = zeros(nx, ny, nz)
        c14 = zeros(nx, ny, nz)
        c15 = zeros(nx, ny, nz)
        c16 = zeros(nx, ny, nz)
        c22 = zeros(nx, ny, nz)
        c23 = zeros(nx, ny, nz)
        c24 = zeros(nx, ny, nz)
        c25 = zeros(nx, ny, nz)
        c26 = zeros(nx, ny, nz)
        c33 = zeros(nx, ny, nz)
        c34 = zeros(nx, ny, nz)
        c35 = zeros(nx, ny, nz)
        c36 = zeros(nx, ny, nz)
        c44 = zeros(nx, ny, nz)
        c45 = zeros(nx, ny, nz)
        c46 = zeros(nx, ny, nz)
        c55 = zeros(nx, ny, nz)
        c56 = zeros(nx, ny, nz)
        c66 = zeros(nx, ny, nz)
        rho = zeros(nx, ny, nz)

        aniso_param = this%anisotropy_type
        select case (aniso_param)

            case ('iso')
                vp = permute(this%vp, 321)
                vs = permute(this%vs, 321)
                rho = permute(this%rho, 321)

                call assert(maxval(abs(vp)) > 0, ' Error: vp is all zero.')
                call assert(maxval(abs(vs)) > 0, ' Error: vs is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

                c11 = rho*vp**2
                c22 = c11
                c33 = c11
                c44 = rho*vs**2
                c55 = c44
                c66 = c44
                c12 = c11 - 2*c55
                c13 = c12
                c23 = c12

                minv = minval(vs)
                maxv = maxval(vp)

            case ('thomsen')
                vp = permute(this%vp, 321)
                vs = permute(this%vs, 321)
                tieps = permute(this%tieps, 321)
                tidel = permute(this%tidel, 321)
                tigam = permute(this%tigam, 321)
                tithe = permute(this%tithe, 321)
                tiphi = permute(this%tiphi, 321)
                rho = permute(this%rho, 321)

                call assert(maxval(abs(vp)) > 0, ' Error: vp is all zero.')
                call assert(maxval(abs(vs)) > 0, ' Error: vs is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

                !$omp parallel do private(i, j, k) collapse(3)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            call thomsen_to_cij(vp(i, j, k), vs(i, j, k), rho(i, j, k), &
                                tieps(i, j, k), tidel(i, j, k), tigam(i, j, k), tithe(i, j, k), tiphi(i, j, k), &
                                c11(i, j, k), c12(i, j, k), c13(i, j, k), c14(i, j, k), c15(i, j, k), c16(i, j, k), &
                                c22(i, j, k), c23(i, j, k), c24(i, j, k), c25(i, j, k), c26(i, j, k), &
                                c33(i, j, k), c34(i, j, k), c35(i, j, k), c36(i, j, k), &
                                c44(i, j, k), c45(i, j, k), c46(i, j, k), &
                                c55(i, j, k), c56(i, j, k), &
                                c66(i, j, k))
                        end do
                    end do
                end do
                !$omp end parallel do

            case ('a-t')
                vp = permute(this%vp, 321)
                vs = permute(this%vs, 321)
                tieps = permute(this%tieps, 321)
                tieta = permute(this%tieta, 321)
                tigam = permute(this%tigam, 321)
                tithe = permute(this%tithe, 321)
                tiphi = permute(this%tiphi, 321)
                rho = permute(this%rho, 321)

                call assert(maxval(abs(vp)) > 0, ' Error: vp is all zero.')
                call assert(maxval(abs(vs)) > 0, ' Error: vs is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

                !$omp parallel do private(i, j, k) collapse(3)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            call alkhalifah_tsvankin_to_cij(vp(i, j, k), vs(i, j, k), rho(i, j, k), &
                                tieps(i, j, k), tieta(i, j, k), tigam(i, j, k), tithe(i, j, k), tiphi(i, j, k), &
                                c11(i, j, k), c12(i, j, k), c13(i, j, k), c14(i, j, k), c15(i, j, k), c16(i, j, k), &
                                c22(i, j, k), c23(i, j, k), c24(i, j, k), c25(i, j, k), c26(i, j, k), &
                                c33(i, j, k), c34(i, j, k), c35(i, j, k), c36(i, j, k), &
                                c44(i, j, k), c45(i, j, k), c46(i, j, k), &
                                c55(i, j, k), c56(i, j, k), &
                                c66(i, j, k))
                        end do
                    end do
                end do
                !$omp end parallel do

            case ('cij')
                c11 = permute(this%c11, 321)
                c12 = permute(this%c12, 321)
                c13 = permute(this%c13, 321)
                c14 = permute(this%c14, 321)
                c15 = permute(this%c15, 321)
                c16 = permute(this%c16, 321)
                c22 = permute(this%c22, 321)
                c23 = permute(this%c23, 321)
                c24 = permute(this%c24, 321)
                c25 = permute(this%c25, 321)
                c26 = permute(this%c26, 321)
                c33 = permute(this%c33, 321)
                c34 = permute(this%c34, 321)
                c35 = permute(this%c35, 321)
                c36 = permute(this%c36, 321)
                c44 = permute(this%c44, 321)
                c45 = permute(this%c45, 321)
                c46 = permute(this%c46, 321)
                c55 = permute(this%c55, 321)
                c56 = permute(this%c56, 321)
                c66 = permute(this%c66, 321)
                rho = permute(this%rho, 321)

                call assert(maxval(abs(c11)) > 0, ' Error: c11 is all zero.')
                call assert(maxval(abs(c22)) > 0, ' Error: c22 is all zero.')
                call assert(maxval(abs(c33)) > 0, ' Error: c33 is all zero.')
                call assert(maxval(abs(c44)) > 0, ' Error: c44 is all zero.')
                call assert(maxval(abs(c55)) > 0, ' Error: c55 is all zero.')
                call assert(maxval(abs(c66)) > 0, ' Error: c66 is all zero.')
                call assert(maxval(abs(rho)) > 0, ' Error: rho is all zero.')

        end select

        ! Get PML Vp at boundary nodes
        pmlvp = -float_large
        rayleigh_vs = 0.0
        select case (aniso_param)

            case ('iso')

                !$omp parallel do private(i, j, k) collapse(3) reduction(max:pmlvp)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            if (i == 1 .or. i == nx .or. j == 1 .or. j == ny .or. k == 1 .or. k == nz) then
                                pmlvp = max(pmlvp, vp(i, j, k))
                            end if
                            if (k == 1) then
                                rayleigh_vs = rayleigh_vs + vs(i, j, k)
                            end if
                        end do
                    end do
                end do
                !$omp end parallel do

            case ('thomsen')

                !$omp parallel do private(i, j, k) collapse(3) reduction(max:pmlvp)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            if (i == 1 .or. i == nx .or. j == 1 .or. j == ny .or. k == 1 .or. k == nz) then
                                pmlvp = max(pmlvp, vp(i, j, k)*sqrt(1.0 + 2*tieps(i, j, k)))
                            end if
                            if (k == 1) then
                                rayleigh_vs = rayleigh_vs + vs(i, j, k)
                            end if
                        end do
                    end do
                end do
                !$omp end parallel do

            case ('a-t')

                !$omp parallel do private(i, j, k) collapse(3) reduction(max:pmlvp)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            if (i == 1 .or. i == nx .or. j == 1 .or. j == ny .or. k == 1 .or. k == nz) then
                                pmlvp = max(pmlvp, vp(i, j, k))
                            end if
                            if (k == 1) then
                                rayleigh_vs = rayleigh_vs + vs(i, j, k)
                            end if
                        end do
                    end do
                end do
                !$omp end parallel do

            case ('cij')

                !$omp parallel do private(i, j, k) collapse(3) reduction(max:pmlvp)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            if (i == 1 .or. i == nx .or. j == 1 .or. j == ny .or. k == 1 .or. k == nz) then
                                call min_max_phase_velocity_3d(c11(i, j, k), c12(i, j, k), c13(i, j, k), c14(i, j, k), c15(i, j, k), c16(i, j, k), &
                                    c22(i, j, k), c23(i, j, k), c24(i, j, k), c25(i, j, k), c26(i, j, k), &
                                    c33(i, j, k), c34(i, j, k), c35(i, j, k), c36(i, j, k), &
                                    c44(i, j, k), c45(i, j, k), c46(i, j, k), &
                                    c55(i, j, k), c56(i, j, k), &
                                    c66(i, j, k), &
                                    rho(i, j, k), 36, 36, tmpmax, tmpmin)
                                pmlvp = max(pmlvp, tmpmax)
                                if (k == 1) then
                                    rayleigh_vs = rayleigh_vs + tmpmin
                                end if
                            end if
                        end do
                    end do
                end do
                !$omp end parallel do

        end select

        rayleigh_vs = rayleigh_vs/(nx*ny)

        ! Pad models with PML
        call pad_array(c11, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c12, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c13, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c14, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c15, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c16, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c22, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c23, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c24, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c25, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c26, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c33, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c34, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c35, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c36, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c44, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c45, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c46, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c55, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c56, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(c66, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(rho, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])

        ! Refine near-surface mesh along z when doing free-surface modeling
        yn_free_surface = this%free_surface
        if (yn_free_surface) then

            ! =================================================================
            ! Setup topography
            call alloc_array(topo_ixiy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(topo_ixhy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(topo_hxiy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(topo_hxhy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(slopex_ixiy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(slopex_ixhy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(slopex_hxiy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(slopex_hxhy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(slopey_ixiy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(slopey_ixhy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(slopey_hxiy, [1, nx, 1, ny], pad=pml + 1)
            call alloc_array(slopey_hxhy, [1, nx, 1, ny], pad=pml + 1)

            if (this%file_topo /= '') then

                ntp = count_nonempty_lines(this%file_topo)
                topo = load(this%file_topo, ntp, 3, ascii=.true.)

                topox = meshgrid([nx + 2*pml + 2, ny + 2*pml + 2]*2 + [1, 1], [dx, dy]/2.0, [-(pml + 1)*dx + ox, -(pml + 1)*dy + oy], dim=1)
                topoy = meshgrid([nx + 2*pml + 2, ny + 2*pml + 2]*2 + [1, 1], [dx, dy]/2.0, [-(pml + 1)*dx + ox, -(pml + 1)*dy + oy], dim=2)

                ! Interpolate given xy-elev to obtain topography
                topoz = reshape(ginterp(topo(:, 1), topo(:, 2), topo(:, 3), topox, topoy), [nx + 2*pml + 2, ny + 2*pml + 2]*2 + [1, 1])

                ! Get integer and half spatial points topography by selecting
                topo_ixiy = topoz(2:size(topoz, 1):2, 2:size(topoz, 2):2)
                topo_hxiy = topoz(1:size(topoz, 1) - 1:2, 2:size(topoz, 2):2)
                topo_ixhy = topoz(2:size(topoz, 1):2, 1:size(topoz, 2) - 1:2)
                topo_hxhy = topoz(1:size(topoz, 1) - 1:2, 1:size(topoz, 2) - 1:2)

                ! Compute slopes in x and y directions
                slopex_ixiy = deriv(topo_ixiy, dim=1)/dx
                slopex_hxiy = deriv(topo_hxiy, dim=1)/dx
                slopex_ixhy = deriv(topo_ixhy, dim=1)/dx
                slopex_hxhy = deriv(topo_hxhy, dim=1)/dx
                slopey_ixiy = deriv(topo_ixiy, dim=2)/dy
                slopey_hxiy = deriv(topo_hxiy, dim=2)/dy
                slopey_ixhy = deriv(topo_ixhy, dim=2)/dy
                slopey_hxhy = deriv(topo_hxhy, dim=2)/dy

                if (rankid_group == 0) then
                    call warn(date_time_compact()//' Number of topography control points = '//num2str(ntp))
                end if

                wmax = max(maxval(topo_ixiy), maxval(topo_hxiy), maxval(topo_ixhy), maxval(topo_hxhy))
                if (rankid_group == 0) then
                    call warn(date_time_compact()//' Max elevation above z=0 = '//num2str(wmax, '(es)'))
                end if

                wmin = min(minval(slopex_ixiy), minval(slopex_hxiy), minval(slopex_ixhy), minval(slopex_hxhy))
                wmax = max(maxval(slopex_ixiy), maxval(slopex_hxiy), maxval(slopex_ixhy), maxval(slopex_hxhy))
                if (rankid_group == 0) then
                    call warn(date_time_compact()//' Slope range in x = '//num2str(wmin, '(es)') &
                        //' ~ '//num2str(wmax, '(es)'))
                end if

                wmin = min(minval(slopey_ixiy), minval(slopey_hxiy), minval(slopey_ixhy), minval(slopey_hxhy))
                wmax = max(maxval(slopey_ixiy), maxval(slopey_hxiy), maxval(slopey_ixhy), maxval(slopey_hxhy))
                if (rankid_group == 0) then
                    call warn(date_time_compact()//' Slope range in y = '//num2str(wmin, '(es)') &
                        //' ~ '//num2str(wmax, '(es)'))
                end if

            end if

            ! =================================================================
            ! Setup mesh

            ! Maximum depth in physical and virtual domains
            topo_max = maxval(topo_ixiy)
            depth_max = (this%nz - 1)*this%dz - topo_max
            eta_max = depth_max

            ! Check stability of mesh
            stable_eta_max = min( &
                minval((topo_ixiy + depth_max)/max(1.0, abs(slopex_ixiy), abs(slopey_ixiy))), &
                minval((topo_ixhy + depth_max)/max(1.0, abs(slopex_ixhy), abs(slopey_ixhy))), &
                minval((topo_hxiy + depth_max)/max(1.0, abs(slopex_hxiy), abs(slopey_hxiy))), &
                minval((topo_hxhy + depth_max)/max(1.0, abs(slopex_hxhy), abs(slopey_hxhy))))

            if (eta_max > stable_eta_max) then
                eta_max = 0.975*stable_eta_max
                if (rankid_group == 0) then
                    call warn(date_time_compact()//' Set eta_max = '//num2str(eta_max, '(es)')//' to ensure stability. ')
                end if
            end if

            dz = eta_max/(nz - 1)

            ! Approximate Rayleigh wavelength at the free surface
            rayleigh_wavelength = rayleigh_vs/minval(this%gmtr%srcr(:)%f0)
            refine_nz = ceiling(0.5*rayleigh_wavelength/dz)

            ! Compute integer grid point positions
            ! Initial dz = dz/dz_refine, and it gradually increases after refine_nz, by 2.5% per grid
            ! 5% can cause very weak artificial reflections, while 1% is just too costly.
            rr = this%free_surface_dz_refine
            dz_max = min(this%dz_max*eta_max/(depth_max + topo_max), dz/rr)

            ! First compute a standard depth varying 1D mesh
            eta_zz_i = regspace(0.0, dz_max, (refine_nz - 1)*dz)
            l = size(eta_zz_i)
            do while (eta_zz_i(l) < eta_max)
                dz_max = min(dz_max*1.025, this%dz_max*eta_max/(depth_max + topo_max))
                eta_zz_i = [eta_zz_i, eta_zz_i(l) + dz_max]
                l = l + 1
            end do

            ! New nz is the length of eta_zz_i
            nz = size(eta_zz_i)

            ! Half integer grid points are interpolated
            eta_zz_h = zeros(nz)
            eta_zz_h(1) = -0.5*dz/rr
            do i = 1, nz - 1
                eta_zz_h(i + 1) = 0.5*(eta_zz_i(i) + eta_zz_i(i + 1))
            end do

            call pad_array(eta_zz_i, [pml + 1, pml + 1])
            call pad_array(eta_zz_h, [pml + 1, pml + 1])

            do i = 1, pml + 1
                eta_zz_i(1 - i) = eta_zz_i(1) - i*dz/rr
                eta_zz_i(nz + i) = eta_zz_i(nz) + i*dz_max
                eta_zz_h(1 - i) = eta_zz_h(1) - i*dz/rr
                eta_zz_h(nz + i) = eta_zz_h(nz) + i*dz_max
            end do

            call alloc_array(eta_dz_i, [1, nz], pad=pml + 1, source=deriv(eta_zz_i, method='center'))
            call alloc_array(eta_dz_h, [1, nz], pad=pml + 1, source=deriv(eta_zz_h, method='center'))

            call alloc_array(eta_dz_scaling_i, [1, nz], pad=pml + 1, source=dz/deriv(eta_zz_i, method='center'))
            call alloc_array(eta_dz_scaling_h, [1, nz], pad=pml + 1, source=dz/deriv(eta_zz_h, method='center'))

            ! Then scale the standard irregular mesh to the entire domain
            call alloc_array(zz_i, [1, nx, 1, ny, 1, nz], pad=pml + 1)
            call alloc_array(zz_h, [1, nx, 1, ny, 1, nz], pad=pml + 1)

            !$omp parallel do private(i, j, k) collapse(2)
            do j = -pml, ny + pml + 1
                do i = -pml, nx + pml + 1

                    zz_i(i, j, 1:nz) = rescale(eta_zz_i(1:nz), [-topo_ixiy(i, j), depth_max])

                end do
            end do
            !$omp end parallel do

            dz_max = mean(zz_i(:, :, nz) - zz_i(:, :, nz - 1))

            !$omp parallel do private(i, j, k) collapse(2)
            do j = -pml, ny + pml + 1
                do i = -pml, nx + pml + 1

                    do k = 1, pml + 1
                        zz_i(i, j, nz + k) = zz_i(i, j, nz) + k*dz_max
                    end do

                    zz_h(i, j, 1) = -topo_ixiy(i, j) - 0.5*(zz_i(i, j, 2) - zz_i(i, j, 1))
                    do k = 1, nz + pml
                        zz_h(i, j, k + 1) = 0.5*(zz_i(i, j, k) + zz_i(i, j, k + 1))
                    end do

                    do k = 1, pml + 1
                        zz_i(i, j, 1 - k) = zz_i(i, j, 1) - k*(zz_i(i, j, 2) - zz_i(i, j, 1))
                        zz_h(i, j, 1 - k) = zz_h(i, j, 1) - k*(zz_h(i, j, 2) - zz_h(i, j, 1))
                    end do

                end do
            end do
            !$omp end parallel do

            ! Depth varying grid size
            call alloc_array(dz_i, [1, nx, 1, ny, 1, nz], pad=pml + 1, source=deriv(zz_i, method='center', dim=3))
            call alloc_array(dz_h, [1, nx, 1, ny, 1, nz], pad=pml + 1, source=deriv(zz_h, method='center', dim=3))

            ! Interpolate medium parameter models
            call map_regular_to_irregular(c11, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c12, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c13, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c14, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c15, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c16, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c22, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c23, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c24, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c25, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c26, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c33, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c34, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c35, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c36, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c44, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c45, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c46, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c55, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c56, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(c66, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])
            call map_regular_to_irregular(rho, this, [-pml, nx + pml + 1, -pml, ny + pml + 1, -pml, nz + pml + 1])

            if (rankid_group == 0) then
                call warn(date_time_compact()//' Original nz = '//num2str(this%nz)//', adjusted nz = '//num2str(nz))
                call warn(date_time_compact()//' min(dz) = '//num2str(min(minval(dz_i), minval(dz_h)), '(es)') &
                    //', max(dz) = '//num2str(max(maxval(dz_i), maxval(dz_h)), '(es)'))
            end if

            ! Clip virtual mesh to make sure it does not exceed eta_max
            eta_zz_i = clip(eta_zz_i, -pml*dz, eta_max)
            eta_zz_h = clip(eta_zz_h, -pml*dz, eta_max)

        end if

        call mpibarrier_group

        ! Domain decomposition
        if (yn_free_surface) then

            call domain_decomp_regular_group(nx + 2*pml, ny + 2*pml, nz + pml, nx1, nx2, ny1, ny2, nz1, nz2, &
                weights1=[ones(pml)*1.25, ones(nx), ones(pml)*1.25], &
                weights2=[ones(pml)*1.25, ones(ny), ones(pml)*1.25], &
                weights3=[ones(nz), ones(pml)*1.25])

            nx1 = nx1 - pml
            nx2 = nx2 - pml
            ny1 = ny1 - pml
            ny2 = ny2 - pml

        else

            call domain_decomp_regular_group(nx + 2*pml, ny + 2*pml, nz + 2*pml, nx1, nx2, ny1, ny2, nz1, nz2, &
                weights1=[ones(pml)*1.25, ones(nx), ones(pml)*1.25], &
                weights2=[ones(pml)*1.25, ones(ny), ones(pml)*1.25], &
                weights3=[ones(pml)*1.25, ones(nz), ones(pml)*1.25])

            nx1 = nx1 - pml
            nx2 = nx2 - pml
            ny1 = ny1 - pml
            ny2 = ny2 - pml
            nz1 = nz1 - pml
            nz2 = nz2 - pml

        end if

        nx1_interior = max(1, nx1)
        nx2_interior = min(nx2, nx)
        ny1_interior = max(1, ny1)
        ny2_interior = min(ny2, ny)
        nz1_interior = max(1, nz1)
        nz2_interior = min(nz2, nz)

        ! Check stability and dispersion
        dt = this%dt
        tmax = this%tmax
        temp = sum(abs(fdcoefs))

        if (yn_free_surface) then
            ! Irregular mesh

            call alloc_array(mdispersion, [1, nx, 1, ny, 1, nz], pad=pml)
            call alloc_array(mstability, [1, nx, 1, ny, 1, nz], pad=pml)

            !$omp parallel do private(i, j, k, minv, maxv, dz0) collapse(3)
            do k = max(1, nz1), nz2
                do j = ny1, ny2
                    do i = nx1, nx2

                        select case (aniso_param)

                            case ('iso')
                                call min_max_phase_velocity_3d(c11(i, j, k), c12(i, j, k), c13(i, j, k), c14(i, j, k), c15(i, j, k), c16(i, j, k), &
                                    c22(i, j, k), c23(i, j, k), c24(i, j, k), c25(i, j, k), c26(i, j, k), &
                                    c33(i, j, k), c34(i, j, k), c35(i, j, k), c36(i, j, k), &
                                    c44(i, j, k), c45(i, j, k), c46(i, j, k), &
                                    c55(i, j, k), c56(i, j, k), &
                                    c66(i, j, k), &
                                    rho(i, j, k), 1, 1, maxv, minv)

                            case default
                                ! Get min/max phase velocities for each point; for each point, the velocities are computed at
                                ! 18 polar angles * 18 azimuth angles to ensure accuracy.
                                call min_max_phase_velocity_3d(c11(i, j, k), c12(i, j, k), c13(i, j, k), c14(i, j, k), c15(i, j, k), c16(i, j, k), &
                                    c22(i, j, k), c23(i, j, k), c24(i, j, k), c25(i, j, k), c26(i, j, k), &
                                    c33(i, j, k), c34(i, j, k), c35(i, j, k), c36(i, j, k), &
                                    c44(i, j, k), c45(i, j, k), c46(i, j, k), &
                                    c55(i, j, k), c56(i, j, k), &
                                    c66(i, j, k), &
                                    rho(i, j, k), 18, 18, maxv, minv)
                        end select

                        dz0 = sqrt((1.0/dx + max(0.0, eta_max - eta_zz_i(k))/eta_max*slopex_ixiy(i, j)*1.0/eta_dz_i(k))**2 &
                            + (1.0/dy + (eta_max - eta_zz_i(k))/eta_max*slopey_ixiy(i, j)*1.0/eta_dz_i(k))**2 &
                            + (1.0/eta_dz_i(k))**2)

                        mstability(i, j, k) = 1.0/(temp*maxv*dz0)
                        ! Empirical relation to ensure sufficient accuracy for surface waves
                        if (sum(dz_i(i, j, 1:k)) <= 0.5*rayleigh_wavelength) then
                            mdispersion(i, j, k) = min(0.9*minv/7.0/max(dx, dy), 0.2*0.9*minv/10.0/max(dz_i(i, j, k), dz_h(i, j, k)))
                        else
                            mdispersion(i, j, k) = minv/7.0/max(dx, dy, max(dz_i(i, j, k), dz_h(i, j, k)))
                        end if

                    end do
                end do
            end do
            !$omp end parallel do

            call allreduce_array(mstability)
            call allreduce_array(mdispersion)

            mstability(:, :, -pml + 1:0) = float_huge
            mdispersion(:, :, -pml + 1:0) = float_huge

            dtstable = minval(mstability)
            f0clean = minval(mdispersion)

        else
            ! Regular mesh

            select case (aniso_param)

                case ('iso')
                    minv = minval(vs)
                    maxv = maxval(vp)

                case ('thomsen')
                    minv = minval(vs)
                    maxv = maxval(vp*sqrt(1.0 + 2.0*tieps))

                case ('a-t')
                    minv = minval(vs)
                    maxv = maxval(vp)

                case ('cij')

                    minv = float_huge
                    maxv = -float_huge

                    !$omp parallel do private(i, j, k) collapse(3) reduction(min:minv) reduction(max:maxv)
                    do k = nz1, nz2
                        do j = ny1, ny2
                            do i = nx1, nx2

                                ! Get min/max phase velocities for each point; for each point, the velocities are computed at
                                ! 18 polar angles * 18 azimuth angles to ensure accuracy.
                                call min_max_phase_velocity_3d(c11(i, j, k), c12(i, j, k), c13(i, j, k), c14(i, j, k), c15(i, j, k), c16(i, j, k), &
                                    c22(i, j, k), c23(i, j, k), c24(i, j, k), c25(i, j, k), c26(i, j, k), &
                                    c33(i, j, k), c34(i, j, k), c35(i, j, k), c36(i, j, k), &
                                    c44(i, j, k), c45(i, j, k), c46(i, j, k), &
                                    c55(i, j, k), c56(i, j, k), &
                                    c66(i, j, k), &
                                    rho(i, j, k), 18, 18, tmpmax, tmpmin)

                                minv = min(minv, tmpmin)
                                maxv = max(maxv, tmpmax)

                            end do
                        end do
                    end do
                    !$omp end parallel do

                    minv = group_min(minv)
                    maxv = group_max(maxv)

            end select

            dtstable = 1.0/(temp*maxv*sqrt(1.0/dx**2 + 1.0/dy**2 + 1.0/dz**2))
            f0clean = minv/max(dx, dy, dz)/7.0

        end if

        if (rankid_group == 0) then
            call warn(date_time_compact()//' Stable dt = '//num2str(dtstable, '(es)') &
                //' s, clean f0 = '//num2str(f0clean, '(es)')//' Hz')
        end if

        call check_dt_f0(dt, dtstable, maxval(this%gmtr%srcr(:)%f0), f0clean)

        nt = nint(tmax/dt + 1)

        ! Save generated mesh
        if (yn_free_surface .and. this%save_mesh .and. rankid_group == 0) then

            open(3, file=tidy(this%dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_mesh.txt')
            do i = -pml + 1, nx + pml
                do j = -pml + 1, ny + pml
                    do k = 1, nz + pml
                        write(3, '(5es)') (i - 1)*dx + ox, (j - 1)*dy + oy, zz_i(i, j, k), c11(i, j, k), mstability(i, j, k)
                    end do
                end do
            end do
            close(3)

        end if

        ! Crop models
        call alloc_array(c11, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c11(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c12, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c12(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c13, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c13(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c14, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c14(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c15, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c15(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c16, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c16(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c22, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c22(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c23, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c23(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c24, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c24(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c25, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c25(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c26, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c26(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c33, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c33(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c34, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c34(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c35, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c35(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c36, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c36(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c44, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c44(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c45, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c45(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c46, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c46(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c55, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c55(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c56, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c56(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(c66, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=c66(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(rho, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, source=rho(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))

        ! Compute topography-related constants
        call compute_topography_coefs

        ! Prepare geometry
        sgmtr = this%gmtr
        sgmtr%nx = nx
        sgmtr%ny = ny
        sgmtr%nz = this%nz
        sgmtr%dx = dx
        sgmtr%dy = dy
        sgmtr%dz = this%dz
        sgmtr%ox = ox
        sgmtr%oy = oy
        sgmtr%oz = oz
        sgmtr%xmin = max(ox, this%gmtr%xmin)
        sgmtr%xmax = min((nx - 1)*dx + ox, this%gmtr%xmax)
        sgmtr%ymin = max(oy, this%gmtr%ymin)
        sgmtr%ymax = min((ny - 1)*dy + oy, this%gmtr%ymax)
        sgmtr%zmin = max(oz, this%gmtr%zmin)
        sgmtr%zmax = min((this%nz - 1)*this%dz + oz, this%gmtr%zmax)
        sgmtr%sxmin = max(ox, this%gmtr%sxmin)
        sgmtr%sxmax = min((nx - 1)*dx + ox, this%gmtr%sxmax)
        sgmtr%symin = max(oy, this%gmtr%symin)
        sgmtr%symax = min((ny - 1)*dy + oy, this%gmtr%symax)
        sgmtr%szmin = max(oz, this%gmtr%szmin)
        sgmtr%szmax = min((this%nz - 1)*this%dz + oz, this%gmtr%szmax)
        sgmtr%rxmin = max(ox, this%gmtr%rxmin)
        sgmtr%rxmax = min((nx - 1)*dx + ox, this%gmtr%rxmax)
        sgmtr%rymin = max(oy, this%gmtr%rymin)
        sgmtr%rymax = min((ny - 1)*dy + oy, this%gmtr%rymax)
        sgmtr%rzmin = max(oz, this%gmtr%rzmin)
        sgmtr%rzmax = min((this%nz - 1)*this%dz + oz, this%gmtr%rzmax)

        if (yn_free_surface) then

            topox = meshgrid([nx + 2*pml + 2, ny + 2*pml + 2], [dx, dy], [-(pml + 1)*dx + ox, -(pml + 1)*dy + oy], dim=1)
            topoy = meshgrid([nx + 2*pml + 2, ny + 2*pml + 2], [dx, dy], [-(pml + 1)*dx + ox, -(pml + 1)*dy + oy], dim=2)

            ! If measure from the free surface then + topo (here, topo has already been sign-corrected)
            ! e.g., sz = 400, topo = 100, then real depth = 500
            stp = ginterp(topox, topoy, flatten(zz_i(:, :, 1)), sgmtr%srcr(:)%x, sgmtr%srcr(:)%y)
            rtp = ginterp(topox, topoy, flatten(zz_i(:, :, 1)), sgmtr%recr(:)%x, sgmtr%recr(:)%y)

            if (this%measure_source_depth_from_surface) then
                sgmtr%srcr(:)%z = sgmtr%srcr(:)%z + stp
            end if
            if (this%measure_receiver_depth_from_surface) then
                sgmtr%recr(:)%z = sgmtr%recr(:)%z + rtp
            end if

            ! Map real location to virtual mesh
            sgmtr%srcr(:)%z = (sgmtr%srcr(:)%z - stp)/(-stp + depth_max)*eta_max
            sgmtr%recr(:)%z = (sgmtr%recr(:)%z - rtp)/(-rtp + depth_max)*eta_max

            ! Avoid placing source on the free surface as otherwise the resulting amplitude
            ! can be problematic
            where (sgmtr%srcr(:)%z < eta_dz_i(1))
                sgmtr%srcr(:)%z = eta_dz_i(1)
                sgmtr%srcr(:)%amp = sgmtr%srcr(:)%amp*1.2
            end where

            call alloc_array(sgmtr%z_i, [1, nx, 1, ny, 1, nz], pad=pml + 1)
            call alloc_array(sgmtr%z_h, [1, nx, 1, ny, 1, nz], pad=pml + 1)
            call alloc_array(sgmtr%dz_i, [1, nx, 1, ny, 1, nz], pad=pml + 1)
            call alloc_array(sgmtr%dz_h, [1, nx, 1, ny, 1, nz], pad=pml + 1)
            !$omp parallel do private(i, j) collapse(2)
            do j = -pml, ny + pml + 1
                do i = -pml, nx + pml + 1
                    sgmtr%z_i(i, j, :) = eta_zz_i
                    sgmtr%z_h(i, j, :) = eta_zz_h
                    sgmtr%dz_i(i, j, :) = eta_dz_i
                    sgmtr%dz_h(i, j, :) = eta_dz_h
                end do
            end do
            !$omp end parallel do

        end if

        call sgmtr%prepare_geometry

        sgmtr%dt = dt
        sgmtr%nt = nt
        do i = 1, sgmtr%ns
            if (sgmtr%srcr(i)%mechanism == 'explosion' .or. sgmtr%srcr(i)%mechanism == 'mt') then
                sgmtr%srcr(i)%time_integration = -1
            end if
        end do
        call sgmtr%prepare_stf

        data_dt = this%data_dt
        data_tmax = this%data_tmax
        data_nt = nint(data_tmax/data_dt + 1)

        cc_step_interval = this%cc_step_interval
        verbose = this%verbose
        yn_reconstruct = this%reconstruct

        dir_synthetic = this%dir_synthetic
        dir_adjoint = this%dir_adjoint
        dir_working = this%dir_working

        if (rankid_group == 0) then
            call make_directory(dir_synthetic)
        end if

        snaps = this%snaps
        np = size(snaps)
        if (np > 0) then
            dir_snapshot = this%dir_snapshot
            if (rankid_group == 0) then
                call make_directory(dir_snapshot)
            end if
        end if

        yn_compx = this%compx
        yn_compy = this%compy
        yn_compz = this%compz

        ! If mt inversion is required
        yn_update_medium = this%yn_update_medium
        yn_update_source = this%yn_update_source

        nc_mt = this%nc_mt
        if (norm2(this%mt) > 0) then

            dstf_dt = zeros(nt, sgmtr%ns)

            do i = 1, sgmtr%ns
                sgmtr%srcr(i)%mechanism = 'mt'
                sgmtr%srcr(i)%moment_tensor(1, 1) = this%mt(1)
                sgmtr%srcr(i)%moment_tensor(2, 2) = this%mt(2)
                sgmtr%srcr(i)%moment_tensor(3, 3) = this%mt(3)
                sgmtr%srcr(i)%moment_tensor(1, 2) = this%mt(4)
                sgmtr%srcr(i)%moment_tensor(1, 3) = this%mt(5)
                sgmtr%srcr(i)%moment_tensor(2, 3) = this%mt(6)

                nbeg = nint(sgmtr%srcr(i)%t0/dt) + 1
                nend = nbeg + sgmtr%srcr(i)%nt - 1
                dstf_dt(nbeg:nend, i) = sgmtr%srcr(i)%stf
                dstf_dt(:, i) = deriv(dstf_dt(:, i))

            end do

        end if

        ! Split traces into group processes
        call alloc_array(trace_range, [0, nrank_group - 1, 1, 2])
        call cut(1, sgmtr%nr, nrank_group, trace_range)

        call mpibarrier_group

    end subroutine prepare_modeling

    !
    !> Allocate memory
    !
    subroutine alloc_forward_wavefield

        ! Stress -- set a
        call alloc_array(stressxx_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyy_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszz_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxy_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxz_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyz_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxvx_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvy_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvz_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvx_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvy_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvz_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvx_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvy_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvz_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Stress -- set b
        call alloc_array(stressxx_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyy_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszz_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxy_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxz_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyz_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxvx_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvy_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvz_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvx_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvy_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvz_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvx_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvy_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvz_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Stress -- set c
        call alloc_array(stressxx_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyy_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszz_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxy_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxz_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyz_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxvx_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvy_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvz_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvx_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvy_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvz_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvx_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvy_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvz_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Stress -- set d
        call alloc_array(stressxx_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyy_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszz_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxy_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxz_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyz_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxvx_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvy_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvz_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvx_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvy_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvz_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvx_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvy_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvz_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Velocity -- set a
        call alloc_array(vx_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vy_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vz_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxxx_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxy_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxz_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyxy_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyy_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyz_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxx_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxy_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxz_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyy_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyz_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzzz_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Velocity -- set b
        call alloc_array(vx_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vy_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vz_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxxx_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxy_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxz_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyxy_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyy_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyz_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxx_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxy_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxz_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyy_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyz_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzzz_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Velocity -- set c
        call alloc_array(vx_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vy_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vz_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxxx_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxy_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxz_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyxy_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyy_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyz_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxx_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxy_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxz_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyy_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyz_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzzz_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Velocity -- set d
        call alloc_array(vx_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vy_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vz_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxxx_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxy_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxz_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyxy_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyy_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyz_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxx_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxy_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxz_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyy_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyz_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzzz_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])

    end subroutine alloc_forward_wavefield

    subroutine alloc_adjoint_wavefield

        ! Stress -- set a
        call alloc_array(stressxxr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyyr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszzr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxyr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxzr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyzr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxvxr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvyr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvzr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvxr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvyr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvzr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvxr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvyr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvzr_ixiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Stress -- set b
        call alloc_array(stressxxr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyyr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszzr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxyr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxzr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyzr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxvxr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvyr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvzr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvxr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvyr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvzr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvxr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvyr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvzr_hxhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Stress -- set c
        call alloc_array(stressxxr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyyr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszzr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxyr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxzr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyzr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxvxr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvyr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvzr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvxr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvyr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvzr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvxr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvyr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvzr_hxiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Stress -- set d
        call alloc_array(stressxxr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyyr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stresszzr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxyr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressxzr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(stressyzr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxvxr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvyr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxvzr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvxr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvyr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyvzr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvxr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvyr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzvzr_ixhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Velocity -- set a
        call alloc_array(vxr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vyr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vzr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxxxr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxyr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxzr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyxyr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyyr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyzr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxxr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxyr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxzr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyyr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyzr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzzzr_hxiyiz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Velocity -- set b
        call alloc_array(vxr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vyr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vzr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxxxr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxyr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxzr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyxyr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyyr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyzr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxxr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxyr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxzr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyyr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyzr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzzzr_ixhyiz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Velocity -- set c
        call alloc_array(vxr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vyr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vzr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxxxr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxyr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxzr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyxyr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyyr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyzr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxxr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxyr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxzr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyyr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyzr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzzzr_ixiyhz, [nx1, nx2, ny1, ny2, nz1, nz2])

        ! Velocity -- set d
        call alloc_array(vxr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vyr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vzr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

        call alloc_array(memory_pdxxxr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxyr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdxxzr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyxyr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyyr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdyyzr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxxr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxyr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzxzr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyyr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzyzr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(memory_pdzzzr_hxhyhz, [nx1, nx2, ny1, ny2, nz1, nz2])

    end subroutine

    !
    !> Check if a point is in a domain decomposition block
    !
    pure function is_in_block(i, j, k) result(y)

        integer, intent(in) :: i, j, k
        logical :: y

        y = (i >= nx1_interior .and. i <= nx2_interior &
            .and. j >= ny1_interior .and. j <= ny2_interior &
            .and. k >= nz1_interior .and. k <= nz2_interior)

    end function

    !
    !> Precompute free-surface related coefficients
    !
    subroutine compute_topography_coefs

        real :: eta, tpcz
        real :: hx, hy
        integer :: i, j, k
        real :: c11_ixiyiz, c12_ixiyiz, c13_ixiyiz, c14_ixiyiz, c15_ixiyiz, c16_ixiyiz
        real :: c22_ixiyiz, c23_ixiyiz, c24_ixiyiz, c25_ixiyiz, c26_ixiyiz
        real :: c33_ixiyiz, c34_ixiyiz, c35_ixiyiz, c36_ixiyiz
        real :: c44_ixiyiz, c45_ixiyiz, c46_ixiyiz
        real :: c55_ixiyiz, c56_ixiyiz
        real :: c66_ixiyiz
        real :: c11_hxhyiz, c12_hxhyiz, c13_hxhyiz, c14_hxhyiz, c15_hxhyiz, c16_hxhyiz
        real :: c22_hxhyiz, c23_hxhyiz, c24_hxhyiz, c25_hxhyiz, c26_hxhyiz
        real :: c33_hxhyiz, c34_hxhyiz, c35_hxhyiz, c36_hxhyiz
        real :: c44_hxhyiz, c45_hxhyiz, c46_hxhyiz
        real :: c55_hxhyiz, c56_hxhyiz
        real :: c66_hxhyiz

        if (.not. (nz1 <= 1 .and. 1 <= nz2) .or. .not. yn_free_surface) then
            return
        end if

        alloc2(coefiix, nx1, nx2, ny1, ny2)
        alloc2(coefiiy, nx1, nx2, ny1, ny2)
        alloc2(coefiiz, nx1, nx2, ny1, ny2)
        alloc2(coefiizr, nx1, nx2, ny1, ny2)
        alloc2(coefhhx, nx1, nx2, ny1, ny2)
        alloc2(coefhhy, nx1, nx2, ny1, ny2)
        alloc2(coefhhz, nx1, nx2, ny1, ny2)
        alloc2(coefhhzr, nx1, nx2, ny1, ny2)

        !$omp parallel do private(i, j) collapse(2)
        do j = ny1, ny2
            do i = nx1, nx2

                coefiix(i, j)%array = zeros(3, 3)
                coefiiy(i, j)%array = zeros(3, 3)
                coefiiz(i, j)%array = zeros(3, 3)
                coefiizr(i, j)%array = zeros(3, 3)
                coefhhx(i, j)%array = zeros(3, 3)
                coefhhy(i, j)%array = zeros(3, 3)
                coefhhz(i, j)%array = zeros(3, 3)
                coefhhzr(i, j)%array = zeros(3, 3)

            end do
        end do
        !$omp end parallel do

        k = 1
        !$omp parallel do private(i, j, eta, hx, hy, tpcz, &
            !$omp c11_ixiyiz, c12_ixiyiz, c13_ixiyiz, c14_ixiyiz, c15_ixiyiz, c16_ixiyiz, &
            !$omp c22_ixiyiz, c23_ixiyiz, c24_ixiyiz, c25_ixiyiz, c26_ixiyiz, &
            !$omp c33_ixiyiz, c34_ixiyiz, c35_ixiyiz, c36_ixiyiz, &
            !$omp c44_ixiyiz, c45_ixiyiz, c46_ixiyiz, &
            !$omp c55_ixiyiz, c56_ixiyiz, &
            !$omp c66_ixiyiz) collapse(2)
        do j = ny1, ny2
            do i = nx1, nx2

                tpcz = eta_max/(topo_ixiy(i, j) + depth_max)

                c11_ixiyiz = c11(i, j, k)
                c12_ixiyiz = c12(i, j, k)
                c13_ixiyiz = c13(i, j, k)
                c14_ixiyiz = c14(i, j, k)
                c15_ixiyiz = c15(i, j, k)
                c16_ixiyiz = c16(i, j, k)
                c22_ixiyiz = c22(i, j, k)
                c23_ixiyiz = c23(i, j, k)
                c24_ixiyiz = c24(i, j, k)
                c25_ixiyiz = c25(i, j, k)
                c26_ixiyiz = c26(i, j, k)
                c33_ixiyiz = c33(i, j, k)
                c34_ixiyiz = c34(i, j, k)
                c35_ixiyiz = c35(i, j, k)
                c36_ixiyiz = c36(i, j, k)
                c44_ixiyiz = c44(i, j, k)
                c45_ixiyiz = c45(i, j, k)
                c46_ixiyiz = c46(i, j, k)
                c55_ixiyiz = c55(i, j, k)
                c56_ixiyiz = c56(i, j, k)
                c66_ixiyiz = c66(i, j, k)

                hx = slopex_ixiy(i, j)
                hy = slopey_ixiy(i, j)

                coefiix(i, j)%array(1, 1) = -c15_ixiyiz - c11_ixiyiz*hx - c16_ixiyiz*hy
                coefiix(i, j)%array(1, 2) = -c56_ixiyiz - c16_ixiyiz*hx - c66_ixiyiz*hy
                coefiix(i, j)%array(1, 3) = -c55_ixiyiz - c15_ixiyiz*hx - c56_ixiyiz*hy
                coefiix(i, j)%array(2, 1) = -c14_ixiyiz - c16_ixiyiz*hx - c12_ixiyiz*hy
                coefiix(i, j)%array(2, 2) = -c46_ixiyiz - c66_ixiyiz*hx - c26_ixiyiz*hy
                coefiix(i, j)%array(2, 3) = -c45_ixiyiz - c56_ixiyiz*hx - c25_ixiyiz*hy
                coefiix(i, j)%array(3, 1) = -c13_ixiyiz - c15_ixiyiz*hx - c14_ixiyiz*hy
                coefiix(i, j)%array(3, 2) = -c36_ixiyiz - c56_ixiyiz*hx - c46_ixiyiz*hy
                coefiix(i, j)%array(3, 3) = -c35_ixiyiz - c55_ixiyiz*hx - c45_ixiyiz*hy

                coefiiy(i, j)%array(1, 1) = -c56_ixiyiz - c16_ixiyiz*hx - c66_ixiyiz*hy
                coefiiy(i, j)%array(1, 2) = -c25_ixiyiz - c12_ixiyiz*hx - c26_ixiyiz*hy
                coefiiy(i, j)%array(1, 3) = -c45_ixiyiz - c14_ixiyiz*hx - c46_ixiyiz*hy
                coefiiy(i, j)%array(2, 1) = -c46_ixiyiz - c66_ixiyiz*hx - c26_ixiyiz*hy
                coefiiy(i, j)%array(2, 2) = -c24_ixiyiz - c26_ixiyiz*hx - c22_ixiyiz*hy
                coefiiy(i, j)%array(2, 3) = -c44_ixiyiz - c46_ixiyiz*hx - c24_ixiyiz*hy
                coefiiy(i, j)%array(3, 1) = -c36_ixiyiz - c56_ixiyiz*hx - c46_ixiyiz*hy
                coefiiy(i, j)%array(3, 2) = -c23_ixiyiz - c25_ixiyiz*hx - c24_ixiyiz*hy
                coefiiy(i, j)%array(3, 3) = -c34_ixiyiz - c45_ixiyiz*hx - c44_ixiyiz*hy

                coefiiz(i, j)%array(1, 1) = c55_ixiyiz*tpcz + 2*c15_ixiyiz*hx*tpcz &
                    + c11_ixiyiz*hx**2*tpcz + 2*c56_ixiyiz*hy*tpcz + 2*c16_ixiyiz*hx*hy*tpcz + c66_ixiyiz*hy**2*tpcz
                coefiiz(i, j)%array(1, 2) = c45_ixiyiz*tpcz + c14_ixiyiz*hx*tpcz &
                    + c56_ixiyiz*hx*tpcz + c16_ixiyiz*hx**2*tpcz + c25_ixiyiz*hy*tpcz &
                    + c46_ixiyiz*hy*tpcz + c12_ixiyiz*hx*hy*tpcz + c66_ixiyiz*hx*hy*tpcz + c26_ixiyiz*hy**2*tpcz
                coefiiz(i, j)%array(1, 3) = c35_ixiyiz*tpcz + c13_ixiyiz*hx*tpcz &
                    + c55_ixiyiz*hx*tpcz + c15_ixiyiz*hx**2*tpcz + c36_ixiyiz*hy*tpcz &
                    + c45_ixiyiz*hy*tpcz + c14_ixiyiz*hx*hy*tpcz + c56_ixiyiz*hx*hy*tpcz + c46_ixiyiz*hy**2*tpcz
                coefiiz(i, j)%array(2, 1) = c45_ixiyiz*tpcz + c14_ixiyiz*hx*tpcz &
                    + c56_ixiyiz*hx*tpcz + c16_ixiyiz*hx**2*tpcz + c25_ixiyiz*hy*tpcz &
                    + c46_ixiyiz*hy*tpcz + c12_ixiyiz*hx*hy*tpcz + c66_ixiyiz*hx*hy*tpcz + c26_ixiyiz*hy**2*tpcz
                coefiiz(i, j)%array(2, 2) = c44_ixiyiz*tpcz + 2*c46_ixiyiz*hx*tpcz &
                    + c66_ixiyiz*hx**2*tpcz + 2*c24_ixiyiz*hy*tpcz + 2*c26_ixiyiz*hx*hy*tpcz + c22_ixiyiz*hy**2*tpcz
                coefiiz(i, j)%array(2, 3) = c34_ixiyiz*tpcz + c36_ixiyiz*hx*tpcz &
                    + c45_ixiyiz*hx*tpcz + c56_ixiyiz*hx**2*tpcz + c23_ixiyiz*hy*tpcz &
                    + c44_ixiyiz*hy*tpcz + c25_ixiyiz*hx*hy*tpcz + c46_ixiyiz*hx*hy*tpcz + c24_ixiyiz*hy**2*tpcz
                coefiiz(i, j)%array(3, 1) = c35_ixiyiz*tpcz + c13_ixiyiz*hx*tpcz &
                    + c55_ixiyiz*hx*tpcz + c15_ixiyiz*hx**2*tpcz + c36_ixiyiz*hy*tpcz &
                    + c45_ixiyiz*hy*tpcz + c14_ixiyiz*hx*hy*tpcz + c56_ixiyiz*hx*hy*tpcz + c46_ixiyiz*hy**2*tpcz
                coefiiz(i, j)%array(3, 2) = c34_ixiyiz*tpcz + c36_ixiyiz*hx*tpcz &
                    + c45_ixiyiz*hx*tpcz + c56_ixiyiz*hx**2*tpcz + c23_ixiyiz*hy*tpcz &
                    + c44_ixiyiz*hy*tpcz + c25_ixiyiz*hx*hy*tpcz + c46_ixiyiz*hx*hy*tpcz + c24_ixiyiz*hy**2*tpcz
                coefiiz(i, j)%array(3, 3) = c33_ixiyiz*tpcz + 2*c35_ixiyiz*hx*tpcz &
                    + c55_ixiyiz*hx**2*tpcz + 2*c34_ixiyiz*hy*tpcz + 2*c45_ixiyiz*hx*hy*tpcz + c44_ixiyiz*hy**2*tpcz

                coefiix(i, j)%array = solve(coefiiz(i, j)%array, coefiix(i, j)%array)
                coefiiy(i, j)%array = solve(coefiiz(i, j)%array, coefiiy(i, j)%array)

            end do
        end do
        !$omp end parallel do

        k = 1
        !$omp parallel do private(i, j, eta, hx, hy, tpcz, &
            !$omp c11_hxhyiz, c12_hxhyiz, c13_hxhyiz, c14_hxhyiz, c15_hxhyiz, c16_hxhyiz, &
            !$omp c22_hxhyiz, c23_hxhyiz, c24_hxhyiz, c25_hxhyiz, c26_hxhyiz, &
            !$omp c33_hxhyiz, c34_hxhyiz, c35_hxhyiz, c36_hxhyiz, &
            !$omp c44_hxhyiz, c45_hxhyiz, c46_hxhyiz, &
            !$omp c55_hxhyiz, c56_hxhyiz, &
            !$omp c66_hxhyiz) collapse(2)
        do j = ny1 - 1, ny2 - 1
            do i = nx1 - 1, nx2 - 1

                tpcz = eta_max/(topo_hxhy(i + 1, j + 1) + depth_max)

                c11_hxhyiz = 0.25*sum(c11(i:i + 1, j:j + 1, k))
                c12_hxhyiz = 0.25*sum(c12(i:i + 1, j:j + 1, k))
                c13_hxhyiz = 0.25*sum(c13(i:i + 1, j:j + 1, k))
                c14_hxhyiz = 0.25*sum(c14(i:i + 1, j:j + 1, k))
                c15_hxhyiz = 0.25*sum(c15(i:i + 1, j:j + 1, k))
                c16_hxhyiz = 0.25*sum(c16(i:i + 1, j:j + 1, k))
                c22_hxhyiz = 0.25*sum(c22(i:i + 1, j:j + 1, k))
                c23_hxhyiz = 0.25*sum(c23(i:i + 1, j:j + 1, k))
                c24_hxhyiz = 0.25*sum(c24(i:i + 1, j:j + 1, k))
                c25_hxhyiz = 0.25*sum(c25(i:i + 1, j:j + 1, k))
                c26_hxhyiz = 0.25*sum(c26(i:i + 1, j:j + 1, k))
                c33_hxhyiz = 0.25*sum(c33(i:i + 1, j:j + 1, k))
                c34_hxhyiz = 0.25*sum(c34(i:i + 1, j:j + 1, k))
                c35_hxhyiz = 0.25*sum(c35(i:i + 1, j:j + 1, k))
                c36_hxhyiz = 0.25*sum(c36(i:i + 1, j:j + 1, k))
                c44_hxhyiz = 0.25*sum(c44(i:i + 1, j:j + 1, k))
                c45_hxhyiz = 0.25*sum(c45(i:i + 1, j:j + 1, k))
                c46_hxhyiz = 0.25*sum(c46(i:i + 1, j:j + 1, k))
                c55_hxhyiz = 0.25*sum(c55(i:i + 1, j:j + 1, k))
                c56_hxhyiz = 0.25*sum(c56(i:i + 1, j:j + 1, k))
                c66_hxhyiz = 0.25*sum(c66(i:i + 1, j:j + 1, k))

                hx = slopex_hxhy(i + 1, j + 1)
                hy = slopey_hxhy(i + 1, j + 1)

                coefhhx(i + 1, j + 1)%array(1, 1) = -c15_hxhyiz - c11_hxhyiz*hx - c16_hxhyiz*hy
                coefhhx(i + 1, j + 1)%array(1, 2) = -c56_hxhyiz - c16_hxhyiz*hx - c66_hxhyiz*hy
                coefhhx(i + 1, j + 1)%array(1, 3) = -c55_hxhyiz - c15_hxhyiz*hx - c56_hxhyiz*hy
                coefhhx(i + 1, j + 1)%array(2, 1) = -c14_hxhyiz - c16_hxhyiz*hx - c12_hxhyiz*hy
                coefhhx(i + 1, j + 1)%array(2, 2) = -c46_hxhyiz - c66_hxhyiz*hx - c26_hxhyiz*hy
                coefhhx(i + 1, j + 1)%array(2, 3) = -c45_hxhyiz - c56_hxhyiz*hx - c25_hxhyiz*hy
                coefhhx(i + 1, j + 1)%array(3, 1) = -c13_hxhyiz - c15_hxhyiz*hx - c14_hxhyiz*hy
                coefhhx(i + 1, j + 1)%array(3, 2) = -c36_hxhyiz - c56_hxhyiz*hx - c46_hxhyiz*hy
                coefhhx(i + 1, j + 1)%array(3, 3) = -c35_hxhyiz - c55_hxhyiz*hx - c45_hxhyiz*hy

                coefhhy(i + 1, j + 1)%array(1, 1) = -c56_hxhyiz - c16_hxhyiz*hx - c66_hxhyiz*hy
                coefhhy(i + 1, j + 1)%array(1, 2) = -c25_hxhyiz - c12_hxhyiz*hx - c26_hxhyiz*hy
                coefhhy(i + 1, j + 1)%array(1, 3) = -c45_hxhyiz - c14_hxhyiz*hx - c46_hxhyiz*hy
                coefhhy(i + 1, j + 1)%array(2, 1) = -c46_hxhyiz - c66_hxhyiz*hx - c26_hxhyiz*hy
                coefhhy(i + 1, j + 1)%array(2, 2) = -c24_hxhyiz - c26_hxhyiz*hx - c22_hxhyiz*hy
                coefhhy(i + 1, j + 1)%array(2, 3) = -c44_hxhyiz - c46_hxhyiz*hx - c24_hxhyiz*hy
                coefhhy(i + 1, j + 1)%array(3, 1) = -c36_hxhyiz - c56_hxhyiz*hx - c46_hxhyiz*hy
                coefhhy(i + 1, j + 1)%array(3, 2) = -c23_hxhyiz - c25_hxhyiz*hx - c24_hxhyiz*hy
                coefhhy(i + 1, j + 1)%array(3, 3) = -c34_hxhyiz - c45_hxhyiz*hx - c44_hxhyiz*hy

                coefhhz(i + 1, j + 1)%array(1, 1) = c55_hxhyiz*tpcz + 2*c15_hxhyiz*hx*tpcz &
                    + c11_hxhyiz*hx**2*tpcz + 2*c56_hxhyiz*hy*tpcz + 2*c16_hxhyiz*hx*hy*tpcz + c66_hxhyiz*hy**2*tpcz
                coefhhz(i + 1, j + 1)%array(1, 2) = c45_hxhyiz*tpcz + c14_hxhyiz*hx*tpcz &
                    + c56_hxhyiz*hx*tpcz + c16_hxhyiz*hx**2*tpcz + c25_hxhyiz*hy*tpcz &
                    + c46_hxhyiz*hy*tpcz + c12_hxhyiz*hx*hy*tpcz + c66_hxhyiz*hx*hy*tpcz + c26_hxhyiz*hy**2*tpcz
                coefhhz(i + 1, j + 1)%array(1, 3) = c35_hxhyiz*tpcz + c13_hxhyiz*hx*tpcz &
                    + c55_hxhyiz*hx*tpcz + c15_hxhyiz*hx**2*tpcz + c36_hxhyiz*hy*tpcz &
                    + c45_hxhyiz*hy*tpcz + c14_hxhyiz*hx*hy*tpcz + c56_hxhyiz*hx*hy*tpcz + c46_hxhyiz*hy**2*tpcz
                coefhhz(i + 1, j + 1)%array(2, 1) = c45_hxhyiz*tpcz + c14_hxhyiz*hx*tpcz &
                    + c56_hxhyiz*hx*tpcz + c16_hxhyiz*hx**2*tpcz + c25_hxhyiz*hy*tpcz &
                    + c46_hxhyiz*hy*tpcz + c12_hxhyiz*hx*hy*tpcz + c66_hxhyiz*hx*hy*tpcz + c26_hxhyiz*hy**2*tpcz
                coefhhz(i + 1, j + 1)%array(2, 2) = c44_hxhyiz*tpcz + 2*c46_hxhyiz*hx*tpcz &
                    + c66_hxhyiz*hx**2*tpcz + 2*c24_hxhyiz*hy*tpcz + 2*c26_hxhyiz*hx*hy*tpcz + c22_hxhyiz*hy**2*tpcz
                coefhhz(i + 1, j + 1)%array(2, 3) = c34_hxhyiz*tpcz + c36_hxhyiz*hx*tpcz &
                    + c45_hxhyiz*hx*tpcz + c56_hxhyiz*hx**2*tpcz + c23_hxhyiz*hy*tpcz &
                    + c44_hxhyiz*hy*tpcz + c25_hxhyiz*hx*hy*tpcz + c46_hxhyiz*hx*hy*tpcz + c24_hxhyiz*hy**2*tpcz
                coefhhz(i + 1, j + 1)%array(3, 1) = c35_hxhyiz*tpcz + c13_hxhyiz*hx*tpcz &
                    + c55_hxhyiz*hx*tpcz + c15_hxhyiz*hx**2*tpcz + c36_hxhyiz*hy*tpcz &
                    + c45_hxhyiz*hy*tpcz + c14_hxhyiz*hx*hy*tpcz + c56_hxhyiz*hx*hy*tpcz + c46_hxhyiz*hy**2*tpcz
                coefhhz(i + 1, j + 1)%array(3, 2) = c34_hxhyiz*tpcz + c36_hxhyiz*hx*tpcz &
                    + c45_hxhyiz*hx*tpcz + c56_hxhyiz*hx**2*tpcz + c23_hxhyiz*hy*tpcz &
                    + c44_hxhyiz*hy*tpcz + c25_hxhyiz*hx*hy*tpcz + c46_hxhyiz*hx*hy*tpcz + c24_hxhyiz*hy**2*tpcz
                coefhhz(i + 1, j + 1)%array(3, 3) = c33_hxhyiz*tpcz + 2*c35_hxhyiz*hx*tpcz &
                    + c55_hxhyiz*hx**2*tpcz + 2*c34_hxhyiz*hy*tpcz + 2*c45_hxhyiz*hx*hy*tpcz + c44_hxhyiz*hy**2*tpcz

                coefhhx(i + 1, j + 1)%array = solve(coefhhz(i + 1, j + 1)%array, coefhhx(i + 1, j + 1)%array)
                coefhhy(i + 1, j + 1)%array = solve(coefhhz(i + 1, j + 1)%array, coefhhy(i + 1, j + 1)%array)

            end do
        end do
        !$omp end parallel do

    end subroutine compute_topography_coefs

end module
