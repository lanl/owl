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

module acoustic_iso_2d_vars

    use libflit
    use mod_fdcoef
    use mod_su
    use mod_source_receiver, only: source_receiver_geometry
    use mod_utility, only: check_dt_f0
    use mod_source_receiver, only: nkw

    use acoustic_iso_2d

    implicit none

    real, allocatable, dimension(:, :) :: p, vx, vz
    real, allocatable, dimension(:, :) :: memory_pdxp, memory_pdzp
    real, allocatable, dimension(:, :) :: memory_pdxvx, memory_pdzvz

    real, allocatable, dimension(:, :) :: pr, vxr, vzr
    real, allocatable, dimension(:, :) :: memory_pdxpr, memory_pdzpr
    real, allocatable, dimension(:, :) :: memory_pdxvxr, memory_pdzvzr

    real, allocatable, dimension(:, :) :: src_p, rec_p
    real, allocatable, dimension(:, :) :: src_vx, src_vz, rec_vx, rec_vz

    real, allocatable, dimension(:, :) :: energy_src_v, energy_rec_v
    real, allocatable, dimension(:, :) :: energy_src_a, energy_rec_a

    real, allocatable, dimension(:, :) :: prev_p
    real, allocatable, dimension(:, :) :: vp, rho, bk
    real, allocatable, dimension(:, :) :: prev_vx, prev_vz

    real, allocatable, dimension(:) :: p_lrsh, p_lrrh, p_lrarh
    real, allocatable, dimension(:) :: p_udsh, p_udrh, p_udarh

    real, allocatable, dimension(:, :) :: grad_vp, grad_rho

    real, allocatable, dimension(:) :: snaps

    integer :: np
    integer :: nomega, iw
    complex :: cmpcoef

    real :: idx, idz

    integer :: nx, nz
    real :: dx, dz
    real :: ox, oz

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

contains

    subroutine prepare_modeling(this)

        type(wave_solver_acoustic_iso_2d), intent(in) :: this

        real :: minv, maxv, dtstable, f0clean, temp

        nx = this%nx
        nz = this%nz
        dx = this%dx
        dz = this%dz
        ox = this%ox
        oz = this%oz

        idx = 1.0d0/dx
        idz = 1.0d0/dz

        pml = this%pml

        ! Medium parameter models
        vp = transpose(this%vp)
        rho = transpose(this%rho)
        bk = rho*vp**2

        call pad_array(vp, [pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(rho, [pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(bk, [pml + 1, pml + 1, pml + 1, pml + 1])

        ! Check stability and dispersion
        dt = this%dt
        tmax = this%tmax

        minv = minval(vp)
        maxv = maxval(vp)
        temp = sum(abs(fdcoefs))

        dtstable = 1.0/(temp*maxv*sqrt(1.0/dx**2 + 1.0/dz**2))
        f0clean = minv/max(dx, dz)/7.0

        if (rankid_group == 0) then
            call warn(date_time_compact()//' Stable dt = '//num2str(dtstable, '(es)') &
                //' s, clean f0 = '//num2str(f0clean, '(es)')//' Hz')
        end if

        call check_dt_f0(dt, dtstable, maxval(this%gmtr%srcr(:)%f0), f0clean)
        nt = nint(tmax/dt + 1)

        ! Prepare geometry
        sgmtr = this%gmtr
        sgmtr%nx = nx
        sgmtr%nz = nz
        sgmtr%dx = dx
        sgmtr%dz = dz
        sgmtr%ox = ox
        sgmtr%oz = oz
        sgmtr%xmin = max(ox, this%gmtr%xmin)
        sgmtr%xmax = min((nx - 1)*dx + ox, this%gmtr%xmax)
        sgmtr%zmin = max(oz, this%gmtr%zmin)
        sgmtr%zmax = min((nz - 1)*dz + oz, this%gmtr%zmax)
        sgmtr%sxmin = max(ox, this%gmtr%sxmin)
        sgmtr%sxmax = min((nx - 1)*dx + ox, this%gmtr%sxmax)
        sgmtr%szmin = max(oz, this%gmtr%szmin)
        sgmtr%szmax = min((nz - 1)*dz + oz, this%gmtr%szmax)
        sgmtr%rxmin = max(ox, this%gmtr%rxmin)
        sgmtr%rxmax = min((nx - 1)*dx + ox, this%gmtr%rxmax)
        sgmtr%rzmin = max(oz, this%gmtr%rzmin)
        sgmtr%rzmax = min((nz - 1)*dz + oz, this%gmtr%rzmax)
        call sgmtr%prepare_geometry

        sgmtr%dt = dt
        sgmtr%nt = nt
        call sgmtr%prepare_stf

        data_dt = this%data_dt
        data_tmax = this%data_tmax
        data_nt = nint(data_tmax/data_dt + 1)

        cc_step_interval = this%cc_step_interval
        verbose = this%verbose
        yn_free_surface = this%free_surface
        yn_reconstruct = this%reconstruct

        dir_synthetic = this%dir_synthetic
        dir_adjoint = this%dir_adjoint
        dir_working = this%dir_working

        call make_directory(dir_synthetic)

        snaps = this%snaps
        np = size(snaps)
        if (np > 0) then
            dir_snapshot = this%dir_snapshot
            call make_directory(dir_snapshot)
        end if

    end subroutine

    subroutine alloc_forward_wavefield

        call alloc_array(vx, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(vz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(p, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxvx, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzvz, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxp, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzp, [1, nx, 1, nz], pad=pml + fdhalf)

    end subroutine

    subroutine alloc_adjoint_wavefield

        call alloc_array(vxr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(vzr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(pr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxvxr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzvzr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdxpr, [1, nx, 1, nz], pad=pml + fdhalf)
        call alloc_array(memory_pdzpr, [1, nx, 1, nz], pad=pml + fdhalf)

    end subroutine

end module
