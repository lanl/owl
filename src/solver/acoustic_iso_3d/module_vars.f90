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

module acoustic_iso_3d_vars

    use libflit
    use mod_fdcoef
    use mod_su
    use mod_source_receiver, only: source_receiver_geometry
    use mod_utility, only: check_dt_f0
    use mod_source_receiver, only: nkw

    use acoustic_iso_3d

    implicit none

    real, allocatable, dimension(:, :, :) :: p, vx, vy, vz
    real, allocatable, dimension(:, :, :) :: memory_pdxp_xmin, memory_pdxp_xmax
    real, allocatable, dimension(:, :, :) :: memory_pdyp_ymin, memory_pdyp_ymax
    real, allocatable, dimension(:, :, :) :: memory_pdzp_zmin, memory_pdzp_zmax
    real, allocatable, dimension(:, :, :) :: memory_pdxvx_xmin, memory_pdxvx_xmax
    real, allocatable, dimension(:, :, :) :: memory_pdyvy_ymin, memory_pdyvy_ymax
    real, allocatable, dimension(:, :, :) :: memory_pdzvz_zmin, memory_pdzvz_zmax

    real, allocatable, dimension(:, :, :) :: pr, vxr, vyr, vzr
    real, allocatable, dimension(:, :, :) :: memory_pdxpr_xmin, memory_pdxpr_xmax
    real, allocatable, dimension(:, :, :) :: memory_pdypr_ymin, memory_pdypr_ymax
    real, allocatable, dimension(:, :, :) :: memory_pdzpr_zmin, memory_pdzpr_zmax
    real, allocatable, dimension(:, :, :) :: memory_pdxvxr_xmin, memory_pdxvxr_xmax
    real, allocatable, dimension(:, :, :) :: memory_pdyvyr_ymin, memory_pdyvyr_ymax
    real, allocatable, dimension(:, :, :) :: memory_pdzvzr_zmin, memory_pdzvzr_zmax

    real, allocatable, dimension(:, :, :) :: energy_src, energy_rec

    real, allocatable, dimension(:, :, :) :: prev_p, prev_vx, prev_vy, prev_vz
    real, allocatable, dimension(:, :, :) :: vp, rho, bk

    real, allocatable, dimension(:, :, :) :: grad_vp, grad_rho

    real, allocatable, dimension(:, :, :) :: src_p, rec_p, rec_pr
    real, allocatable, dimension(:, :, :) :: src_vx, src_vy, src_vz
    real, allocatable, dimension(:, :, :) :: rec_vx, rec_vy, rec_vz

    real, allocatable, dimension(:, :, :) :: src_hilbert, rec_hilbert

    real, allocatable, dimension(:, :, :) :: energy_src_v, energy_rec_v
    real, allocatable, dimension(:, :, :) :: energy_src_a, energy_rec_a

    real, allocatable, dimension(:) :: snaps

    integer :: np
    integer :: nomega, iw
    complex :: cmpcoef

    real :: idx, idy, idz

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

    integer :: nx1, nx2, ny1, ny2, nz1, nz2
    integer :: nx1_interior, nx2_interior, ny1_interior, ny2_interior, nz1_interior, nz2_interior

    real :: pmlvp

    real, allocatable, dimension(:, :, :) :: snapp

    integer, parameter :: htlen = 30

contains

    subroutine prepare_modeling(this)

        type(wave_solver_acoustic_iso_3d), intent(in) :: this

        real :: minv, maxv, dtstable, f0clean, temp

        nx = this%nx
        ny = this%ny
        nz = this%nz
        dx = this%dx
        dy = this%dy
        dz = this%dz
        ox = this%ox
        oy = this%oy
        oz = this%oz

        idx = 1.0d0/dx
        idy = 1.0d0/dy
        idz = 1.0d0/dz

        pml = this%pml

        ! Domain decomposition
        if (yn_free_surface) then
            call domain_decomp_regular_group(nx + 2*pml, ny + 2*pml, nz + pml, nx1, nx2, ny1, ny2, nz1, nz2, &
                weights1=[ones(pml)*1.25, ones(nx), ones(pml)*1.25], &
                weights2=[ones(pml)*1.25, ones(ny), ones(pml)*1.25], &
                weights3=[ones(nz), ones(pml)*1.25])
        else
            call domain_decomp_regular_group(nx + 2*pml, ny + 2*pml, nz + 2*pml, nx1, nx2, ny1, ny2, nz1, nz2, &
                weights1=[ones(pml)*1.25, ones(nx), ones(pml)*1.25], &
                weights2=[ones(pml)*1.25, ones(ny), ones(pml)*1.25], &
                weights3=[ones(pml)*1.25, ones(nz), ones(pml)*1.25])
        end if
        nx1 = nx1 - pml
        nx2 = nx2 - pml
        ny1 = ny1 - pml
        ny2 = ny2 - pml
        if (.not. yn_free_surface) then
            nz1 = nz1 - pml
            nz2 = nz2 - pml
        end if
        nx1_interior = max(1, nx1)
        nx2_interior = min(nx2, nx)
        ny1_interior = max(1, ny1)
        ny2_interior = min(ny2, ny)
        nz1_interior = max(1, nz1)
        nz2_interior = min(nz2, nz)

        ! Medium property models
        vp = permute(this%vp, 321)
        rho = permute(this%rho, 321)
        bk = rho*vp**2

        pmlvp = max(maxval(vp(1, :, :)), maxval(vp(nx, :, :)), &
            maxval(vp(:, 1, :)), maxval(vp(:, ny, :)), &
            maxval(vp(:, :, 1)), maxval(vp(:, :, nz)))

        call pad_array(vp, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(rho, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])
        call pad_array(bk, [pml + 1, pml + 1, pml + 1, pml + 1, pml + 1, pml + 1])

        ! Check stability and dispersion
        dt = this%dt
        tmax = this%tmax

        temp = sum(abs(fdcoefs))
        minv = minval(vp)
        maxv = maxval(vp)

        dtstable = 1.0/(temp*maxv*sqrt(1.0/dx**2 + 1.0/dy**2 + 1.0/dz**2))
        f0clean = minv/max(dx, dy, dz)/7.0

        if (rankid_group == 0) then
            call warn(date_time_compact()//' Stable dt = '//num2str(dtstable, '(es)') &
                //' s, clean f0 = '//num2str(f0clean, '(es)')//' Hz')
        end if

        call check_dt_f0(dt, dtstable, maxval(this%gmtr%srcr(:)%f0), f0clean)

        nt = nint(tmax/dt + 1)

        ! Crop models
        call alloc_array(vp, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=vp(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(rho, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=rho(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))
        call alloc_array(bk, [nx1, nx2, ny1, ny2, nz1, nz2], pad=1, &
            source=bk(nx1 - 1:nx2 + 1, ny1 - 1:ny2 + 1, nz1 - 1:nz2 + 1))

        ! Prepare geometry
        sgmtr = this%gmtr
        sgmtr%nx = nx
        sgmtr%ny = ny
        sgmtr%nz = nz
        sgmtr%dx = dx
        sgmtr%dy = dy
        sgmtr%dz = dz
        sgmtr%ox = ox
        sgmtr%oy = oy
        sgmtr%oz = oz
        sgmtr%xmin = max(ox, this%gmtr%xmin)
        sgmtr%xmax = min((nx - 1)*dx + ox, this%gmtr%xmax)
        sgmtr%ymin = max(oy, this%gmtr%ymin)
        sgmtr%ymax = min((ny - 1)*dy + oy, this%gmtr%ymax)
        sgmtr%zmin = max(oz, this%gmtr%zmin)
        sgmtr%zmax = min((nz - 1)*dz + oz, this%gmtr%zmax)
        sgmtr%sxmin = max(ox, this%gmtr%sxmin)
        sgmtr%sxmax = min((nx - 1)*dx + ox, this%gmtr%sxmax)
        sgmtr%symin = max(oy, this%gmtr%symin)
        sgmtr%symax = min((ny - 1)*dy + oy, this%gmtr%symax)
        sgmtr%szmin = max(oz, this%gmtr%szmin)
        sgmtr%szmax = min((nz - 1)*dz + oz, this%gmtr%szmax)
        sgmtr%rxmin = max(ox, this%gmtr%rxmin)
        sgmtr%rxmax = min((nx - 1)*dx + ox, this%gmtr%rxmax)
        sgmtr%rymin = max(oy, this%gmtr%rymin)
        sgmtr%rymax = min((ny - 1)*dy + oy, this%gmtr%rymax)
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

        call mpibarrier_group

    end subroutine prepare_modeling

    !
    !> Source wavefields
    !
    subroutine alloc_forward_wavefield

        call alloc_array(vx, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vy, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vz, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(p, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdxvx_xmin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdxvx_xmax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdyvy_ymin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdyvy_ymax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdzvz_zmin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdzvz_zmax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdxp_xmin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdxp_xmax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdyp_ymin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdyp_ymax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdzp_zmin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdzp_zmax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

    end subroutine alloc_forward_wavefield

    !
    !> Receiver wavefields
    !
    subroutine alloc_adjoint_wavefield

        call alloc_array(vxr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vyr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(vzr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(pr, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdxvxr_xmin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdxvxr_xmax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdyvyr_ymin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdyvyr_ymax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdzvzr_zmin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdzvzr_zmax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdxpr_xmin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdxpr_xmax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdypr_ymin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdypr_ymax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdzpr_zmin, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)
        call alloc_array(memory_pdzpr_zmax, [nx1, nx2, ny1, ny2, nz1, nz2], pad=fdhalf)

    end subroutine alloc_adjoint_wavefield

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

end module
